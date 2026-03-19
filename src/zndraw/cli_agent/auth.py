"""CLI auth commands: login, status, logout."""

from __future__ import annotations

import time
import webbrowser
from datetime import UTC, datetime

import httpx
import typer

from zndraw.server_manager import TokenEntry

from .connection import (
    PasswordOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_token_store,
    resolve_url,
)
from .output import json_print

auth_app = typer.Typer()


@auth_app.command("login")
def login(
    url: UrlOpt = None,
    code: bool = typer.Option(
        False,
        "--code",
        help="Print URL only, don't open browser",
    ),
) -> None:
    """Login via browser approval (device-code flow)."""
    with cli_error_handler():
        resolved_url = resolve_url(url)
        store = get_token_store()

        with httpx.Client(base_url=resolved_url, timeout=30.0) as client:
            # 1. Create challenge
            resp = client.post("/v1/auth/cli-login")
            resp.raise_for_status()
            challenge = resp.json()

            code_str = challenge["code"]
            secret = challenge["secret"]
            approve_url = f"{resolved_url}/auth/cli?code={code_str}"

            typer.echo(f"\n  Your code: {code_str}\n")

            if code:
                typer.echo(f"  Visit: {approve_url}")
            else:
                typer.echo(f"  Opening browser... (or visit: {approve_url})")
                webbrowser.open(approve_url)

            typer.echo("  Waiting for approval...\n")

            # 2. Poll for approval
            for _ in range(300):  # 5 min max (1s intervals)
                time.sleep(1)
                poll = client.get(
                    f"/v1/auth/cli-login/{code_str}",
                    params={"secret": secret},
                )

                if poll.status_code == 404:
                    typer.echo("Login rejected.", err=True)
                    raise typer.Exit(code=1)

                if poll.status_code == 410:
                    typer.echo("Login challenge expired.", err=True)
                    raise typer.Exit(code=1)

                data = poll.json()
                if data["status"] == "approved" and data["token"]:
                    token = data["token"]

                    # Fetch user info for storage
                    me_resp = client.get(
                        "/v1/auth/users/me",
                        headers={"Authorization": f"Bearer {token}"},
                    )
                    email = me_resp.json().get("email", "unknown")

                    store.set(
                        resolved_url,
                        TokenEntry(
                            access_token=token,
                            email=email,
                            stored_at=datetime.now(UTC),
                        ),
                    )
                    typer.echo(f"Logged in as {email}")
                    typer.echo(f"Token saved to {store.path}")
                    return

            typer.echo("Login timed out.", err=True)
            raise typer.Exit(code=1)


@auth_app.command("status")
def status(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
) -> None:
    """Show current authentication identity."""
    with cli_error_handler():
        resolved_url = resolve_url(url)
        store = get_token_store()

        # Validate mutual exclusion
        if token is not None and (user is not None or password is not None):
            msg = "Cannot combine --token with --user/--password"
            raise ValueError(msg)
        if user is not None and password is None:
            msg = "Missing --password (required when --user is provided)"
            raise ValueError(msg)
        if password is not None and user is None:
            msg = "Missing --user (required when --password is provided)"
            raise ValueError(msg)

        # Determine token and source — no guest fallback
        if token is not None:
            active_token = token
            token_source = "flag"
        elif user is not None and password is not None:
            with httpx.Client(base_url=resolved_url, timeout=10.0) as client:
                resp = client.post(
                    "/v1/auth/jwt/login",
                    data={"username": user, "password": password},
                )
                resp.raise_for_status()
                active_token = resp.json()["access_token"]
            token_source = "login"
        else:
            entry = store.get(resolved_url)
            if entry is not None:
                active_token = entry.access_token
                token_source = "stored"
            else:
                json_print(
                    {
                        "server": resolved_url,
                        "user_id": None,
                        "email": None,
                        "is_superuser": False,
                        "token_source": "none",
                    }
                )
                return

        with httpx.Client(
            base_url=resolved_url,
            headers={"Authorization": f"Bearer {active_token}"},
            timeout=10.0,
        ) as client:
            resp = client.get("/v1/auth/users/me")
            if resp.status_code != 200:
                # Only delete stored token on auth failures, not server errors
                if resp.status_code in (401, 403) and token_source == "stored":
                    store.delete(resolved_url)
                json_print(
                    {
                        "server": resolved_url,
                        "user_id": None,
                        "email": None,
                        "is_superuser": False,
                        "token_source": "expired"
                        if resp.status_code in (401, 403)
                        else "error",
                    }
                )
                return
            user_data = resp.json()

        json_print(
            {
                "server": resolved_url,
                "user_id": user_data.get("id"),
                "email": user_data.get("email"),
                "is_superuser": user_data.get("is_superuser", False),
                "token_source": token_source,
            }
        )


@auth_app.command("logout")
def logout(
    url: UrlOpt = None,
) -> None:
    """Remove stored token for the current server."""
    with cli_error_handler():
        resolved_url = resolve_url(url)
        store = get_token_store()
        store.delete(resolved_url)
        typer.echo(f"Logged out from {resolved_url}")
