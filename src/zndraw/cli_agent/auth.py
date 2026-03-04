"""CLI auth commands: login, status, logout."""

from __future__ import annotations

import time
import webbrowser
from datetime import UTC, datetime

import httpx
import typer

from zndraw.server_manager import TokenEntry, TokenStore

from .connection import TokenOpt, UrlOpt, cli_error_handler, resolve_url
from .output import json_print

auth_app = typer.Typer()


def _get_token_store() -> TokenStore:
    """Return the default TokenStore (testable seam)."""
    return TokenStore()


@auth_app.command("login")
def login(
    url: UrlOpt = None,
    code: bool = typer.Option(
        False,
        "--code",
        help="Print URL only, don't open browser",  # noqa: FBT003
    ),
) -> None:
    """Login via browser approval (device-code flow)."""
    with cli_error_handler():
        resolved_url = resolve_url(url)
        store = _get_token_store()

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
) -> None:
    """Show current authentication identity."""
    with cli_error_handler():
        resolved_url = resolve_url(url)
        store = _get_token_store()

        # Determine token and source — no guest fallback
        if token is not None:
            active_token = token
            token_source = "flag"
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
                store.delete(resolved_url)
                json_print(
                    {
                        "server": resolved_url,
                        "user_id": None,
                        "email": None,
                        "is_superuser": False,
                        "token_source": "expired",
                    }
                )
                return
            user = resp.json()

        json_print(
            {
                "server": resolved_url,
                "user_id": user.get("id"),
                "email": user.get("email"),
                "is_superuser": user.get("is_superuser", False),
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
        store = _get_token_store()
        store.delete(resolved_url)
        typer.echo(f"Logged out from {resolved_url}")
