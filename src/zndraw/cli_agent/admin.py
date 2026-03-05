"""CLI admin commands: users list, users login."""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Annotated

import httpx
import typer

from zndraw.server_manager import TokenEntry

from .connection import (
    TokenOpt,
    UrlOpt,
    get_token_store,
    cli_error_handler,
    resolve_token,
    resolve_url,
)
from .output import json_print

admin_app = typer.Typer(name="admin", help="Admin operations (superuser only)")
users_app = typer.Typer(name="users", help="User management")
admin_app.add_typer(users_app, name="users")


@users_app.command("list")
def list_users(
    url: UrlOpt = None,
    token: TokenOpt = None,
    limit: Annotated[int, typer.Option(help="Max results")] = 100,
    offset: Annotated[int, typer.Option(help="Pagination offset")] = 0,
) -> None:
    """List all users (superuser only)."""
    with cli_error_handler():
        resolved_url = resolve_url(url)
        resolved_token = resolve_token(resolved_url, token)

        with httpx.Client(
            base_url=resolved_url,
            headers={"Authorization": f"Bearer {resolved_token}"},
            timeout=10.0,
        ) as client:
            resp = client.get(
                "/v1/admin/users", params={"limit": limit, "offset": offset}
            )
            resp.raise_for_status()
            json_print(resp.json())


@users_app.command("login")
def login_as_user(
    user_id: Annotated[str, typer.Argument(help="User ID to mint token for")],
    url: UrlOpt = None,
    token: TokenOpt = None,
) -> None:
    """Mint a token as another user (superuser only)."""
    with cli_error_handler():
        resolved_url = resolve_url(url)
        resolved_token = resolve_token(resolved_url, token)
        store = get_token_store()

        with httpx.Client(
            base_url=resolved_url,
            headers={"Authorization": f"Bearer {resolved_token}"},
            timeout=10.0,
        ) as client:
            # Mint token for target user
            resp = client.post(f"/v1/admin/users/{user_id}/token")
            resp.raise_for_status()
            data = resp.json()
            new_token = data["access_token"]

            # Fetch target user's email
            me_resp = client.get(
                "/v1/auth/users/me",
                headers={"Authorization": f"Bearer {new_token}"},
            )
            email = me_resp.json().get("email", "unknown")

        store.set(
            resolved_url,
            TokenEntry(
                access_token=new_token,
                email=email,
                stored_at=datetime.now(UTC),
            ),
        )
        typer.echo(f"Logged in as {email}")
        typer.echo(f"Token saved to {store.path}")
