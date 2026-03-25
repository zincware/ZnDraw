"""CLI admin commands: users list, users login."""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Annotated

import typer

from zndraw.state_file import StateFile, TokenEntry

from .connection import (
    PasswordOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_connection,
)
from .output import json_print

admin_app = typer.Typer(name="admin", help="Admin operations (superuser only)")
users_app = typer.Typer(name="users", help="User management")
admin_app.add_typer(users_app, name="users")


@users_app.command("list")
def list_users(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    limit: Annotated[int, typer.Option(help="Max results")] = 100,
    offset: Annotated[int, typer.Option(help="Pagination offset")] = 0,
) -> None:
    """List all users (superuser only)."""
    with cli_error_handler():
        conn = get_connection(url, token, user, password)
        try:
            resp = conn.get(
                "/v1/admin/users", params={"limit": limit, "offset": offset}
            )
            json_print(resp.json())
        finally:
            conn.close()


@users_app.command("login")
def login_as_user(
    user_id: Annotated[str, typer.Argument(help="User ID to mint token for")],
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
) -> None:
    """Mint a token as another user (superuser only)."""
    with cli_error_handler():
        conn = get_connection(url, token, user, password)
        state = StateFile()

        try:
            # Mint token for target user
            resp = conn.post(f"/v1/admin/users/{user_id}/token")
            data = resp.json()
            new_token = data["access_token"]

            # Fetch target user's email
            me_resp = conn.client.get(
                "/v1/auth/users/me",
                headers={"Authorization": f"Bearer {new_token}"},
            )
            email = me_resp.json().get("email", "unknown")
        finally:
            conn.close()

        state.add_token(
            conn.base_url,
            TokenEntry(
                access_token=new_token,
                email=email,
                stored_at=datetime.now(UTC),
            ),
        )
        typer.echo(f"Logged in as {email}")
        typer.echo(f"Token saved to {state.path}")
