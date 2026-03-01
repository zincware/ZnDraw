from __future__ import annotations

from typing import Annotated

import typer

from zndraw import ZnDraw

from .connection import cli_error_handler, resolve_url
from .output import json_print

auth_app = typer.Typer()


@auth_app.command("login")
def login(
    ctx: typer.Context,
    user: Annotated[str | None, typer.Option(help="Username")] = None,
    password: Annotated[str | None, typer.Option(help="Password")] = None,
) -> None:
    """Login with credentials or as a guest."""
    with cli_error_handler():
        url = resolve_url(ctx.obj["url"])
        if user is not None and password is not None:
            token = ZnDraw.login(url=url, username=user, password=password)
        else:
            # Guest auth
            import httpx

            response = httpx.post(f"{url}/v1/auth/guest")
            response.raise_for_status()
            token = response.json()["access_token"]
        json_print({"access_token": token})
