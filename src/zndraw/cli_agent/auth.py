from __future__ import annotations

from typing import Annotated

import httpx
import typer

from .connection import resolve_url
from .output import json_print

auth_app = typer.Typer()


@auth_app.command("login")
def login(
    ctx: typer.Context,
    user: Annotated[str | None, typer.Option(help="Username")] = None,
    password: Annotated[str | None, typer.Option(help="Password")] = None,
) -> None:
    """Login with credentials or as a guest."""
    url = resolve_url(ctx.obj["url"])

    if user is not None and password is not None:
        response = httpx.post(
            f"{url}/v1/auth/jwt/login",
            data={"username": user, "password": password},
        )
    else:
        response = httpx.post(f"{url}/v1/auth/guest")

    response.raise_for_status()
    json_print(response.json())
