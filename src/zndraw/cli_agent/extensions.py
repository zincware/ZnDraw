from __future__ import annotations

import json
from typing import Annotated

import typer
from zndraw_joblib.schemas import (
    JobResponse,
    JobSummary,
    PaginatedResponse,
)

from .connection import (
    PasswordOpt,
    RoomOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

extensions_app = typer.Typer(name="extensions", help="Extension operations")


@extensions_app.command("list")
def list_extensions(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
) -> None:
    """List available extensions for a room."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        data = vis.api.list_extensions()
        json_print(PaginatedResponse[JobSummary].model_validate(data))
        vis.disconnect()


@extensions_app.command("describe")
def describe(
    name: Annotated[
        str | None, typer.Argument(help="Fully qualified extension name")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
) -> None:
    """Describe an extension by its fully qualified name."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Extension name is required")
        vis = get_zndraw(url, token, room, user, password)
        data = vis.api.get_extension(name)
        json_print(JobResponse.model_validate(data))
        vis.disconnect()


@extensions_app.command(
    "run",
    context_settings={
        "allow_extra_args": True,
        "allow_interspersed_args": True,
        "ignore_unknown_options": True,
    },
)
def run_extension(
    ctx: typer.Context,
    extension_name: Annotated[
        str | None, typer.Argument(help="Fully qualified extension name")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    wait: Annotated[
        bool, typer.Option("--wait", help="Wait for the task to complete")
    ] = False,
    timeout: Annotated[
        float, typer.Option("--timeout", help="Timeout in seconds (with --wait)")
    ] = 300,
) -> None:
    """Run an extension with optional extra parameters passed as --key value pairs."""
    with cli_error_handler():
        room = resolve_room(room)
        if extension_name is None:
            raise typer.BadParameter("Extension name is required")
        vis = get_zndraw(url, token, room, user, password)

        extra_args = ctx.args
        payload: dict = {}
        it = iter(extra_args)
        for arg_token in it:
            key = arg_token.lstrip("-")
            try:
                raw_value = next(it)
            except StopIteration:
                raise typer.BadParameter(
                    f"Missing value for argument --{key}"
                ) from None
            try:
                value = json.loads(raw_value)
            except json.JSONDecodeError:
                value = raw_value
            payload[key] = value

        task = vis.run(extension_name, **payload)
        if wait:
            task.wait(timeout=timeout)
        json_print(task.fetch())
        vis.disconnect()
