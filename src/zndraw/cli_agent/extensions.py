from __future__ import annotations

import json
from typing import Annotated

import typer
from zndraw_joblib.schemas import (
    JobResponse,
    JobSummary,
    PaginatedResponse,
    TaskResponse,
    TaskSubmitRequest,
)

from .connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
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
) -> None:
    """List available extensions for a room."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
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
) -> None:
    """Describe an extension by its fully qualified name."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Extension name is required")
        vis = get_zndraw(url, token, room)
        data = vis.api.get_extension(name)
        json_print(JobResponse.model_validate(data))
        vis.disconnect()


@extensions_app.command(
    "run",
    context_settings={"allow_extra_args": True, "allow_interspersed_args": False},
)
def run_extension(
    ctx: typer.Context,
    extension_name: Annotated[
        str | None, typer.Argument(help="Fully qualified extension name")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Run an extension with optional extra parameters passed as --key value pairs."""
    with cli_error_handler():
        room = resolve_room(room)
        if extension_name is None:
            raise typer.BadParameter("Extension name is required")
        vis = get_zndraw(url, token, room)

        extra_args = ctx.args
        payload: dict = {}
        it = iter(extra_args)
        for arg_token in it:
            key = arg_token.lstrip("-")
            try:
                raw_value = next(it)
            except StopIteration:
                raise typer.BadParameter(f"Missing value for argument --{key}")
            try:
                value = json.loads(raw_value)
            except json.JSONDecodeError:
                value = raw_value
            payload[key] = value

        data = vis.api.submit_task(extension_name, payload)
        json_print(TaskResponse.model_validate(data))
        vis.disconnect()
