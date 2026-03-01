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

from .connection import cli_error_handler, get_zndraw
from .output import json_print

extensions_app = typer.Typer(name="extensions", help="Extension operations")


@extensions_app.command("list")
def list_extensions(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List available extensions for a room."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.list_extensions()
        json_print(PaginatedResponse[JobSummary].model_validate(data))
        vis.disconnect()


@extensions_app.command("describe")
def describe(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    name: Annotated[str, typer.Argument(help="Fully qualified extension name")],
) -> None:
    """Describe an extension by its fully qualified name."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.get_extension(name)
        json_print(JobResponse.model_validate(data))
        vis.disconnect()


@extensions_app.command(
    "run",
    context_settings={"allow_extra_args": True, "allow_interspersed_args": False},
)
def run_extension(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    extension_name: Annotated[
        str, typer.Argument(help="Fully qualified extension name")
    ],
) -> None:
    """Run an extension with optional extra parameters passed as --key value pairs."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)

        extra_args = ctx.args
        payload: dict = {}
        it = iter(extra_args)
        for token in it:
            key = token.lstrip("-")
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
