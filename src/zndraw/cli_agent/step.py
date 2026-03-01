from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import StepResponse, StepUpdateRequest, StepUpdateResponse

from .connection import cli_error_handler, get_zndraw
from .output import json_print

step_app = typer.Typer()


@step_app.command("get")
def get_step(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Get the current step."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.get_step()
        json_print(StepResponse.model_validate(data))
        vis.disconnect()


@step_app.command("set")
def set_step(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    index: Annotated[int, typer.Argument(help="Step index")],
) -> None:
    """Set the current step."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.update_step(index)
        json_print(StepUpdateResponse.model_validate(data))
        vis.disconnect()
