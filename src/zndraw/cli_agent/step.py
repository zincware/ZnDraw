from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import StepResponse, StepUpdateResponse

from .connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

step_app = typer.Typer()


@step_app.command("get")
def get_step(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get the current step."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        data = vis.api.get_step()
        json_print(StepResponse.model_validate(data))
        vis.disconnect()


@step_app.command("set")
def set_step(
    index: Annotated[int | None, typer.Argument(help="Step index")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Set the current step."""
    with cli_error_handler():
        room = resolve_room(room)
        if index is None:
            raise typer.BadParameter("Step index is required")
        vis = get_zndraw(url, token, room)
        data = vis.api.update_step(index)
        json_print(StepUpdateResponse.model_validate(data))
        vis.disconnect()
