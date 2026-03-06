from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import StepResponse, StepUpdateRequest, StepUpdateResponse

from .connection import get_connection
from .output import json_print

step_app = typer.Typer()


@step_app.command("get")
def get_step(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Get the current step."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}/step")
    json_print(StepResponse.model_validate(response.json()))


@step_app.command("set")
def set_step(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    index: Annotated[int, typer.Argument(help="Step index")],
) -> None:
    """Set the current step."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    request = StepUpdateRequest(step=index)
    response = conn.put(f"/v1/rooms/{room}/step", json=request.model_dump())
    json_print(StepUpdateResponse.model_validate(response.json()))
