from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import PresenceResponse

from .connection import cli_error_handler, get_zndraw
from .output import json_print

sessions_app = typer.Typer()


@sessions_app.command("list")
def list_sessions(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List all sessions (presence) in a room."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        sids = vis.api.list_sessions()
        json_print(PresenceResponse(items=sids))
        vis.disconnect()
