from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import PresenceResponse, SessionSettingsResponse

from .connection import get_connection
from .output import json_print

sessions_app = typer.Typer()


@sessions_app.command("list")
def list_sessions(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List all sessions (presence) in a room."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}/presence")
    json_print(PresenceResponse.model_validate(response.json()))


@sessions_app.command("settings")
def session_settings(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    sid: Annotated[str, typer.Argument(help="Session ID")],
) -> None:
    """Get settings for a session."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}/sessions/{sid}/settings")
    json_print(SessionSettingsResponse.model_validate(response.json()))
