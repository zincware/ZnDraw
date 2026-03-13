from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import ActiveCameraResponse, SessionsListResponse

from .connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

sessions_app = typer.Typer()


@sessions_app.command("list")
def list_sessions(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """List all active frontend sessions in a room."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        items = vis.api.list_sessions()
        json_print(SessionsListResponse(items=items))
        vis.disconnect()


@sessions_app.command("get-camera")
def get_camera(
    session_id: Annotated[str | None, typer.Argument(help="Session SID")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get the active camera key for a session."""
    with cli_error_handler():
        room = resolve_room(room)
        if session_id is None:
            raise typer.BadParameter("Session ID is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/sessions/{session_id}/active-camera",
            headers=vis.api.get_headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(ActiveCameraResponse.model_validate(resp.json()))
        vis.disconnect()


@sessions_app.command("set-camera")
def set_camera(
    session_id: Annotated[str | None, typer.Argument(help="Session SID")] = None,
    camera_key: Annotated[
        str | None, typer.Argument(help="Camera geometry key")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Set the active camera for a session."""
    with cli_error_handler():
        room = resolve_room(room)
        if session_id is None:
            raise typer.BadParameter("Session ID is required")
        if camera_key is None:
            raise typer.BadParameter("Camera key is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.put(
            f"/v1/rooms/{vis.room}/sessions/{session_id}/active-camera",
            json={"active_camera": camera_key},
            headers=vis.api.get_headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(ActiveCameraResponse.model_validate(resp.json()))
        vis.disconnect()
