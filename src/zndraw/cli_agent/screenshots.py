from __future__ import annotations

import base64
import pathlib
from typing import Annotated

import typer

from zndraw.schemas import (
    OffsetPage,
    ScreenshotCaptureCreate,
    ScreenshotListItem,
    ScreenshotResponse,
    SessionsListResponse,
)

from .connection import get_connection
from .output import json_print

screenshots_app = typer.Typer(name="screenshots", help="Screenshot operations")


@screenshots_app.command("list")
def list_screenshots(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List screenshots for a room."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/rooms/{room}/screenshots")
    json_print(OffsetPage[ScreenshotListItem].model_validate(resp.json()))


@screenshots_app.command("request")
def request_screenshot(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Request a screenshot from the first active session in the room."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])

    sessions_resp = conn.get(f"/v1/rooms/{room}/sessions")
    sessions = SessionsListResponse.model_validate(sessions_resp.json())
    if not sessions.items:
        typer.echo("No active sessions available for this room.", err=True)
        raise typer.Exit(code=1)

    session_id = sessions.items[0]
    request = ScreenshotCaptureCreate(session_id=session_id)
    resp = conn.post(f"/v1/rooms/{room}/screenshots", json=request.model_dump())
    json_print(ScreenshotResponse.model_validate(resp.json()))


@screenshots_app.command("get")
def get(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    id: Annotated[str, typer.Argument(help="Screenshot ID")],
    output: Annotated[
        str | None, typer.Option(help="Path to save screenshot file")
    ] = None,
) -> None:
    """Get a screenshot by ID."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/rooms/{room}/screenshots/{id}")
    data = ScreenshotResponse.model_validate(resp.json())

    if output:
        if data.data is None:
            typer.echo("Screenshot has no image data.", err=True)
            raise typer.Exit(code=1)
        image_data = base64.b64decode(data.data)
        pathlib.Path(output).write_bytes(image_data)
    else:
        json_print(data)
