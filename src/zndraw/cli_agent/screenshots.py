from __future__ import annotations

import base64
import pathlib
from typing import Annotated

import typer

from zndraw.schemas import (
    OffsetPage,
    ScreenshotListItem,
    ScreenshotResponse,
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

screenshots_app = typer.Typer(name="screenshots", help="Screenshot operations")


@screenshots_app.command("list")
def list_screenshots(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """List screenshots for a room."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        items = vis.api.list_screenshots()
        json_print(
            OffsetPage[ScreenshotListItem].model_validate(
                {"items": items, "total": len(items), "limit": len(items), "offset": 0}
            )
        )
        vis.disconnect()


@screenshots_app.command("request")
def request_screenshot(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    session: Annotated[
        str | None, typer.Option("--session", help="Target session SID")
    ] = None,
) -> None:
    """Request a screenshot from the current user's first active session."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)

        if session is not None:
            target_sid = session
        else:
            # Get current user's email, then find their session
            me = vis.api.get_me()
            my_email = me.get("email", "")
            all_sessions = vis.api.list_sessions()
            own_sessions = [s for s in all_sessions if s.email == my_email]
            if not own_sessions:
                typer.echo(
                    "No active browser sessions for your user in this room.",
                    err=True,
                )
                raise typer.Exit(code=1)
            target_sid = own_sessions[0].sid

        result = vis.api.create_screenshot_capture(target_sid)
        json_print(ScreenshotResponse.model_validate(result))
        vis.disconnect()


@screenshots_app.command("get")
def get(
    id: Annotated[str | None, typer.Argument(help="Screenshot ID")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    output: Annotated[
        str | None, typer.Option(help="Path to save screenshot file")
    ] = None,
) -> None:
    """Get a screenshot by ID."""
    with cli_error_handler():
        room = resolve_room(room)
        if id is None:
            raise typer.BadParameter("Screenshot ID is required")
        vis = get_zndraw(url, token, room)
        data = vis.api.get_screenshot(int(id))
        resp = ScreenshotResponse.model_validate(data)

        if output:
            if resp.data is None:
                typer.echo("Screenshot has no image data.", err=True)
                raise typer.Exit(code=1)
            image_data = base64.b64decode(resp.data)
            pathlib.Path(output).write_bytes(image_data)
        else:
            json_print(resp)
        vis.disconnect()
