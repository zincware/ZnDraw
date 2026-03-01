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

from .connection import cli_error_handler, get_zndraw
from .output import json_print

screenshots_app = typer.Typer(name="screenshots", help="Screenshot operations")


@screenshots_app.command("list")
def list_screenshots(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List screenshots for a room."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        items = vis.api.list_screenshots()
        json_print(
            OffsetPage[ScreenshotListItem].model_validate(
                {"items": items, "total": len(items), "limit": len(items), "offset": 0}
            )
        )
        vis.disconnect()


@screenshots_app.command("request")
def request_screenshot(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Request a screenshot from the first active session in the room."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        sids = vis.api.list_sessions()
        if not sids:
            typer.echo("No active sessions available for this room.", err=True)
            raise typer.Exit(code=1)

        result = vis.api.create_screenshot_capture(sids[0])
        json_print(ScreenshotResponse.model_validate(result))
        vis.disconnect()


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
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
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
