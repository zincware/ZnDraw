from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import (
    BookmarkCreateRequest,
    BookmarkResponse,
    BookmarksResponse,
    StatusResponse,
)

from .connection import cli_error_handler, get_zndraw
from .output import json_print

bookmarks_app = typer.Typer()


@bookmarks_app.command("list")
def list_bookmarks(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List all bookmarks."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/bookmarks", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(BookmarksResponse.model_validate(resp.json()))
        vis.disconnect()


@bookmarks_app.command("set")
def set_bookmark(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
    label: Annotated[str, typer.Option("--label", help="Bookmark label")] = "",
) -> None:
    """Set a bookmark."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        if index is None:
            index = vis.step
        if not label:
            label = f"Frame {index}"
        request = BookmarkCreateRequest(label=label)
        data = vis.api.set_bookmark(index, request.label)
        json_print(BookmarkResponse.model_validate(data))
        vis.disconnect()


@bookmarks_app.command("delete")
def delete_bookmark(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
) -> None:
    """Delete a bookmark."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        if index is None:
            index = vis.step
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/bookmarks/{index}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
