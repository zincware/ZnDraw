from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import (
    BookmarkCreateRequest,
    BookmarkResponse,
    BookmarksResponse,
    StatusResponse,
)

from .connection import get_connection, get_current_step
from .output import json_print

bookmarks_app = typer.Typer()


@bookmarks_app.command("list")
def list_bookmarks(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List all bookmarks."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}/bookmarks")
    json_print(BookmarksResponse.model_validate(response.json()))


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
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    if index is None:
        index = get_current_step(conn, room)
    if not label:
        label = f"Frame {index}"
    request = BookmarkCreateRequest(label=label)
    response = conn.put(
        f"/v1/rooms/{room}/bookmarks/{index}", json=request.model_dump()
    )
    json_print(BookmarkResponse.model_validate(response.json()))


@bookmarks_app.command("delete")
def delete_bookmark(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
) -> None:
    """Delete a bookmark."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    if index is None:
        index = get_current_step(conn, room)
    response = conn.delete(f"/v1/rooms/{room}/bookmarks/{index}")
    json_print(StatusResponse.model_validate(response.json()))
