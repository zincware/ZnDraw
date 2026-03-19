from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import (
    BookmarkCreateRequest,
    BookmarkResponse,
    BookmarksResponse,
    StatusResponse,
)

from .connection import (
    PasswordOpt,
    RoomOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

bookmarks_app = typer.Typer()


@bookmarks_app.command("list")
def list_bookmarks(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
) -> None:
    """List all bookmarks."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/bookmarks", headers=vis.api.get_headers()
        )
        vis.api.raise_for_status(resp)
        json_print(BookmarksResponse.model_validate(resp.json()))
        vis.disconnect()


@bookmarks_app.command("set")
def set_bookmark(
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
    label: Annotated[
        str | None, typer.Argument(help="Bookmark label (auto-generated if omitted)")
    ] = None,
) -> None:
    """Set a bookmark."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        if index is None:
            index = vis.step
        if label is None:
            label = f"Frame {index}"
        request = BookmarkCreateRequest(label=label)
        data = vis.api.set_bookmark(index, request.label)
        json_print(BookmarkResponse.model_validate(data))
        vis.disconnect()


@bookmarks_app.command("delete")
def delete_bookmark(
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
) -> None:
    """Delete a bookmark."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        if index is None:
            index = vis.step
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/bookmarks/{index}", headers=vis.api.get_headers()
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
