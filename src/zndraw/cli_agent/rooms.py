from __future__ import annotations

import uuid
import webbrowser
from typing import Annotated

import typer

from zndraw.schemas import (
    CollectionResponse,
    RoomCreate,
    RoomCreateResponse,
    RoomPatchResponse,
    RoomResponse,
)

from .connection import (
    PasswordOpt,
    RoomOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_connection,
    get_zndraw,
    resolve_room,
)
from .output import json_print

rooms_app = typer.Typer()


@rooms_app.command("list")
def list_rooms(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    search: Annotated[str | None, typer.Option(help="Search query")] = None,
) -> None:
    """List all rooms."""
    with cli_error_handler():
        conn = get_connection(url, token, user, password)
        try:
            params: dict[str, str] = {}
            if search is not None:
                params["search"] = search
            resp = conn.get("/v1/rooms", params=params)
            json_print(CollectionResponse[RoomResponse].model_validate(resp.json()))
        finally:
            conn.close()


@rooms_app.command("create")
def create_room(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
    copy_from: Annotated[
        str | None, typer.Option("--copy-from", help="Copy from existing room ID")
    ] = None,
) -> None:
    """Create a new room."""
    with cli_error_handler():
        conn = get_connection(url, token, user, password)
        request = RoomCreate(
            room_id=room if room is not None else str(uuid.uuid4()),
            copy_from=copy_from,
        )
        response = conn.post("/v1/rooms", json=request.model_dump())
        json_print(RoomCreateResponse.model_validate(response.json()))


@rooms_app.command("info")
def room_info(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get room info."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        json_print(RoomResponse.model_validate(vis.api.get_room_info()))
        vis.disconnect()


@rooms_app.command("lock")
def lock_room(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
) -> None:
    """Lock a room (prevent edits by non-admin users)."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        vis.locked = True
        json_print(RoomPatchResponse.model_validate(vis.api.get_room_info()))
        vis.disconnect()


@rooms_app.command("unlock")
def unlock_room(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
) -> None:
    """Unlock a room (allow edits again)."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        vis.locked = False
        json_print(RoomPatchResponse.model_validate(vis.api.get_room_info()))
        vis.disconnect()


@rooms_app.command("open")
def open_room(
    url: UrlOpt = None,
    _token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Open a room in the browser."""
    with cli_error_handler():
        from zndraw.client.settings import ClientSettings

        from .connection import EXIT_CLIENT_ERROR, EXIT_CONNECTION_ERROR, die

        room = resolve_room(room)
        overrides = {"url": url} if url is not None else {}
        try:
            settings = ClientSettings(**overrides)
        except Exception as exc:
            die("Configuration Error", str(exc), 400, EXIT_CLIENT_ERROR)

        if settings.url is None:
            die(
                "No Server Found",
                "No running zndraw server found. Start one with `uv run zndraw` or pass `--url`.",
                503,
                EXIT_CONNECTION_ERROR,
            )
        room_url = f"{settings.url}/rooms/{room}"
        typer.echo(room_url)
        webbrowser.open(room_url)


@rooms_app.command("set-default")
def set_default_room(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
) -> None:
    """Set a room as the default template for new rooms."""
    with cli_error_handler():
        room = resolve_room(room)
        conn = get_connection(url, token, user, password)
        response = conn.put("/v1/server-settings/default-room", json={"room_id": room})
        json_print(response.json())
