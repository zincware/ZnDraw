from __future__ import annotations

import uuid
import webbrowser
from typing import Annotated
from uuid import UUID

import typer

from zndraw.schemas import (
    CollectionResponse,
    RoomCreate,
    RoomCreateResponse,
    RoomResponse,
)

from .connection import (
    EXIT_CLIENT_ERROR,
    EXIT_CONNECTION_ERROR,
    PasswordOpt,
    RoomOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    die,
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
    name: Annotated[
        str | None, typer.Option("--name", help="Room name (within your namespace)")
    ] = None,
    copy_from: Annotated[
        str | None, typer.Option("--copy-from", help="Copy from existing room")
    ] = None,
) -> None:
    """Create a new room in the caller's user namespace."""
    with cli_error_handler():
        conn = get_connection(url, token, user, password)
        me = conn.get("/v1/auth/users/me").json()
        request = RoomCreate(
            owner_id=UUID(me["id"]),
            name=name if name is not None else str(uuid.uuid4()),
        )
        if copy_from is not None:
            request = request.model_copy(update={"copy_from": copy_from})
        response = conn.post("/v1/rooms", json=request.model_dump(mode="json"))
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


@rooms_app.command("open")
def open_room(
    url: UrlOpt = None,
    _token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Open a room in the browser."""
    with cli_error_handler():
        from zndraw.client.settings import ClientSettings

        room = resolve_room(room)
        if "/" not in room:
            die(
                "Invalid room",
                "Room must be in '<owner_uuid>/<name>' form.",
                400,
                EXIT_CLIENT_ERROR,
            )
        overrides = {"url": url} if url is not None else {}
        try:
            settings = ClientSettings(**overrides)
        except (ValueError, TypeError) as exc:
            die("Configuration Error", str(exc), 400, EXIT_CLIENT_ERROR)

        if settings.url is None:
            die(
                "No Server Found",
                "No running zndraw server found. "
                "Start one with `uv run zndraw` or pass `--url`.",
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
