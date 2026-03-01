from __future__ import annotations

import uuid
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
    cli_error_handler,
    get_connection,
    get_zndraw,
    resolve_token,
    resolve_url,
)
from .output import json_print

rooms_app = typer.Typer()


@rooms_app.command("list")
def list_rooms(
    ctx: typer.Context,
    search: Annotated[str | None, typer.Option(help="Search query")] = None,
) -> None:
    """List all rooms."""
    with cli_error_handler():
        from zndraw.client import APIManager

        url = resolve_url(ctx.obj["url"])
        token = resolve_token(url, ctx.obj["token"])
        api = APIManager(url=url, room_id="", token=token)
        try:
            params: dict[str, str] = {}
            if search is not None:
                params["search"] = search
            resp = api.http.get("/v1/rooms", params=params, headers=api._headers())
            api.raise_for_status(resp)
            json_print(CollectionResponse[RoomResponse].model_validate(resp.json()))
        finally:
            api.close()


@rooms_app.command("create")
def create_room(
    ctx: typer.Context,
    room_id: Annotated[
        str | None, typer.Option(help="Room ID (generated if not given)")
    ] = None,
    copy_from: Annotated[
        str | None, typer.Option("--copy-from", help="Copy from existing room ID")
    ] = None,
) -> None:
    """Create a new room."""
    with cli_error_handler():
        conn = get_connection(ctx.obj["url"], ctx.obj["token"])
        request = RoomCreate(room_id=room_id or str(uuid.uuid4()), copy_from=copy_from)
        response = conn.post("/v1/rooms", json=request.model_dump())
        json_print(RoomCreateResponse.model_validate(response.json()))


@rooms_app.command("info")
def room_info(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Get room info."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        json_print(RoomResponse.model_validate(vis.api.get_room_info()))
        vis.disconnect()


@rooms_app.command("lock")
def lock_room(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Lock a room (prevent edits by non-admin users)."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        vis.locked = True
        json_print(RoomPatchResponse.model_validate(vis.api.get_room_info()))
        vis.disconnect()


@rooms_app.command("unlock")
def unlock_room(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Unlock a room (allow edits again)."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        vis.locked = False
        json_print(RoomPatchResponse.model_validate(vis.api.get_room_info()))
        vis.disconnect()


@rooms_app.command("set-default")
def set_default_room(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID to set as default template")],
) -> None:
    """Set a room as the default template for new rooms."""
    with cli_error_handler():
        conn = get_connection(ctx.obj["url"], ctx.obj["token"])
        response = conn.put("/v1/server-settings/default-room", json={"room_id": room})
        json_print(response.json())
