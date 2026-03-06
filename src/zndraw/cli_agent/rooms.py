from __future__ import annotations

import uuid
from typing import Annotated

import typer

from zndraw.schemas import (
    CollectionResponse,
    RoomCreate,
    RoomCreateResponse,
    RoomResponse,
)

from .connection import get_connection
from .output import json_print

rooms_app = typer.Typer()


@rooms_app.command("list")
def list_rooms(
    ctx: typer.Context,
    search: Annotated[str | None, typer.Option(help="Search query")] = None,
) -> None:
    """List all rooms."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    params = {"search": search} if search is not None else {}
    response = conn.get("/v1/rooms", params=params)
    json_print(CollectionResponse[RoomResponse].model_validate(response.json()))


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
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}")
    json_print(RoomResponse.model_validate(response.json()))
