from __future__ import annotations

import json
from typing import Annotated

import typer

from zndraw.schemas import (
    GeometriesResponse,
    GeometryCreateRequest,
    GeometryResponse,
    StatusResponse,
)

from .connection import get_connection
from .output import json_print

geometries_app = typer.Typer(name="geometries", help="Geometry operations")


@geometries_app.command("list")
def list_geometries(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List geometry keys for a room."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/rooms/{room}/geometries")
    json_print(GeometriesResponse.model_validate(resp.json()))


@geometries_app.command("get")
def get(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Geometry key")],
) -> None:
    """Get a geometry by key."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/rooms/{room}/geometries/{key}")
    json_print(GeometryResponse.model_validate(resp.json()))


@geometries_app.command("set")
def set_geometry(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Geometry key")],
    type_name: Annotated[str, typer.Option("--type", help="Geometry type name")],
    data: Annotated[str, typer.Option(help="Geometry data as JSON string")],
) -> None:
    """Set a geometry by key."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    request = GeometryCreateRequest(type=type_name, data=json.loads(data))
    resp = conn.put(f"/v1/rooms/{room}/geometries/{key}", json=request.model_dump())
    json_print(StatusResponse.model_validate(resp.json()))


@geometries_app.command("delete")
def delete(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Geometry key")],
) -> None:
    """Delete a geometry by key."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.delete(f"/v1/rooms/{room}/geometries/{key}")
    json_print(StatusResponse.model_validate(response.json()))
