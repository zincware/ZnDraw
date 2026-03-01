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

from .connection import cli_error_handler, get_zndraw
from .output import json_print

geometries_app = typer.Typer(name="geometries", help="Geometry operations")


@geometries_app.command("list")
def list_geometries(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List geometry keys for a room."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/geometries", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(GeometriesResponse.model_validate(resp.json()))
        vis.disconnect()


@geometries_app.command("get")
def get(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Geometry key")],
) -> None:
    """Get a geometry by key."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/geometries/{key}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(GeometryResponse.model_validate(resp.json()))
        vis.disconnect()


@geometries_app.command("set")
def set_geometry(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Geometry key")],
    type_name: Annotated[str, typer.Option("--type", help="Geometry type name")],
    data: Annotated[str, typer.Option(help="Geometry data as JSON string")],
) -> None:
    """Set a geometry by key."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        request = GeometryCreateRequest(type=type_name, data=json.loads(data))
        resp = vis.api.http.put(
            f"/v1/rooms/{vis.room}/geometries/{key}",
            json=request.model_dump(),
            headers=vis.api._headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()


@geometries_app.command("delete")
def delete(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Geometry key")],
) -> None:
    """Delete a geometry by key."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/geometries/{key}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
