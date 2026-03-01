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
    """List geometries for a room (compact summary)."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/geometries", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        full = GeometriesResponse.model_validate(resp.json())
        summary = [
            {
                "key": key,
                "type": entry.type,
                "active": entry.data.get("active", True),
                "owner": entry.data.get("owner"),
            }
            for key, entry in full.items.items()
        ]
        json_print(summary)
        vis.disconnect()


@geometries_app.command("types")
def types() -> None:
    """List available geometry type names."""
    from zndraw.geometries import geometries

    json_print(list(geometries.keys()))


@geometries_app.command("describe")
def describe(
    type_name: Annotated[str, typer.Argument(help="Geometry type name")],
) -> None:
    """Show schema and defaults for a geometry type."""
    from zndraw.geometries import geometries

    model = geometries.get(type_name)
    if model is None:
        raise typer.BadParameter(
            f"Unknown geometry type {type_name!r}. "
            f"Available: {', '.join(geometries.keys())}"
        )
    json_print(
        {
            "name": type_name,
            "schema": model.model_json_schema(),
            "defaults": model().model_dump(),
        }
    )


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
