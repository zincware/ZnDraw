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

from .connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

geometries_app = typer.Typer(name="geometries", help="Geometry operations")


@geometries_app.command("list")
def list_geometries(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """List geometries for a room (compact summary)."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
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
    key: Annotated[str | None, typer.Argument(help="Geometry key")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get a geometry by key."""
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Geometry key is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/geometries/{key}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(GeometryResponse.model_validate(resp.json()))
        vis.disconnect()


@geometries_app.command("set")
def set_geometry(
    key: Annotated[str | None, typer.Argument(help="Geometry key")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    data: Annotated[
        str | None, typer.Option(help="Geometry data as JSON string")
    ] = None,
    type_name: Annotated[
        str | None,
        typer.Option("--type", help="Geometry type (omit for partial update)"),
    ] = None,
) -> None:
    """Set a geometry by key.

    With ``--type``: full replace (PUT).
    Without ``--type``: partial update / deep merge (PATCH).
    """
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Geometry key is required")
        if data is None:
            raise typer.BadParameter("--data is required")
        vis = get_zndraw(url, token, room)
        parsed = json.loads(data)
        if type_name is not None:
            request = GeometryCreateRequest(type=type_name, data=parsed)
            resp = vis.api.http.put(
                f"/v1/rooms/{vis.room}/geometries/{key}",
                json=request.model_dump(),
                headers=vis.api._headers(),
            )
        else:
            resp = vis.api.http.patch(
                f"/v1/rooms/{vis.room}/geometries/{key}",
                json={"data": parsed},
                headers=vis.api._headers(),
            )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()


@geometries_app.command("toggle")
def toggle_geometry(
    key: Annotated[str | None, typer.Argument(help="Geometry key")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Toggle a geometry's active state."""
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Geometry key is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/geometries/{key}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        geom = resp.json()["geometry"]
        active = not geom["data"].get("active", True)
        patch_resp = vis.api.http.patch(
            f"/v1/rooms/{vis.room}/geometries/{key}",
            json={"data": {"active": active}},
            headers=vis.api._headers(),
        )
        vis.api.raise_for_status(patch_resp)
        json_print({"key": key, "active": active})
        vis.disconnect()


@geometries_app.command("set-prop")
def set_prop(
    key: Annotated[str | None, typer.Argument(help="Geometry key")] = None,
    prop: Annotated[
        str | None, typer.Argument(help="Property path (dot notation)")
    ] = None,
    value: Annotated[str | None, typer.Argument(help="Value (JSON-parsed)")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Set a single property on a geometry using dot notation.

    Example: ``geometries set-prop particles material.opacity 0.5``
    """
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Geometry key is required")
        if prop is None:
            raise typer.BadParameter("Property path is required")
        if value is None:
            raise typer.BadParameter("Value is required")
        try:
            parsed_value = json.loads(value)
        except json.JSONDecodeError:
            parsed_value = value

        # Build nested dict from dot path
        parts = prop.split(".")
        data: dict = {}
        current = data
        for part in parts[:-1]:
            current[part] = {}
            current = current[part]
        current[parts[-1]] = parsed_value

        vis = get_zndraw(url, token, room)
        resp = vis.api.http.patch(
            f"/v1/rooms/{vis.room}/geometries/{key}",
            json={"data": data},
            headers=vis.api._headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()


@geometries_app.command("delete")
def delete(
    key: Annotated[str | None, typer.Argument(help="Geometry key")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Delete a geometry by key."""
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Geometry key is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/geometries/{key}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
