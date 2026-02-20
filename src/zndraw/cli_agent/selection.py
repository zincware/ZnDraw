from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import (
    GeometrySelectionResponse,
    SelectionUpdateRequest,
    StatusResponse,
)

from .connection import get_connection
from .output import json_print

selection_app = typer.Typer()


@selection_app.command("get")
def get_selection(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Get the current selection for a geometry."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}/geometries/{geometry}/selection")
    json_print(GeometrySelectionResponse.model_validate(response.json()))


@selection_app.command("set")
def set_selection(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    indices: Annotated[list[int], typer.Argument(help="Indices to select")],
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Set the selection for a geometry."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    request = SelectionUpdateRequest(indices=indices)
    response = conn.put(
        f"/v1/rooms/{room}/geometries/{geometry}/selection",
        json=request.model_dump(),
    )
    json_print(StatusResponse.model_validate(response.json()))


@selection_app.command("clear")
def clear_selection(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Clear the selection for a geometry."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    request = SelectionUpdateRequest(indices=[])
    response = conn.put(
        f"/v1/rooms/{room}/geometries/{geometry}/selection",
        json=request.model_dump(),
    )
    json_print(StatusResponse.model_validate(response.json()))
