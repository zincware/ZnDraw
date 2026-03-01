from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import (
    GeometrySelectionResponse,
    SelectionUpdateRequest,
    StatusResponse,
)

from .connection import cli_error_handler, get_zndraw
from .output import json_print

selection_app = typer.Typer()


@selection_app.command("get")
def get_selection(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Get the current selection for a geometry."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.get_selection(geometry)
        json_print(GeometrySelectionResponse.model_validate(data))
        vis.disconnect()


@selection_app.command("set")
def set_selection(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    indices: Annotated[list[int], typer.Argument(help="Indices to select")],
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Set the selection for a geometry."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.update_selection(geometry, list(indices))
        json_print(StatusResponse.model_validate(data))
        vis.disconnect()


@selection_app.command("clear")
def clear_selection(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Clear the selection for a geometry."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.update_selection(geometry, [])
        json_print(StatusResponse.model_validate(data))
        vis.disconnect()
