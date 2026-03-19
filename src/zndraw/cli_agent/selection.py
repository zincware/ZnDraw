from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import (
    GeometrySelectionResponse,
    StatusResponse,
)

from .connection import (
    PasswordOpt,
    RoomOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

selection_app = typer.Typer()


@selection_app.command("get")
def get_selection(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Get the current selection for a geometry."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        data = vis.api.get_selection(geometry)
        json_print(GeometrySelectionResponse.model_validate(data))
        vis.disconnect()


@selection_app.command("set")
def set_selection(
    indices: Annotated[
        list[int] | None, typer.Argument(help="Indices to select")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Set the selection for a geometry."""
    with cli_error_handler():
        room = resolve_room(room)
        if indices is None:
            raise typer.BadParameter("Indices are required")
        vis = get_zndraw(url, token, room, user, password)
        data = vis.api.update_selection(geometry, list(indices))
        json_print(StatusResponse.model_validate(data))
        vis.disconnect()


@selection_app.command("clear")
def clear_selection(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
    geometry: Annotated[str, typer.Option(help="Geometry key")] = "particles",
) -> None:
    """Clear the selection for a geometry."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        data = vis.api.update_selection(geometry, [])
        json_print(StatusResponse.model_validate(data))
        vis.disconnect()
