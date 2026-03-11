from __future__ import annotations

import pathlib
from typing import Annotated

import typer

from zndraw.schemas import (
    CollectionResponse,
    FigureCreateResponse,
    FigureData,
    FigureResponse,
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

figures_app = typer.Typer(name="figures", help="Figure operations")


@figures_app.command("list")
def list_figures(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """List figures for a room."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/figures", headers=vis.api.get_headers()
        )
        vis.api.raise_for_status(resp)
        json_print(CollectionResponse[str].model_validate(resp.json()))
        vis.disconnect()


@figures_app.command("get")
def get(
    key: Annotated[str | None, typer.Argument(help="Figure key")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get a figure by key."""
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Figure key is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/figures/{key}", headers=vis.api.get_headers()
        )
        vis.api.raise_for_status(resp)
        json_print(FigureResponse.model_validate(resp.json()))
        vis.disconnect()


@figures_app.command("set")
def set_figure(
    key: Annotated[str | None, typer.Argument(help="Figure key")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    file: Annotated[
        str | None, typer.Option(help="Path to JSON file with plotly figure data")
    ] = None,
    data: Annotated[str | None, typer.Option(help="Inline JSON figure data")] = None,
) -> None:
    """Set a figure by key from a JSON file or inline JSON data."""
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Figure key is required")
        if file is None and data is None:
            raise typer.BadParameter("Either --file or --data is required")
        vis = get_zndraw(url, token, room)
        if data is not None:
            file_contents = data
        else:
            file_contents = pathlib.Path(file).read_text()  # type: ignore[arg-type]
        figure_data = FigureData(data=file_contents)
        result = vis.api.set_figure(key, figure_data.model_dump())
        json_print(FigureCreateResponse.model_validate(result))
        vis.disconnect()


@figures_app.command("delete")
def delete(
    key: Annotated[str | None, typer.Argument(help="Figure key")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Delete a figure by key."""
    with cli_error_handler():
        room = resolve_room(room)
        if key is None:
            raise typer.BadParameter("Figure key is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/figures/{key}", headers=vis.api.get_headers()
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
