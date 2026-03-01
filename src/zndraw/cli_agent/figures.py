from __future__ import annotations

import pathlib
from typing import Annotated

import typer

from zndraw.schemas import (
    CollectionResponse,
    FigureCreateRequest,
    FigureCreateResponse,
    FigureData,
    FigureResponse,
    StatusResponse,
)

from .connection import cli_error_handler, get_zndraw
from .output import json_print

figures_app = typer.Typer(name="figures", help="Figure operations")


@figures_app.command("list")
def list_figures(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List figures for a room."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/figures", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(CollectionResponse[str].model_validate(resp.json()))
        vis.disconnect()


@figures_app.command("get")
def get(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Figure key")],
) -> None:
    """Get a figure by key."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/figures/{key}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(FigureResponse.model_validate(resp.json()))
        vis.disconnect()


@figures_app.command("set")
def set_figure(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Figure key")],
    file: Annotated[
        str, typer.Option(help="Path to JSON file with plotly figure data")
    ],
) -> None:
    """Set a figure by key from a JSON file."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        file_contents = pathlib.Path(file).read_text()
        figure_data = FigureData(data=file_contents)
        vis.api.set_figure(key, figure_data.model_dump())
        json_print(FigureCreateResponse(key=key, status="ok"))
        vis.disconnect()


@figures_app.command("delete")
def delete(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Figure key")],
) -> None:
    """Delete a figure by key."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/figures/{key}", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
