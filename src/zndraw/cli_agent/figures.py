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

from .connection import get_connection
from .output import json_print

figures_app = typer.Typer(name="figures", help="Figure operations")


@figures_app.command("list")
def list_figures(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List figures for a room."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/rooms/{room}/figures")
    json_print(CollectionResponse[str].model_validate(resp.json()))


@figures_app.command("get")
def get(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Figure key")],
) -> None:
    """Get a figure by key."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/rooms/{room}/figures/{key}")
    json_print(FigureResponse.model_validate(resp.json()))


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
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    file_contents = pathlib.Path(file).read_text()
    request = FigureCreateRequest(figure=FigureData(data=file_contents))
    resp = conn.post(f"/v1/rooms/{room}/figures/{key}", json=request.model_dump())
    json_print(FigureCreateResponse.model_validate(resp.json()))


@figures_app.command("delete")
def delete(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    key: Annotated[str, typer.Argument(help="Figure key")],
) -> None:
    """Delete a figure by key."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.delete(f"/v1/rooms/{room}/figures/{key}")
    json_print(StatusResponse.model_validate(response.json()))
