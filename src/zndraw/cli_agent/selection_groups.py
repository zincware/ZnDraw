from __future__ import annotations

import json
from typing import Annotated

import typer

from zndraw.schemas import (
    SelectionGroupResponse,
    SelectionGroupsListResponse,
    SelectionGroupUpdateRequest,
    StatusResponse,
)

from .connection import cli_error_handler, get_zndraw
from .output import json_print

selection_groups_app = typer.Typer()


@selection_groups_app.command("list")
def list_selection_groups(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List all selection groups."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/selection-groups", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(SelectionGroupsListResponse.model_validate(resp.json()))
        vis.disconnect()


@selection_groups_app.command("get")
def get_selection_group(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    name: Annotated[str, typer.Argument(help="Selection group name")],
) -> None:
    """Get a selection group by name."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/selection-groups/{name}",
            headers=vis.api._headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(SelectionGroupResponse.model_validate(resp.json()))
        vis.disconnect()


@selection_groups_app.command("set")
def set_selection_group(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    name: Annotated[str, typer.Argument(help="Selection group name")],
    data: Annotated[str, typer.Option("--data", help="Selection data as JSON string")],
) -> None:
    """Set a selection group."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        result = vis.api.set_selection_group(name, json.loads(data))
        json_print(StatusResponse.model_validate(result))
        vis.disconnect()


@selection_groups_app.command("delete")
def delete_selection_group(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    name: Annotated[str, typer.Argument(help="Selection group name")],
) -> None:
    """Delete a selection group."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/selection-groups/{name}",
            headers=vis.api._headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
