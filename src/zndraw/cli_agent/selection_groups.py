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

from .connection import get_connection
from .output import json_print

selection_groups_app = typer.Typer()


@selection_groups_app.command("list")
def list_selection_groups(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """List all selection groups."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}/selection-groups")
    json_print(SelectionGroupsListResponse.model_validate(response.json()))


@selection_groups_app.command("get")
def get_selection_group(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    name: Annotated[str, typer.Argument(help="Selection group name")],
) -> None:
    """Get a selection group by name."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.get(f"/v1/rooms/{room}/selection-groups/{name}")
    json_print(SelectionGroupResponse.model_validate(response.json()))


@selection_groups_app.command("set")
def set_selection_group(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    name: Annotated[str, typer.Argument(help="Selection group name")],
    data: Annotated[str, typer.Option("--data", help="Selection data as JSON string")],
) -> None:
    """Set a selection group."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    request = SelectionGroupUpdateRequest(selections=json.loads(data))
    response = conn.put(
        f"/v1/rooms/{room}/selection-groups/{name}", json=request.model_dump()
    )
    json_print(StatusResponse.model_validate(response.json()))


@selection_groups_app.command("delete")
def delete_selection_group(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    name: Annotated[str, typer.Argument(help="Selection group name")],
) -> None:
    """Delete a selection group."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    response = conn.delete(f"/v1/rooms/{room}/selection-groups/{name}")
    json_print(StatusResponse.model_validate(response.json()))
