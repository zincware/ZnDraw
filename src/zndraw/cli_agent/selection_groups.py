from __future__ import annotations

import json
from typing import Annotated

import typer

from zndraw.schemas import (
    SelectionGroupResponse,
    SelectionGroupsListResponse,
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

selection_groups_app = typer.Typer()


@selection_groups_app.command("list")
def list_selection_groups(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """List all selection groups."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/selection-groups", headers=vis.api._headers()
        )
        vis.api.raise_for_status(resp)
        json_print(SelectionGroupsListResponse.model_validate(resp.json()))
        vis.disconnect()


@selection_groups_app.command("get")
def get_selection_group(
    name: Annotated[str | None, typer.Argument(help="Selection group name")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get a selection group by name."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Selection group name is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/selection-groups/{name}",
            headers=vis.api._headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(SelectionGroupResponse.model_validate(resp.json()))
        vis.disconnect()


@selection_groups_app.command("set")
def set_selection_group(
    name: Annotated[str | None, typer.Argument(help="Selection group name")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    data: Annotated[
        str | None, typer.Option("--data", help="Selection data as JSON string")
    ] = None,
) -> None:
    """Set a selection group.

    ``--data`` accepts a JSON object mapping geometry keys to index lists,
    e.g. ``{"particles": [0,1,2]}``.
    """
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Selection group name is required")
        if data is None:
            raise typer.BadParameter("--data is required")
        vis = get_zndraw(url, token, room)
        parsed = json.loads(data)
        result = vis.api.set_selection_group(name, parsed)
        json_print(StatusResponse.model_validate(result))
        vis.disconnect()


@selection_groups_app.command("delete")
def delete_selection_group(
    name: Annotated[str | None, typer.Argument(help="Selection group name")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Delete a selection group."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Selection group name is required")
        vis = get_zndraw(url, token, room)
        resp = vis.api.http.delete(
            f"/v1/rooms/{vis.room}/selection-groups/{name}",
            headers=vis.api._headers(),
        )
        vis.api.raise_for_status(resp)
        json_print(StatusResponse.model_validate(resp.json()))
        vis.disconnect()
