from __future__ import annotations

from typing import Annotated

import typer
from zndraw_joblib.schemas import PaginatedResponse, TaskResponse

from .connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    cli_error_handler,
    get_connection,
    get_zndraw,
    resolve_room,
)
from .output import json_print

jobs_app = typer.Typer(name="jobs", help="Job operations")


@jobs_app.command("status")
def status(
    task_id: Annotated[str, typer.Argument(help="Task ID")],
    url: UrlOpt = None,
    token: TokenOpt = None,
) -> None:
    """Get the status of a task by ID."""
    with cli_error_handler():
        conn = get_connection(url, token)
        resp = conn.get(f"/v1/joblib/tasks/{task_id}")
        json_print(TaskResponse.model_validate(resp.json()))


@jobs_app.command("list")
def list_jobs(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    status: Annotated[str | None, typer.Option(help="Filter by status")] = None,
) -> None:
    """List tasks for a room."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        data = vis.api.list_tasks(status=status)
        json_print(PaginatedResponse[TaskResponse].model_validate(data))
        vis.disconnect()
