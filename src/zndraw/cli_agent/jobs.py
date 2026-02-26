from __future__ import annotations

from typing import Annotated

import typer
from zndraw_joblib.schemas import PaginatedResponse, TaskResponse

from .connection import get_connection
from .output import json_print

jobs_app = typer.Typer(name="jobs", help="Job operations")


@jobs_app.command("status")
def status(
    ctx: typer.Context,
    task_id: Annotated[str, typer.Argument(help="Task ID")],
) -> None:
    """Get the status of a task by ID."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/joblib/tasks/{task_id}")
    json_print(TaskResponse.model_validate(resp.json()))


@jobs_app.command("list")
def list_jobs(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    status: Annotated[str | None, typer.Option(help="Filter by status")] = None,
) -> None:
    """List tasks for a room."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    params: dict = {}
    if status:
        params["status"] = status
    resp = conn.get(f"/v1/joblib/rooms/{room}/tasks", params=params)
    json_print(PaginatedResponse[TaskResponse].model_validate(resp.json()))
