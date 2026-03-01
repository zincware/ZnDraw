from __future__ import annotations

from typing import Annotated

import typer
from zndraw_joblib.schemas import PaginatedResponse, TaskResponse

from .connection import cli_error_handler, get_zndraw
from .output import json_print

jobs_app = typer.Typer(name="jobs", help="Job operations")


@jobs_app.command("status")
def status(
    ctx: typer.Context,
    task_id: Annotated[str, typer.Argument(help="Task ID")],
) -> None:
    """Get the status of a task by ID."""
    with cli_error_handler():
        # Task lookup is room-agnostic but we need a ZnDraw for auth.
        # Use a dummy room — get_task doesn't need room_id.
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room="_")
        data = vis.api.get_task(task_id)
        json_print(TaskResponse.model_validate(data))
        vis.disconnect()


@jobs_app.command("list")
def list_jobs(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    status: Annotated[str | None, typer.Option(help="Filter by status")] = None,
) -> None:
    """List tasks for a room."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.list_tasks(status=status)
        json_print(PaginatedResponse[TaskResponse].model_validate(data))
        vis.disconnect()
