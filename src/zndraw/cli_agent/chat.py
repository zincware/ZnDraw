from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import MessageResponse, MessagesResponse

from .connection import cli_error_handler, get_zndraw
from .output import json_print

chat_app = typer.Typer()


@chat_app.command("list")
def list_messages(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    limit: Annotated[
        int | None, typer.Option(help="Maximum number of messages")
    ] = None,
    before: Annotated[
        str | None, typer.Option(help="Return messages before this timestamp")
    ] = None,
) -> None:
    """List chat messages."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.list_chat_messages(limit=limit, before=before)
        json_print(MessagesResponse.model_validate(data))
        vis.disconnect()


@chat_app.command("send")
def send_message(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    message: Annotated[str, typer.Argument(help="Message content")],
) -> None:
    """Send a chat message."""
    with cli_error_handler():
        vis = get_zndraw(ctx.obj["url"], ctx.obj["token"], room)
        data = vis.api.create_chat_message(message)
        json_print(MessageResponse.model_validate(data))
        vis.disconnect()
