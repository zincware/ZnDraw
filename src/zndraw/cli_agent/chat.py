from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import MessageCreate, MessageResponse, MessagesResponse

from .connection import get_connection
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
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    params: dict = {}
    if limit is not None:
        params["limit"] = limit
    if before is not None:
        params["before"] = before
    response = conn.get(f"/v1/rooms/{room}/chat/messages", params=params)
    json_print(MessagesResponse.model_validate(response.json()))


@chat_app.command("send")
def send_message(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    message: Annotated[str, typer.Argument(help="Message content")],
) -> None:
    """Send a chat message."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    request = MessageCreate(content=message)
    response = conn.post(f"/v1/rooms/{room}/chat/messages", json=request.model_dump())
    json_print(MessageResponse.model_validate(response.json()))
