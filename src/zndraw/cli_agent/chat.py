from __future__ import annotations

from typing import Annotated

import typer

from zndraw.schemas import MessageResponse, MessagesResponse

from .connection import (
    PasswordOpt,
    RoomOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

chat_app = typer.Typer()


@chat_app.command("list")
def list_messages(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    limit: Annotated[
        int | None, typer.Option(help="Maximum number of messages")
    ] = None,
    before: Annotated[
        str | None, typer.Option(help="Return messages before this timestamp")
    ] = None,
) -> None:
    """List chat messages."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room, user, password)
        data = vis.api.list_chat_messages(limit=limit, before=before)
        json_print(MessagesResponse.model_validate(data))
        vis.disconnect()


@chat_app.command("send")
def send_message(
    message: Annotated[str | None, typer.Argument(help="Message content")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
) -> None:
    """Send a chat message."""
    with cli_error_handler():
        room = resolve_room(room)
        if message is None:
            raise typer.BadParameter("Message content is required")
        # Reverse zsh shell escaping of ! → \!
        message = message.replace("\\!", "!")
        # Standard escape processing
        message = message.replace("\\n", "\n").replace("\\t", "\t")
        vis = get_zndraw(url, token, room, user, password)
        data = vis.api.create_chat_message(message)
        json_print(MessageResponse.model_validate(data))
        vis.disconnect()
