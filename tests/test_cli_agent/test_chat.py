"""Tests for CLI chat commands."""

from __future__ import annotations

from typer.testing import CliRunner

from zndraw.schemas import MessageResponse, MessagesResponse

from .conftest import invoke_cli


def test_chat_list_empty(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """chat list on a fresh room should return empty items."""
    data = invoke_cli(cli_runner, server_url, auth_token, ["chat", "list", test_room])
    resp = MessagesResponse.model_validate(data)
    assert isinstance(resp.items, list)


def test_chat_send(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """chat send should create a message and return it."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["chat", "send", test_room, "Hello from CLI"],
    )
    resp = MessageResponse.model_validate(data)
    assert resp.content == "Hello from CLI"


def test_chat_send_then_list(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """Sent message should appear in chat list."""
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["chat", "send", test_room, "Test message"],
    )
    data = invoke_cli(cli_runner, server_url, auth_token, ["chat", "list", test_room])
    resp = MessagesResponse.model_validate(data)
    contents = [m.content for m in resp.items]
    assert "Test message" in contents
