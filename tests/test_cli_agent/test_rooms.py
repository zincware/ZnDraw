"""Tests for CLI room commands."""

from __future__ import annotations

from typer.testing import CliRunner

from zndraw.schemas import CollectionResponse, RoomCreateResponse, RoomResponse

from .conftest import invoke_cli


def test_rooms_create(cli_runner: CliRunner, server_url: str, auth_token: str) -> None:
    """rooms create should return a valid RoomCreateResponse."""
    data = invoke_cli(cli_runner, server_url, auth_token, ["rooms", "create"])
    resp = RoomCreateResponse.model_validate(data)
    assert resp.status == "ok"
    assert resp.created is True
    assert len(resp.room_id) > 0


def test_rooms_list(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """rooms list should include the test room."""
    data = invoke_cli(cli_runner, server_url, auth_token, ["rooms", "list"])
    resp = CollectionResponse[RoomResponse].model_validate(data)
    room_ids = [r.id for r in resp.items]
    assert test_room in room_ids


def test_rooms_info(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """rooms info should return valid RoomResponse for the test room."""
    data = invoke_cli(cli_runner, server_url, auth_token, ["rooms", "info", test_room])
    resp = RoomResponse.model_validate(data)
    assert resp.id == test_room
    assert resp.frame_count >= 0


def test_rooms_create_returns_different_ids(
    cli_runner: CliRunner, server_url: str, auth_token: str
) -> None:
    """Each rooms create should produce a unique room ID."""
    data1 = invoke_cli(cli_runner, server_url, auth_token, ["rooms", "create"])
    data2 = invoke_cli(cli_runner, server_url, auth_token, ["rooms", "create"])
    resp1 = RoomCreateResponse.model_validate(data1)
    resp2 = RoomCreateResponse.model_validate(data2)
    assert resp1.room_id != resp2.room_id
