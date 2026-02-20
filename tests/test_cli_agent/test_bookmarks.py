"""Tests for CLI bookmarks commands."""

from __future__ import annotations

from typer.testing import CliRunner

from zndraw.schemas import BookmarksResponse, StatusResponse

from .conftest import invoke_cli


def test_bookmarks_list_empty(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """bookmarks list should return empty initially."""
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["bookmarks", "list", test_room]
    )
    resp = BookmarksResponse.model_validate(data)
    assert isinstance(resp.items, dict)


def test_bookmarks_set_with_index(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """bookmarks set with explicit index should create a bookmark."""
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["bookmarks", "set", test_room, "0", "--label", "Start"],
    )
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["bookmarks", "list", test_room]
    )
    resp = BookmarksResponse.model_validate(data)
    assert "0" in resp.items
    assert resp.items["0"] == "Start"


def test_bookmarks_set_default_index(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """bookmarks set without index should bookmark current step."""
    # Ensure step is 0
    invoke_cli(cli_runner, server_url, auth_token, ["step", "set", test_room, "0"])
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["bookmarks", "set", test_room, "--label", "Current"],
    )
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["bookmarks", "list", test_room]
    )
    resp = BookmarksResponse.model_validate(data)
    assert "0" in resp.items
    assert resp.items["0"] == "Current"


def test_bookmarks_set_default_label(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """bookmarks set without --label should use auto-generated label."""
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["bookmarks", "set", test_room, "0"],
    )
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["bookmarks", "list", test_room]
    )
    resp = BookmarksResponse.model_validate(data)
    assert "0" in resp.items
    assert "Frame 0" == resp.items["0"]


def test_bookmarks_delete_with_index(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """bookmarks delete should remove a bookmark."""
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["bookmarks", "set", test_room, "0", "--label", "ToDelete"],
    )
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["bookmarks", "delete", test_room, "0"],
    )
    StatusResponse.model_validate(data)

    after = invoke_cli(
        cli_runner, server_url, auth_token, ["bookmarks", "list", test_room]
    )
    resp = BookmarksResponse.model_validate(after)
    assert "0" not in resp.items
