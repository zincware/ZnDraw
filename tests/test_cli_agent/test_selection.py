"""Tests for CLI selection and selection-groups commands."""

from __future__ import annotations

from typer.testing import CliRunner

from zndraw.schemas import (
    GeometrySelectionResponse,
    SelectionGroupResponse,
    SelectionGroupsListResponse,
    StatusResponse,
)

from .conftest import invoke_cli


def test_selection_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """selection get should return a GeometrySelectionResponse."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection", "get", test_room],
    )
    resp = GeometrySelectionResponse.model_validate(data)
    assert isinstance(resp.selection, list)


def test_selection_set_and_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """selection set then get should reflect the change."""
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection", "set", test_room, "0", "1", "2"],
    )
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection", "get", test_room],
    )
    resp = GeometrySelectionResponse.model_validate(data)
    assert resp.selection == [0, 1, 2]


def test_selection_clear(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """selection clear should empty the selection."""
    # Set something first
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection", "set", test_room, "0", "1"],
    )
    # Clear it
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection", "clear", test_room],
    )
    StatusResponse.model_validate(data)

    # Verify it's cleared
    get_data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection", "get", test_room],
    )
    resp = GeometrySelectionResponse.model_validate(get_data)
    assert resp.selection == []


def test_selection_with_geometry_option(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """selection get with --geometry should use the specified key."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection", "get", test_room, "--geometry", "particles"],
    )
    resp = GeometrySelectionResponse.model_validate(data)
    assert resp.key == "particles"


def test_selection_groups_list_empty(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """selection-groups list should return empty initially."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection-groups", "list", test_room],
    )
    resp = SelectionGroupsListResponse.model_validate(data)
    assert isinstance(resp.items, dict)


def test_selection_groups_set_and_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """selection-groups set then get should persist the group."""
    # Read selection_groups.py to understand the set command args
    # set command takes room, name, and --selections as JSON string
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        [
            "selection-groups",
            "set",
            test_room,
            "my-group",
            "--data",
            '{"particles": [0, 1, 2]}',
        ],
    )
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection-groups", "get", test_room, "my-group"],
    )
    resp = SelectionGroupResponse.model_validate(data)
    assert resp.group["particles"] == [0, 1, 2]


def test_selection_groups_delete(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """selection-groups delete should remove the group."""
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        [
            "selection-groups",
            "set",
            test_room,
            "to-delete",
            "--data",
            '{"particles": [0]}',
        ],
    )
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["selection-groups", "delete", test_room, "to-delete"],
    )
    StatusResponse.model_validate(data)
