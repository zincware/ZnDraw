"""Tests for CLI geometries commands."""

from __future__ import annotations

from typer.testing import CliRunner

from zndraw.schemas import GeometriesResponse, GeometryResponse, StatusResponse

from .conftest import invoke_cli


def test_geometries_list(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """geometries list should return a list of geometry key strings."""
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["geometries", "list", test_room]
    )
    resp = GeometriesResponse.model_validate(data)
    assert "particles" in resp.items


def test_geometries_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """geometries get should return a GeometryResponse."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["geometries", "get", test_room, "particles"],
    )
    resp = GeometryResponse.model_validate(data)
    assert resp.key == "particles"


def test_geometries_set_and_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """geometries set should create a new geometry."""
    set_data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        [
            "geometries",
            "set",
            test_room,
            "test-geo",
            "--type",
            "particles",
            "--data",
            '{"positions": [[0, 0, 0]]}',
        ],
    )
    StatusResponse.model_validate(set_data)
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["geometries", "get", test_room, "test-geo"],
    )
    resp = GeometryResponse.model_validate(data)
    assert resp.key == "test-geo"


def test_geometries_delete(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """geometries delete should remove a geometry."""
    invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        [
            "geometries",
            "set",
            test_room,
            "to-delete",
            "--type",
            "particles",
            "--data",
            '{"positions": [[0, 0, 0]]}',
        ],
    )
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["geometries", "delete", test_room, "to-delete"],
    )
    StatusResponse.model_validate(data)
