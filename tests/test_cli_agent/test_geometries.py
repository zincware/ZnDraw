"""Tests for CLI geometries commands."""

from __future__ import annotations

import json

from typer.testing import CliRunner

from zndraw.cli_agent import app
from zndraw.schemas import GeometryResponse, StatusResponse

from .conftest import invoke_cli


def test_geometries_list(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """geometries list should return a compact summary."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["geometries", "list", "--room", test_room],
    )
    assert isinstance(data, list)
    assert len(data) > 0
    keys = {entry["key"] for entry in data}
    assert "particles" in keys
    for entry in data:
        assert "key" in entry
        assert "type" in entry
        assert "active" in entry
        assert "owner" in entry


def test_geometries_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """geometries get should return a GeometryResponse."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["geometries", "get", "--room", test_room, "particles"],
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
            "--room",
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
        ["geometries", "get", "--room", test_room, "test-geo"],
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
            "--room",
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
        ["geometries", "delete", "--room", test_room, "to-delete"],
    )
    StatusResponse.model_validate(data)


def test_geometries_types(cli_runner: CliRunner) -> None:
    """geometries types should list available type names."""
    result = cli_runner.invoke(app, ["geometries", "types"])
    assert result.exit_code == 0, result.stdout
    data = json.loads(result.stdout)
    assert isinstance(data, list)
    assert "Sphere" in data
    assert "Bond" in data


def test_geometries_describe(cli_runner: CliRunner) -> None:
    """geometries describe should return schema and defaults for a type."""
    result = cli_runner.invoke(app, ["geometries", "describe", "Sphere"])
    assert result.exit_code == 0, result.stdout
    data = json.loads(result.stdout)
    assert data["name"] == "Sphere"
    assert "schema" in data
    assert "defaults" in data


def test_geometries_describe_nonexistent(cli_runner: CliRunner) -> None:
    """geometries describe should fail for unknown type."""
    result = cli_runner.invoke(app, ["geometries", "describe", "NoSuchType"])
    assert result.exit_code != 0
