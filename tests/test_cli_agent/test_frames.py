"""Tests for CLI frames commands."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from zndraw.cli_agent import app
from zndraw.schemas import FrameBulkResponse, StatusResponse, StepResponse

from .conftest import invoke_cli

if TYPE_CHECKING:
    from typer.testing import CliRunner


def test_frames_count(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames count should report at least 1 frame."""
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["frames", "count", "--room", test_room]
    )
    resp = StepResponse.model_validate(data)
    assert resp.total_frames >= 1


def test_frames_get_with_index(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames get --room <room> 0 should return frame data."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["frames", "get", "--room", test_room, "0"],
    )
    assert "arrays.numbers" in data
    assert "arrays.positions" in data
    assert "cell" in data
    assert "pbc" in data


def test_frames_get_default_index(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames get --room <room> (no index) should use current step."""
    invoke_cli(
        cli_runner, server_url, auth_token, ["step", "set", "--room", test_room, "0"]
    )
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["frames", "get", "--room", test_room]
    )
    assert "arrays.numbers" in data
    assert "arrays.positions" in data


def test_frames_get_xyz_format(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames get with --format xyz should output XYZ text."""
    result = cli_runner.invoke(
        app,
        [
            "frames",
            "get",
            "--url",
            server_url,
            "--token",
            auth_token,
            "--room",
            test_room,
            "0",
            "--format",
            "xyz",
        ],
    )
    assert result.exit_code == 0, result.stderr
    lines = result.stdout.strip().split("\n")
    assert len(lines) >= 2


def test_frames_list(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames list should return a list of frames."""
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["frames", "list", "--room", test_room]
    )
    assert isinstance(data, list)
    assert len(data) >= 1
    assert "arrays.numbers" in data[0]


def test_frames_extend_and_count(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """Extending with a file should increase frame count."""
    before = invoke_cli(
        cli_runner, server_url, auth_token, ["frames", "count", "--room", test_room]
    )
    before_resp = StepResponse.model_validate(before)

    with tempfile.NamedTemporaryFile(suffix=".xyz", mode="w", delete=False) as f:
        f.write("2\nproperties=species:S:1:pos:R:3\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
        xyz_path = f.name

    try:
        data = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "extend", "--room", test_room, xyz_path],
        )
        FrameBulkResponse.model_validate(data)

        after = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "count", "--room", test_room],
        )
        after_resp = StepResponse.model_validate(after)
        assert after_resp.total_frames == before_resp.total_frames + 1
    finally:
        Path(xyz_path).unlink(missing_ok=True)


def test_frames_delete_with_index(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames delete should remove a frame."""
    with tempfile.NamedTemporaryFile(suffix=".xyz", mode="w", delete=False) as f:
        f.write("1\nproperties=species:S:1:pos:R:3\nH 0.0 0.0 0.0\n")
        xyz_path = f.name

    try:
        invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "extend", "--room", test_room, xyz_path],
        )
        before = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "count", "--room", test_room],
        )
        before_resp = StepResponse.model_validate(before)
        last_idx = before_resp.total_frames - 1
        data = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "delete", "--room", test_room, str(last_idx)],
        )
        StatusResponse.model_validate(data)

        after = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "count", "--room", test_room],
        )
        after_resp = StepResponse.model_validate(after)
        assert after_resp.total_frames == before_resp.total_frames - 1
    finally:
        Path(xyz_path).unlink(missing_ok=True)


def test_frames_delete_default_index(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames delete with no index should delete the frame at current step."""
    with tempfile.NamedTemporaryFile(suffix=".xyz", mode="w", delete=False) as f:
        f.write("1\nproperties=species:S:1:pos:R:3\nH 0.0 0.0 0.0\n")
        xyz_path = f.name

    try:
        invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "extend", "--room", test_room, xyz_path],
        )
        before = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "count", "--room", test_room],
        )
        before_resp = StepResponse.model_validate(before)
        invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            [
                "step",
                "set",
                "--room",
                test_room,
                str(before_resp.total_frames - 1),
            ],
        )
        data = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "delete", "--room", test_room],
        )
        StatusResponse.model_validate(data)

        after = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "count", "--room", test_room],
        )
        after_resp = StepResponse.model_validate(after)
        assert after_resp.total_frames == before_resp.total_frames - 1
    finally:
        Path(xyz_path).unlink(missing_ok=True)
