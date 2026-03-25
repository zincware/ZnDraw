"""Tests for bug fixes from the 2026-03-03 SKILL.md design."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from zndraw.cli_agent import app
from zndraw_joblib.schemas import TaskResponse

from .conftest import invoke_cli

if TYPE_CHECKING:
    from typer.testing import CliRunner


def test_jobs_status_without_room(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """jobs status should work without a room (uses get_connection)."""
    # First, submit a job to get a task_id
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["extensions", "run", "--room", test_room, "@internal:selections:All"],
    )
    task_id = data["id"]

    # Now query jobs status — this used to fail with 404 because room="_"
    status_data = invoke_cli(
        cli_runner, server_url, auth_token, ["jobs", "status", task_id]
    )
    resp = TaskResponse.model_validate(status_data)
    assert str(resp.id) == task_id


def test_chat_send_exclamation(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    r"""chat send should strip zsh-style \! escaping."""
    # Simulate what zsh does: the shell turns ! into \!
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["chat", "send", "--room", test_room, "Hello\\!"],
    )
    assert data["content"] == "Hello!"


def test_chat_send_no_false_positive(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """chat send should not mangle messages without backslash-bang."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["chat", "send", "--room", test_room, "Normal message"],
    )
    assert data["content"] == "Normal message"


def test_frames_extend_positional_file(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames extend should accept FILE as a positional argument."""
    with tempfile.NamedTemporaryFile(suffix=".xyz", mode="w", delete=False) as f:
        f.write("1\nproperties=species:S:1:pos:R:3\nH 0.0 0.0 0.0\n")
        xyz_path = f.name

    try:
        # Positional (no --file flag)
        from zndraw.schemas import FrameBulkResponse

        data = invoke_cli(
            cli_runner,
            server_url,
            auth_token,
            ["frames", "extend", "--room", test_room, xyz_path],
        )
        resp = FrameBulkResponse.model_validate(data)
        assert resp.total >= 1
    finally:
        Path(xyz_path).unlink(missing_ok=True)


def test_figures_set_with_data_option(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """figures set --data should accept inline JSON."""
    figure_json = json.dumps({"data": [{"x": [1, 2], "y": [3, 4], "type": "scatter"}]})
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["figures", "set", "--room", test_room, "test-fig", "--data", figure_json],
    )
    assert data["key"] == "test-fig"


def test_figures_set_requires_file_or_data(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """figures set without --file or --data should fail."""
    result = cli_runner.invoke(
        app,
        [
            "figures",
            "set",
            "--url",
            server_url,
            "--token",
            auth_token,
            "--room",
            test_room,
            "test-fig",
        ],
    )
    assert result.exit_code != 0


def test_preset_list_json_output(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """preset list should output valid JSON (not Python repr)."""
    result = cli_runner.invoke(
        app,
        [
            "presets",
            "list",
            "--url",
            server_url,
            "--token",
            auth_token,
            "--room",
            test_room,
        ],
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    data = json.loads(result.stdout)
    assert isinstance(data, list)
    # Each entry should be a compact summary with name + description
    for entry in data:
        assert isinstance(entry, dict)
        assert "name" in entry
        assert "description" in entry
        assert "rules" not in entry  # compact, not full dump
