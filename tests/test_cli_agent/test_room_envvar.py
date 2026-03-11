"""Tests that ZNDRAW_ROOM env var works as a substitute for --room."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest
from zndraw_joblib.schemas import JobSummary, PaginatedResponse, TaskResponse

from zndraw.cli_agent import app
from zndraw.schemas import StepResponse

if TYPE_CHECKING:
    from typer.testing import CliRunner


def _invoke_with_env_room(
    cli_runner: CliRunner,
    server_url: str,
    auth_token: str,
    room: str,
    args: list[str],
):
    """Invoke a CLI command with ZNDRAW_ROOM set instead of --room."""
    cmd_path = args[:2]
    cmd_args = args[2:]
    return cli_runner.invoke(
        app,
        [*cmd_path, "--url", server_url, "--token", auth_token, *cmd_args],
        env={"ZNDRAW_ROOM": room},
    )


def test_extensions_list_via_envvar(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions list should work with ZNDRAW_ROOM instead of --room."""
    result = _invoke_with_env_room(
        cli_runner, server_url, auth_token, test_room, ["extensions", "list"]
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    data = json.loads(result.stdout)
    resp = PaginatedResponse[JobSummary].model_validate(data)
    assert resp.total > 0


def test_extensions_run_via_envvar(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions run should work with ZNDRAW_ROOM instead of --room."""
    result = _invoke_with_env_room(
        cli_runner,
        server_url,
        auth_token,
        test_room,
        ["extensions", "run", "@internal:selections:All"],
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    data = json.loads(result.stdout)
    resp = TaskResponse.model_validate(data)
    assert resp.job_name == "@internal:selections:All"
    assert resp.room_id == test_room


def test_frames_count_via_envvar(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """frames count should work with ZNDRAW_ROOM instead of --room."""
    result = _invoke_with_env_room(
        cli_runner, server_url, auth_token, test_room, ["frames", "count"]
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    data = json.loads(result.stdout)
    resp = StepResponse.model_validate(data)
    assert resp.total_frames >= 1


def test_missing_room_errors(
    cli_runner: CliRunner, server_url: str, auth_token: str
) -> None:
    """Command should fail with helpful error when neither --room nor ZNDRAW_ROOM is set."""
    result = cli_runner.invoke(
        app,
        ["extensions", "list", "--url", server_url, "--token", auth_token],
    )
    assert result.exit_code != 0


@pytest.mark.parametrize(
    "command",
    [
        ["extensions", "list"],
        ["frames", "count"],
        ["selection", "get"],
        ["bookmarks", "list"],
    ],
)
def test_envvar_works_across_commands(
    cli_runner: CliRunner,
    server_url: str,
    auth_token: str,
    test_room: str,
    command: list[str],
) -> None:
    """ZNDRAW_ROOM env var should work across multiple command groups."""
    result = _invoke_with_env_room(
        cli_runner, server_url, auth_token, test_room, command
    )
    assert result.exit_code == 0, (
        f"cmd={command} exit={result.exit_code}\n"
        f"stderr={result.stderr}\nstdout={result.stdout}"
    )
