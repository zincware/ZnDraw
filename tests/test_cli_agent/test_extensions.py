"""Tests for CLI extension commands."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

from zndraw_joblib.schemas import (
    JobResponse,
    JobSummary,
    PaginatedResponse,
    TaskResponse,
)

from zndraw.cli_agent import app

from .conftest import invoke_cli

if TYPE_CHECKING:
    from typer.testing import CliRunner


def test_extensions_list(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions list should return a PaginatedResponse of JobSummary."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["extensions", "list", "--room", test_room],
    )
    resp = PaginatedResponse[JobSummary].model_validate(data)
    assert resp.total > 0
    assert len(resp.items) > 0


def test_extensions_list_contains_builtin(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions list should contain built-in extensions like Delete and All."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        ["extensions", "list", "--room", test_room],
    )
    resp = PaginatedResponse[JobSummary].model_validate(data)
    names = {item.name for item in resp.items}
    assert "Delete" in names
    assert "All" in names


def test_extensions_describe(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions describe should return a JobResponse with schema."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        [
            "extensions",
            "describe",
            "--room",
            test_room,
            "@internal:selections:All",
        ],
    )
    resp = JobResponse.model_validate(data)
    assert resp.full_name == "@internal:selections:All"
    assert resp.name == "All"
    assert resp.category == "selections"
    assert isinstance(resp.schema_, dict)


def test_extensions_describe_nonexistent(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions describe with a nonexistent name should fail."""
    result = cli_runner.invoke(
        app,
        [
            "extensions",
            "describe",
            "--url",
            server_url,
            "--token",
            auth_token,
            "--room",
            test_room,
            "@internal:modifiers:DoesNotExist",
        ],
    )
    assert result.exit_code != 0


def test_extensions_run(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions run should return a TaskResponse."""
    data = invoke_cli(
        cli_runner,
        server_url,
        auth_token,
        [
            "extensions",
            "run",
            "--room",
            test_room,
            "@internal:selections:All",
        ],
    )
    resp = TaskResponse.model_validate(data)
    assert resp.job_name == "@internal:selections:All"
    assert resp.room_id == test_room


def test_extensions_run_with_args(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions run with --key value kwargs should submit them as payload."""
    result = cli_runner.invoke(
        app,
        [
            "extensions",
            "run",
            "--url",
            server_url,
            "--token",
            auth_token,
            "--room",
            test_room,
            "@internal:selections:Range",
            "--start",
            "0",
            "--stop",
            "5",
            "--step",
            "1",
        ],
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    data = json.loads(result.stdout)
    resp = TaskResponse.model_validate(data)
    assert resp.job_name == "@internal:selections:Range"
    assert resp.payload["start"] == 0
    assert resp.payload["stop"] == 5


def test_extensions_run_with_smiles(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """extensions run accepts SMILES with special chars like ( ) =."""
    result = cli_runner.invoke(
        app,
        [
            "extensions",
            "run",
            "--url",
            server_url,
            "--token",
            auth_token,
            "--room",
            test_room,
            "@internal:modifiers:AddFromSMILES",
            "--smiles",
            "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
        ],
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    data = json.loads(result.stdout)
    resp = TaskResponse.model_validate(data)
    assert resp.job_name == "@internal:modifiers:AddFromSMILES"
    assert resp.payload["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"


def test_extensions_run_wait_after_extension_name(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """--wait and --timeout placed after the extension name must be handled by Typer,
    not consumed as extension payload kwargs."""
    result = cli_runner.invoke(
        app,
        [
            "extensions",
            "run",
            "--url",
            server_url,
            "--token",
            auth_token,
            "--room",
            test_room,
            "@internal:selections:All",
            "--wait",
            "--timeout",
            "30",
        ],
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    data = json.loads(result.stdout)
    resp = TaskResponse.model_validate(data)
    assert resp.status == "completed"
    assert resp.payload == {}
