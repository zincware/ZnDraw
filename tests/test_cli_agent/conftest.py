"""Fixtures for CLI agent tests.

Uses the existing server_factory from the main conftest to get a real
running server, then provides a CliRunner and helpers for invoking
zndraw-cli commands.
"""

from __future__ import annotations

import json
from typing import Any

import httpx
import pytest
from typer.testing import CliRunner

from zndraw.cli_agent import app


@pytest.fixture(name="cli_runner")
def cli_runner_fixture() -> CliRunner:
    """Provide a Typer CliRunner."""
    return CliRunner()


@pytest.fixture(name="server_url")
def server_url_fixture(server: str) -> str:
    """Provide the test server URL."""
    return server


@pytest.fixture(name="auth_token")
def auth_token_fixture(server_url: str) -> str:
    """Get a guest auth token via the REST guest endpoint."""
    resp = httpx.post(f"{server_url}/v1/auth/guest", timeout=10.0)
    resp.raise_for_status()
    return resp.json()["access_token"]


@pytest.fixture(name="test_room")
def test_room_fixture(server_url: str, cli_runner: CliRunner, auth_token: str) -> str:
    """Create a test room with one empty frame and return its ID."""
    result = cli_runner.invoke(
        app,
        ["rooms", "create", "--url", server_url, "--token", auth_token],
    )
    assert result.exit_code == 0, result.stderr
    data = json.loads(result.stdout)
    return data["room_id"]


def invoke_cli(
    cli_runner: CliRunner,
    server_url: str,
    auth_token: str,
    args: list[str],
) -> Any:
    """Invoke a CLI command and return parsed JSON output.

    Parameters
    ----------
    cli_runner
        The Typer CliRunner.
    server_url
        Server URL for --url option.
    auth_token
        JWT for --token option.
    args
        Command args: ``["subapp", "command", ...]``.
        Connection options (--url, --token) are inserted after the
        command path (first two elements).

    Returns
    -------
    Any
        Parsed JSON from stdout.
    """
    cmd_path = args[:2]
    cmd_args = args[2:]
    connection_opts = ["--url", server_url, "--token", auth_token]
    result = cli_runner.invoke(
        app,
        [*cmd_path, *connection_opts, *cmd_args],
    )
    assert result.exit_code == 0, (
        f"exit={result.exit_code}\nstderr={result.stderr}\nstdout={result.stdout}"
    )
    return json.loads(result.stdout)
