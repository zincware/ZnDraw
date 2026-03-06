"""Tests for CLI auth commands."""

from __future__ import annotations

import json

from typer.testing import CliRunner

from zndraw.cli_agent import app


def test_auth_login_returns_access_token(
    cli_runner: CliRunner, server_url: str
) -> None:
    """Guest login should return an access_token."""
    result = cli_runner.invoke(app, ["--url", server_url, "auth", "login"])
    assert result.exit_code == 0, result.stderr
    data = json.loads(result.stdout)
    assert "access_token" in data
    assert isinstance(data["access_token"], str)
    assert len(data["access_token"]) > 0


def test_auth_login_tokens_are_unique(cli_runner: CliRunner, server_url: str) -> None:
    """Each guest login should produce a different token."""
    result1 = cli_runner.invoke(app, ["--url", server_url, "auth", "login"])
    result2 = cli_runner.invoke(app, ["--url", server_url, "auth", "login"])
    assert result1.exit_code == 0
    assert result2.exit_code == 0
    token1 = json.loads(result1.stdout)["access_token"]
    token2 = json.loads(result2.stdout)["access_token"]
    assert token1 != token2


def test_auth_login_invalid_url(cli_runner: CliRunner) -> None:
    """Login to a bad URL should fail with non-zero exit code."""
    result = cli_runner.invoke(app, ["--url", "http://127.0.0.1:1", "auth", "login"])
    assert result.exit_code != 0
