"""Tests for CLI auth commands."""

from __future__ import annotations

from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

from zndraw.cli_agent import app

if TYPE_CHECKING:
    from typer.testing import CliRunner


def test_auth_login_opens_browser(cli_runner: CliRunner, server_url: str) -> None:
    """Login (without --code) should open the browser."""
    challenge_resp = MagicMock()
    challenge_resp.status_code = 200
    challenge_resp.raise_for_status = MagicMock()
    challenge_resp.json.return_value = {
        "code": "ABCD1234",
        "secret": "test-secret",
    }

    approved_resp = MagicMock()
    approved_resp.status_code = 200
    approved_resp.json.return_value = {
        "status": "approved",
        "token": "approved.jwt.token",
    }

    me_resp = MagicMock()
    me_resp.status_code = 200
    me_resp.json.return_value = {"id": "u1", "email": "test@example.com"}

    def mock_get(path, **kwargs):
        if "cli-login" in path:
            return approved_resp
        return me_resp

    with (
        patch("zndraw.cli_agent.auth.httpx.Client") as mock_cls,
        patch("zndraw.cli_agent.auth.webbrowser.open") as mock_browser,
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        mock_client = MagicMock()
        mock_client.__enter__ = MagicMock(return_value=mock_client)
        mock_client.__exit__ = MagicMock(return_value=False)
        mock_client.post.return_value = challenge_resp
        mock_client.get.side_effect = mock_get
        mock_cls.return_value = mock_client

        result = cli_runner.invoke(app, ["auth", "login", "--url", server_url])

    assert result.exit_code == 0, result.output
    mock_browser.assert_called_once()
    assert "Logged in" in result.output


def test_auth_login_code_flag_does_not_open_browser(
    cli_runner: CliRunner, server_url: str
) -> None:
    """Login with --code should print URL instead of opening browser."""
    challenge_resp = MagicMock()
    challenge_resp.status_code = 200
    challenge_resp.raise_for_status = MagicMock()
    challenge_resp.json.return_value = {
        "code": "ABCD1234",
        "secret": "test-secret",
    }

    approved_resp = MagicMock()
    approved_resp.status_code = 200
    approved_resp.json.return_value = {
        "status": "approved",
        "token": "approved.jwt.token",
    }

    me_resp = MagicMock()
    me_resp.status_code = 200
    me_resp.json.return_value = {"id": "u1", "email": "test@example.com"}

    def mock_get(path, **kwargs):
        if "cli-login" in path:
            return approved_resp
        return me_resp

    with (
        patch("zndraw.cli_agent.auth.httpx.Client") as mock_cls,
        patch("zndraw.cli_agent.auth.webbrowser.open") as mock_browser,
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        mock_client = MagicMock()
        mock_client.__enter__ = MagicMock(return_value=mock_client)
        mock_client.__exit__ = MagicMock(return_value=False)
        mock_client.post.return_value = challenge_resp
        mock_client.get.side_effect = mock_get
        mock_cls.return_value = mock_client

        result = cli_runner.invoke(
            app, ["auth", "login", "--url", server_url, "--code"]
        )

    assert result.exit_code == 0, result.output
    mock_browser.assert_not_called()
    assert "Visit:" in result.output


def test_auth_login_invalid_url(cli_runner: CliRunner) -> None:
    """Login to a bad URL should fail with non-zero exit code."""
    with patch("zndraw.cli_agent.auth.webbrowser.open"):
        result = cli_runner.invoke(
            app, ["auth", "login", "--url", "http://127.0.0.1:1"]
        )
    assert result.exit_code != 0
