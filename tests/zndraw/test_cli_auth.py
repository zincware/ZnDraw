"""Tests for CLI auth commands (login, status, logout)."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

import httpx
import pytest
from typer.testing import CliRunner

from zndraw.cli_agent import app
from zndraw.state_file import StateFile, TokenEntry

runner = CliRunner()


@pytest.fixture
def state_file(tmp_path):
    return StateFile(directory=tmp_path)


@pytest.fixture
def stored_entry():
    return TokenEntry(
        access_token="stored.jwt.token",
        email="user@example.com",
        stored_at=datetime(2026, 3, 1, tzinfo=UTC),
    )


# -- auth status --------------------------------------------------------------


def test_auth_status_with_stored_token(server: str, state_file, stored_entry, monkeypatch):
    """auth status should show identity from stored token against a real server."""
    # Get a real guest token so /v1/auth/users/me will succeed
    resp = httpx.post(f"{server}/v1/auth/guest", timeout=10.0)
    resp.raise_for_status()
    real_token = resp.json()["access_token"]

    real_entry = TokenEntry(
        access_token=real_token,
        email="guest",
        stored_at=datetime(2026, 3, 1, tzinfo=UTC),
    )
    state_file.add_token(server, real_entry)

    monkeypatch.setattr("zndraw.cli_agent.auth.StateFile", lambda: state_file)
    monkeypatch.setattr("zndraw.cli_agent.auth._resolve_url", lambda url: server)

    result = runner.invoke(app, ["auth", "status"])

    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout)
    assert "server" in data


def test_auth_status_not_logged_in(server: str, state_file, monkeypatch):
    """auth status with no stored/explicit token should report not logged in."""
    monkeypatch.setattr("zndraw.cli_agent.auth.StateFile", lambda: state_file)
    monkeypatch.setattr("zndraw.cli_agent.auth._resolve_url", lambda url: server)

    result = runner.invoke(app, ["auth", "status"])

    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout)
    assert data["token_source"] == "none"
    assert data["email"] is None


# -- auth logout ---------------------------------------------------------------


def test_auth_logout_removes_token(server: str, state_file, stored_entry, monkeypatch):
    """auth logout should remove the stored token for the server."""
    state_file.add_token(server, stored_entry)

    monkeypatch.setattr("zndraw.cli_agent.auth.StateFile", lambda: state_file)
    monkeypatch.setattr("zndraw.cli_agent.auth._resolve_url", lambda url: server)

    result = runner.invoke(app, ["auth", "logout"])

    assert result.exit_code == 0
    assert state_file.get_token(server) is None


def test_auth_logout_no_token(server: str, state_file, monkeypatch):
    """auth logout when no token is stored should not error."""
    monkeypatch.setattr("zndraw.cli_agent.auth.StateFile", lambda: state_file)
    monkeypatch.setattr("zndraw.cli_agent.auth._resolve_url", lambda url: server)

    result = runner.invoke(app, ["auth", "logout"])

    assert result.exit_code == 0


# -- auth login ----------------------------------------------------------------


def test_auth_login_approved(server: str, state_file, monkeypatch):
    """auth login should store token on successful approval."""
    challenge_resp = MagicMock()
    challenge_resp.status_code = 200
    challenge_resp.raise_for_status = MagicMock()
    challenge_resp.json.return_value = {
        "code": "ABCD1234",
        "secret": "test-secret",
        "approve_url": "/auth/cli-login/ABCD1234",
    }

    # First poll: pending, second: approved
    pending_resp = MagicMock()
    pending_resp.status_code = 200
    pending_resp.json.return_value = {"status": "pending", "token": None}

    approved_resp = MagicMock()
    approved_resp.status_code = 200
    approved_resp.json.return_value = {
        "status": "approved",
        "token": "approved.jwt.token",
    }

    me_resp = MagicMock()
    me_resp.status_code = 200
    me_resp.json.return_value = {
        "id": "user-123",
        "email": "user@example.com",
        "is_superuser": False,
    }

    call_count = 0

    def mock_get(path, **kwargs):
        nonlocal call_count
        if "cli-login" in path:
            call_count += 1
            return pending_resp if call_count == 1 else approved_resp
        return me_resp

    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_client.post.return_value = challenge_resp
    mock_client.get.side_effect = mock_get

    monkeypatch.setattr("zndraw.cli_agent.auth.StateFile", lambda: state_file)
    monkeypatch.setattr("zndraw.cli_agent.auth._resolve_url", lambda url: server)

    with (
        patch("zndraw.cli_agent.auth.httpx.Client", return_value=mock_client),
        # why: webbrowser.open is a real OS side-effect that cannot run in CI
        patch("zndraw.cli_agent.auth.webbrowser.open"),
        # why: time.sleep(1) x 300 iterations would make tests take minutes
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        result = runner.invoke(app, ["auth", "login"])

    assert result.exit_code == 0
    stored = state_file.get_token(server)
    assert stored is not None
    assert stored.access_token == "approved.jwt.token"
    assert stored.email == "user@example.com"


def test_auth_login_rejected(server: str, state_file, monkeypatch):
    """auth login should show error when challenge is rejected (404)."""
    challenge_resp = MagicMock()
    challenge_resp.status_code = 200
    challenge_resp.raise_for_status = MagicMock()
    challenge_resp.json.return_value = {
        "code": "ABCD1234",
        "secret": "test-secret",
        "approve_url": "/auth/cli-login/ABCD1234",
    }

    rejected_resp = MagicMock()
    rejected_resp.status_code = 404

    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_client.post.return_value = challenge_resp
    mock_client.get.return_value = rejected_resp

    monkeypatch.setattr("zndraw.cli_agent.auth.StateFile", lambda: state_file)
    monkeypatch.setattr("zndraw.cli_agent.auth._resolve_url", lambda url: server)

    with (
        patch("zndraw.cli_agent.auth.httpx.Client", return_value=mock_client),
        # why: webbrowser.open is a real OS side-effect that cannot run in CI
        patch("zndraw.cli_agent.auth.webbrowser.open"),
        # why: time.sleep(1) x 300 iterations would make tests take minutes
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        result = runner.invoke(app, ["auth", "login"])

    assert result.exit_code != 0


def test_auth_login_expired(server: str, state_file, monkeypatch):
    """auth login should show error when challenge expires (410)."""
    challenge_resp = MagicMock()
    challenge_resp.status_code = 200
    challenge_resp.raise_for_status = MagicMock()
    challenge_resp.json.return_value = {
        "code": "ABCD1234",
        "secret": "test-secret",
        "approve_url": "/auth/cli-login/ABCD1234",
    }

    expired_resp = MagicMock()
    expired_resp.status_code = 410

    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_client.post.return_value = challenge_resp
    mock_client.get.return_value = expired_resp

    monkeypatch.setattr("zndraw.cli_agent.auth.StateFile", lambda: state_file)
    monkeypatch.setattr("zndraw.cli_agent.auth._resolve_url", lambda url: server)

    with (
        patch("zndraw.cli_agent.auth.httpx.Client", return_value=mock_client),
        # why: webbrowser.open is a real OS side-effect that cannot run in CI
        patch("zndraw.cli_agent.auth.webbrowser.open"),
        # why: time.sleep(1) x 300 iterations would make tests take minutes
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        result = runner.invoke(app, ["auth", "login"])

    assert result.exit_code != 0
