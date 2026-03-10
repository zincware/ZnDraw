"""Tests for CLI auth commands (login, status, logout)."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from zndraw.cli_agent import app
from zndraw.server_manager import TokenEntry, TokenStore

runner = CliRunner()


@pytest.fixture
def token_store(tmp_path):
    return TokenStore(directory=tmp_path)


@pytest.fixture
def stored_entry():
    return TokenEntry(
        access_token="stored.jwt.token",
        email="user@example.com",
        stored_at=datetime(2026, 3, 1, tzinfo=UTC),
    )


@pytest.fixture
def mock_httpx_client():
    """Yield a mock httpx.Client context manager for auth module."""
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    with patch("zndraw.cli_agent.auth.httpx.Client", return_value=mock_client):
        yield mock_client


# ── auth status ──────────────────────────────────────────────────────


def test_auth_status_with_stored_token(token_store, stored_entry, mock_httpx_client):
    """auth status should show identity from stored token."""
    token_store.set("http://localhost:8000", stored_entry)

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "id": "abc-123",
        "email": "user@example.com",
        "is_superuser": False,
    }
    mock_httpx_client.get.return_value = mock_response

    with (
        patch("zndraw.cli_agent.auth.get_token_store", return_value=token_store),
        patch(
            "zndraw.cli_agent.auth.resolve_url", return_value="http://localhost:8000"
        ),
    ):
        result = runner.invoke(app, ["auth", "status"])

    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout)
    assert data["email"] == "user@example.com"
    assert "server" in data


def test_auth_status_not_logged_in(token_store):
    """auth status with no stored/explicit token should report not logged in."""
    with (
        patch("zndraw.cli_agent.auth.get_token_store", return_value=token_store),
        patch(
            "zndraw.cli_agent.auth.resolve_url", return_value="http://localhost:8000"
        ),
    ):
        result = runner.invoke(app, ["auth", "status"])

    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout)
    assert data["token_source"] == "none"
    assert data["email"] is None


# ── auth logout ──────────────────────────────────────────────────────


def test_auth_logout_removes_token(token_store, stored_entry):
    """auth logout should remove the stored token for the server."""
    token_store.set("http://localhost:8000", stored_entry)

    with (
        patch("zndraw.cli_agent.auth.get_token_store", return_value=token_store),
        patch(
            "zndraw.cli_agent.auth.resolve_url", return_value="http://localhost:8000"
        ),
    ):
        result = runner.invoke(app, ["auth", "logout"])

    assert result.exit_code == 0
    assert token_store.get("http://localhost:8000") is None


def test_auth_logout_no_token(token_store):
    """auth logout when no token is stored should not error."""
    with (
        patch("zndraw.cli_agent.auth.get_token_store", return_value=token_store),
        patch(
            "zndraw.cli_agent.auth.resolve_url", return_value="http://localhost:8000"
        ),
    ):
        result = runner.invoke(app, ["auth", "logout"])

    assert result.exit_code == 0


# ── auth login ───────────────────────────────────────────────────────


def test_auth_login_approved(token_store, mock_httpx_client):
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

    mock_httpx_client.post.return_value = challenge_resp
    mock_httpx_client.get.side_effect = mock_get

    with (
        patch("zndraw.cli_agent.auth.get_token_store", return_value=token_store),
        patch(
            "zndraw.cli_agent.auth.resolve_url", return_value="http://localhost:8000"
        ),
        patch("zndraw.cli_agent.auth.webbrowser.open"),
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        result = runner.invoke(app, ["auth", "login"])

    assert result.exit_code == 0
    stored = token_store.get("http://localhost:8000")
    assert stored is not None
    assert stored.access_token == "approved.jwt.token"
    assert stored.email == "user@example.com"


def test_auth_login_rejected(token_store, mock_httpx_client):
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

    mock_httpx_client.post.return_value = challenge_resp
    mock_httpx_client.get.return_value = rejected_resp

    with (
        patch("zndraw.cli_agent.auth.get_token_store", return_value=token_store),
        patch(
            "zndraw.cli_agent.auth.resolve_url", return_value="http://localhost:8000"
        ),
        patch("zndraw.cli_agent.auth.webbrowser.open"),
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        result = runner.invoke(app, ["auth", "login"])

    assert result.exit_code != 0


def test_auth_login_expired(token_store, mock_httpx_client):
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

    mock_httpx_client.post.return_value = challenge_resp
    mock_httpx_client.get.return_value = expired_resp

    with (
        patch("zndraw.cli_agent.auth.get_token_store", return_value=token_store),
        patch(
            "zndraw.cli_agent.auth.resolve_url", return_value="http://localhost:8000"
        ),
        patch("zndraw.cli_agent.auth.webbrowser.open"),
        patch("zndraw.cli_agent.auth.time.sleep"),
    ):
        result = runner.invoke(app, ["auth", "login"])

    assert result.exit_code != 0
