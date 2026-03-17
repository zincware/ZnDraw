"""Tests for resolve_token with TokenStore integration."""

from __future__ import annotations

from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

import httpx
import pytest

from zndraw.server_manager import TokenEntry, TokenStore


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
    """Yield a (mock_client, patcher) for httpx.Client as context manager."""
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    with patch("zndraw.cli_agent.connection.httpx.Client", return_value=mock_client):
        yield mock_client


def test_explicit_token_takes_priority(token_store, stored_entry):
    """--token flag / ZNDRAW_TOKEN env should always win."""
    from zndraw.cli_agent.connection import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    with patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000", "explicit.flag.token")

    assert result == "explicit.flag.token"


def test_stored_token_used_when_valid(token_store, stored_entry, mock_httpx_client):
    """Stored token should be used when GET /v1/auth/users/me returns 200."""
    from zndraw.cli_agent.connection import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_httpx_client.get.return_value = mock_response

    with patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000", None)

    assert result == "stored.jwt.token"


def test_stored_token_deleted_on_401(token_store, stored_entry, mock_httpx_client):
    """Stored token should be removed on 401 and fall through to guest."""
    from zndraw.cli_agent.connection import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_401 = MagicMock()
    mock_401.status_code = 401

    mock_guest = MagicMock()
    mock_guest.status_code = 200
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}

    mock_httpx_client.get.return_value = mock_401
    mock_httpx_client.post.return_value = mock_guest

    with patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000", None)

    assert result == "guest.token"
    assert token_store.get("http://localhost:8000") is None


def test_no_stored_token_falls_through_to_guest(token_store, mock_httpx_client):
    """When no stored token exists, should create guest session."""
    from zndraw.cli_agent.connection import resolve_token

    mock_guest = MagicMock()
    mock_guest.status_code = 200
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}

    mock_httpx_client.post.return_value = mock_guest

    with patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000", None)

    assert result == "guest.token"


def test_email_password_login_success(token_store, mock_httpx_client, monkeypatch):
    """ZNDRAW_EMAIL + ZNDRAW_PASSWORD should POST to /v1/auth/jwt/login."""
    from zndraw.cli_agent.connection import resolve_token

    monkeypatch.setenv("ZNDRAW_EMAIL", "admin@example.com")
    monkeypatch.setenv("ZNDRAW_PASSWORD", "secret")

    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.json.return_value = {"access_token": "email.jwt.token"}
    mock_httpx_client.post.return_value = mock_resp

    with patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000", None)

    assert result == "email.jwt.token"
    mock_httpx_client.post.assert_called_once_with(
        "/v1/auth/jwt/login",
        data={"username": "admin@example.com", "password": "secret"},
    )


def test_email_password_login_failure_exits(
    token_store, mock_httpx_client, monkeypatch
):
    """Failed email/password login should exit, not fall through to guest."""
    from zndraw.cli_agent.connection import resolve_token

    monkeypatch.setenv("ZNDRAW_EMAIL", "admin@example.com")
    monkeypatch.setenv("ZNDRAW_PASSWORD", "wrong")

    mock_resp = MagicMock()
    mock_resp.raise_for_status.side_effect = httpx.HTTPStatusError(
        "401", request=MagicMock(), response=MagicMock()
    )
    mock_httpx_client.post.return_value = mock_resp

    with (
        patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store),
        pytest.raises(SystemExit, match="1"),
    ):
        resolve_token("http://localhost:8000", None)


@pytest.mark.parametrize(
    ("email", "password"),
    [("admin@example.com", None), (None, "secret")],
    ids=["email-only", "password-only"],
)
def test_partial_credentials_exits(token_store, monkeypatch, email, password):
    """Setting only one of ZNDRAW_EMAIL/ZNDRAW_PASSWORD should exit."""
    from zndraw.cli_agent.connection import resolve_token

    if email:
        monkeypatch.setenv("ZNDRAW_EMAIL", email)
    else:
        monkeypatch.delenv("ZNDRAW_EMAIL", raising=False)
    if password:
        monkeypatch.setenv("ZNDRAW_PASSWORD", password)
    else:
        monkeypatch.delenv("ZNDRAW_PASSWORD", raising=False)

    with (
        patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store),
        pytest.raises(SystemExit, match="1"),
    ):
        resolve_token("http://localhost:8000", None)


def test_email_password_error_does_not_leak_password(
    token_store, mock_httpx_client, monkeypatch, capsys
):
    """Error message should contain email but not password."""
    from zndraw.cli_agent.connection import resolve_token

    monkeypatch.setenv("ZNDRAW_EMAIL", "admin@example.com")
    monkeypatch.setenv("ZNDRAW_PASSWORD", "super-secret-pass")

    mock_resp = MagicMock()
    mock_resp.raise_for_status.side_effect = httpx.HTTPStatusError(
        "401", request=MagicMock(), response=MagicMock()
    )
    mock_httpx_client.post.return_value = mock_resp

    with (
        patch("zndraw.cli_agent.connection.get_token_store", return_value=token_store),
        pytest.raises(SystemExit),
    ):
        resolve_token("http://localhost:8000", None)

    stderr = capsys.readouterr().err
    assert "admin@example.com" in stderr
    assert "super-secret-pass" not in stderr
