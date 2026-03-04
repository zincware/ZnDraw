"""Tests for resolve_token with TokenStore integration."""

from __future__ import annotations

from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

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


def test_explicit_token_takes_priority(token_store, stored_entry):
    """--token flag / ZNDRAW_TOKEN env should always win."""
    from zndraw.cli_agent.connection import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    with patch(
        "zndraw.cli_agent.connection._get_token_store", return_value=token_store
    ):
        result = resolve_token("http://localhost:8000", "explicit.flag.token")

    assert result == "explicit.flag.token"


def test_stored_token_used_when_valid(token_store, stored_entry):
    """Stored token should be used when GET /v1/auth/users/me returns 200."""
    from zndraw.cli_agent.connection import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_response = MagicMock()
    mock_response.status_code = 200

    with (
        patch("zndraw.cli_agent.connection._get_token_store", return_value=token_store),
        patch("zndraw.cli_agent.connection.httpx.Client") as mock_client_cls,
    ):
        mock_client = MagicMock()
        mock_client.__enter__ = MagicMock(return_value=mock_client)
        mock_client.__exit__ = MagicMock(return_value=False)
        mock_client.get.return_value = mock_response
        mock_client_cls.return_value = mock_client

        result = resolve_token("http://localhost:8000", None)

    assert result == "stored.jwt.token"


def test_stored_token_deleted_on_401(token_store, stored_entry):
    """Stored token should be removed on 401 and fall through to guest."""
    from zndraw.cli_agent.connection import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_401 = MagicMock()
    mock_401.status_code = 401

    mock_guest = MagicMock()
    mock_guest.status_code = 200
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}

    with (
        patch("zndraw.cli_agent.connection._get_token_store", return_value=token_store),
        patch("zndraw.cli_agent.connection.httpx.Client") as mock_client_cls,
    ):
        mock_client = MagicMock()
        mock_client.__enter__ = MagicMock(return_value=mock_client)
        mock_client.__exit__ = MagicMock(return_value=False)

        # First call: GET /v1/auth/users/me → 401
        # Second call (guest): POST /v1/auth/guest → 200
        mock_client.get.return_value = mock_401
        mock_client.post.return_value = mock_guest
        mock_client_cls.return_value = mock_client

        result = resolve_token("http://localhost:8000", None)

    assert result == "guest.token"
    assert token_store.get("http://localhost:8000") is None


def test_no_stored_token_falls_through_to_guest(token_store):
    """When no stored token exists, should create guest session."""
    from zndraw.cli_agent.connection import resolve_token

    mock_guest = MagicMock()
    mock_guest.status_code = 200
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}

    with (
        patch("zndraw.cli_agent.connection._get_token_store", return_value=token_store),
        patch("zndraw.cli_agent.connection.httpx.Client") as mock_client_cls,
    ):
        mock_client = MagicMock()
        mock_client.__enter__ = MagicMock(return_value=mock_client)
        mock_client.__exit__ = MagicMock(return_value=False)
        mock_client.post.return_value = mock_guest
        mock_client_cls.return_value = mock_client

        result = resolve_token("http://localhost:8000", None)

    assert result == "guest.token"
