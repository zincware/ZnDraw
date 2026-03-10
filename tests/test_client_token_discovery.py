"""Tests for ZnDraw Python client token auto-discovery."""

from __future__ import annotations

from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

from zndraw.server_manager import TokenEntry, TokenStore


def test_try_stored_token_valid(tmp_path):
    """_try_stored_token should use stored token when valid."""
    store = TokenStore(directory=tmp_path)
    store.set(
        "http://localhost:8000",
        TokenEntry(
            access_token="stored.jwt",
            email="user@test.com",
            stored_at=datetime(2026, 3, 1, tzinfo=UTC),
        ),
    )

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {"email": "user@test.com"}

    from zndraw.client import ZnDraw

    instance = object.__new__(ZnDraw)
    instance.url = "http://localhost:8000"
    instance.api = MagicMock()

    with (
        patch("zndraw.server_manager.TokenStore", return_value=store),
        patch("httpx.get", return_value=mock_response),
    ):
        result = instance._try_stored_token()

    assert result is True
    assert instance.api.token == "stored.jwt"
    assert instance.user == "user@test.com"


def test_try_stored_token_invalid(tmp_path):
    """_try_stored_token should delete invalid token and return False."""
    store = TokenStore(directory=tmp_path)
    store.set(
        "http://localhost:8000",
        TokenEntry(
            access_token="expired.jwt",
            email="user@test.com",
            stored_at=datetime(2026, 3, 1, tzinfo=UTC),
        ),
    )

    mock_response = MagicMock()
    mock_response.status_code = 401

    from zndraw.client import ZnDraw

    instance = object.__new__(ZnDraw)
    instance.url = "http://localhost:8000"
    instance.api = MagicMock()

    with (
        patch("zndraw.server_manager.TokenStore", return_value=store),
        patch("httpx.get", return_value=mock_response),
    ):
        result = instance._try_stored_token()

    assert result is False
    assert store.get("http://localhost:8000") is None


def test_try_stored_token_no_entry(tmp_path):
    """_try_stored_token should return False when no stored token exists."""
    store = TokenStore(directory=tmp_path)

    from zndraw.client import ZnDraw

    instance = object.__new__(ZnDraw)
    instance.url = "http://localhost:8000"
    instance.api = MagicMock()

    with patch("zndraw.server_manager.TokenStore", return_value=store):
        result = instance._try_stored_token()

    assert result is False
