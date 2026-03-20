"""Tests for shared resolve_token utility in auth_utils."""

from __future__ import annotations

from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

import httpx
import pytest
from pydantic import SecretStr

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
    """Yield a mock for httpx.Client as context manager."""
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        yield mock_client


# --- Validation tests ---


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"token": "t", "user": "u", "password": "p"}, "Cannot combine"),
        ({"token": "t", "password": "p"}, "Cannot combine"),
        ({"user": "u"}, "Missing --password"),
        ({"password": "p"}, "Missing --user"),
    ],
    ids=["token+user+password", "token+password", "user-only", "password-only"],
)
def test_invalid_credential_combinations_raise(kwargs: dict, match: str):
    """Invalid credential combinations should fail fast."""
    from zndraw.auth_utils import resolve_token

    with pytest.raises(ValueError, match=match):
        resolve_token("http://localhost:8000", **kwargs)


def test_secretstr_password_accepted(mock_httpx_client):
    """SecretStr password should be unwrapped automatically."""
    from zndraw.auth_utils import resolve_token

    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.token"}
    mock_httpx_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.get_token_store"):
        result = resolve_token(
            "http://localhost:8000",
            user="admin@example.com",
            password=SecretStr("secret"),
        )

    assert result == "login.token"
    mock_httpx_client.post.assert_called_once_with(
        "/v1/auth/jwt/login",
        data={"username": "admin@example.com", "password": "secret"},
    )


# --- Resolution chain tests ---


def test_explicit_token_takes_priority(token_store, stored_entry):
    """--token flag should always win over stored token."""
    from zndraw.auth_utils import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000", token="explicit.flag.token")

    assert result == "explicit.flag.token"


def test_user_password_login_success(token_store, mock_httpx_client):
    """--user/--password should POST to /v1/auth/jwt/login."""
    from zndraw.auth_utils import resolve_token

    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.jwt.token"}
    mock_httpx_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token(
            "http://localhost:8000", user="admin@example.com", password="secret"
        )

    assert result == "login.jwt.token"
    mock_httpx_client.post.assert_called_once_with(
        "/v1/auth/jwt/login",
        data={"username": "admin@example.com", "password": "secret"},
    )


def test_stored_token_used_when_valid(token_store, stored_entry, mock_httpx_client):
    """Stored token should be used when GET /v1/auth/users/me returns 200."""
    from zndraw.auth_utils import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_httpx_client.get.return_value = mock_response

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000")

    assert result == "stored.jwt.token"


def test_stored_token_deleted_on_401(token_store, stored_entry, mock_httpx_client):
    """Stored token should be removed on 401 and fall through to guest."""
    from zndraw.auth_utils import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_401 = MagicMock()
    mock_401.status_code = 401

    mock_guest = MagicMock()
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}

    mock_httpx_client.get.return_value = mock_401
    mock_httpx_client.post.return_value = mock_guest

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000")

    assert result == "guest.token"
    assert token_store.get("http://localhost:8000") is None


def test_stored_token_preserved_on_5xx(token_store, stored_entry, mock_httpx_client):
    """Transient 5xx should not delete stored token — raise instead."""
    from zndraw.auth_utils import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_500 = MagicMock()
    mock_500.status_code = 500
    mock_500.raise_for_status.side_effect = httpx.HTTPStatusError(
        "500", request=MagicMock(), response=mock_500
    )
    mock_httpx_client.get.return_value = mock_500

    with (
        patch("zndraw.auth_utils.get_token_store", return_value=token_store),
        pytest.raises(httpx.HTTPStatusError),
    ):
        resolve_token("http://localhost:8000")

    assert token_store.get("http://localhost:8000") is not None


def test_no_stored_token_falls_through_to_guest(token_store, mock_httpx_client):
    """When no stored token exists, should create guest session."""
    from zndraw.auth_utils import resolve_token

    mock_guest = MagicMock()
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}
    mock_httpx_client.post.return_value = mock_guest

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000")

    assert result == "guest.token"


def test_user_password_login_failure_raises(mock_httpx_client):
    """Failed login should raise, not fall through to guest."""
    from zndraw.auth_utils import resolve_token

    mock_resp = MagicMock()
    mock_resp.raise_for_status.side_effect = httpx.HTTPStatusError(
        "401", request=MagicMock(), response=MagicMock()
    )
    mock_httpx_client.post.return_value = mock_resp

    with (
        patch("zndraw.auth_utils.get_token_store"),
        pytest.raises(httpx.HTTPStatusError),
    ):
        resolve_token(
            "http://localhost:8000", user="admin@example.com", password="wrong"
        )
