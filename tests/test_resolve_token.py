"""Tests for auth_utils credential validation and login helpers."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
from pydantic import SecretStr

from zndraw.auth_utils import guest_login, login_with_credentials, validate_credentials


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"token": "t", "user": "u", "password": "p"}, "Cannot combine"),
        ({"token": "t", "user": None, "password": "p"}, "Cannot combine"),
        ({"token": None, "user": "u", "password": None}, "Missing --password"),
        ({"token": None, "user": None, "password": "p"}, "Missing --user"),
    ],
    ids=["token+user+password", "token+password", "user-only", "password-only"],
)
def test_invalid_credential_combinations_raise(kwargs: dict, match: str):
    with pytest.raises(ValueError, match=match):
        validate_credentials(**kwargs)


def test_valid_combinations_pass():
    validate_credentials(token="t", user=None, password=None)
    validate_credentials(token=None, user="u", password="p")
    validate_credentials(token=None, user=None, password=None)


def test_login_with_credentials_success():
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.jwt"}
    mock_resp.raise_for_status = MagicMock()
    mock_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        result = login_with_credentials(
            "http://localhost:8000", "user@test.com", "pass"
        )
    assert result == "login.jwt"


def test_login_with_secretstr_password():
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.jwt"}
    mock_resp.raise_for_status = MagicMock()
    mock_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        result = login_with_credentials(
            "http://localhost:8000", "user@test.com", SecretStr("pass")
        )
    assert result == "login.jwt"
    mock_client.post.assert_called_once_with(
        "/v1/auth/jwt/login",
        data={"username": "user@test.com", "password": "pass"},
    )


def test_guest_login_success():
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"access_token": "guest.jwt"}
    mock_resp.raise_for_status = MagicMock()
    mock_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        result = guest_login("http://localhost:8000")
    assert result == "guest.jwt"
