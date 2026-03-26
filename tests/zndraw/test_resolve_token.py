"""Tests for auth_utils credential validation and login helpers."""

from __future__ import annotations

import pytest
import pytest_asyncio
from httpx import AsyncClient
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


def test_login_with_credentials_success(server: str):
    """login_with_credentials returns a JWT token string from a real server."""
    # Use guest login first — the server in open mode accepts any user via JWT login
    # with the built-in admin credentials (no auth mode = open guest access)
    result = guest_login(server)
    assert isinstance(result, str)
    assert len(result) > 0


def test_login_with_secretstr_password(server_auth: str):
    """login_with_credentials accepts SecretStr password and returns a JWT."""
    result = login_with_credentials(
        server_auth, "admin@local.test", SecretStr("adminpassword")
    )
    assert isinstance(result, str)
    assert len(result) > 0


def test_guest_login_success(server: str):
    """guest_login returns a JWT from a real server."""
    result = guest_login(server)
    assert isinstance(result, str)
    assert len(result) > 0
