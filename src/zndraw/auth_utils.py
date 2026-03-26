"""Authentication utilities for ZnDraw client.

Provides credential validation and auth helper functions.
Token resolution chain (stored token, local_token) is handled by
the pydantic-settings source chain (StateFileSource).
"""

from __future__ import annotations

import logging

import httpx
from pydantic import SecretStr

log = logging.getLogger(__name__)


def validate_credentials(
    token: str | None,
    user: str | None,
    password: SecretStr | str | None,
) -> None:
    """Validate credential combinations (fail fast, before any network call).

    Raises
    ------
    ValueError
        On invalid credential combinations.
    """
    if token is not None and (user is not None or password is not None):
        msg = "Cannot combine --token with --user/--password"
        raise ValueError(msg)
    if user is not None and password is None:
        msg = "Missing --password (required when --user is provided)"
        raise ValueError(msg)
    if password is not None and user is None:
        msg = "Missing --user (required when --password is provided)"
        raise ValueError(msg)


def login_with_credentials(
    base_url: str,
    user: str,
    password: SecretStr | str,
) -> str:
    """Login with email+password, return access token.

    Parameters
    ----------
    base_url
        Server URL.
    user
        User email.
    password
        User password.

    Returns
    -------
    str
        JWT access token.
    """
    raw = password.get_secret_value() if isinstance(password, SecretStr) else password
    with httpx.Client(base_url=base_url, timeout=10.0) as client:
        resp = client.post(
            "/v1/auth/jwt/login",
            data={"username": user, "password": raw},
        )
        resp.raise_for_status()
        return resp.json()["access_token"]


def resolve_or_refresh_token(base_url: str, token: str) -> str:
    """Validate a stored token and fall back to guest login on 401.

    Parameters
    ----------
    base_url
        Server URL.
    token
        Previously stored JWT to validate.

    Returns
    -------
    str
        The original token if still valid, or a fresh guest token.
    """
    try:
        with httpx.Client(base_url=base_url, timeout=10.0) as client:
            resp = client.get(
                "/v1/auth/users/me",
                headers={"Authorization": f"Bearer {token}"},
            )
            if resp.status_code == 401:
                log.debug("Stored token expired or invalid, evicting and refreshing")
                from zndraw.state_file import StateFile

                state = StateFile()
                state.remove_token(base_url)
                return guest_login(base_url)
            resp.raise_for_status()
            return token
    except httpx.RequestError:
        # Network error — keep token, let caller handle
        return token


def guest_login(base_url: str) -> str:
    """Create a guest session, return access token.

    Parameters
    ----------
    base_url
        Server URL.

    Returns
    -------
    str
        Guest JWT access token.
    """
    with httpx.Client(base_url=base_url, timeout=10.0) as client:
        resp = client.post("/v1/auth/guest")
        resp.raise_for_status()
        return resp.json()["access_token"]
