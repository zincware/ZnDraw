"""Authentication utilities for ZnDraw client.

Provides credential validation and auth helper functions.
Token resolution chain (stored token, local_token) is handled by
the pydantic-settings source chain (StateFileSource).
"""

from __future__ import annotations

import httpx
from pydantic import SecretStr


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


# ---------------------------------------------------------------------------
# Backwards-compatible shims (used by cli_agent — removed in Task 8/9)
# ---------------------------------------------------------------------------


def get_token_store():
    """Return the default TokenStore (testable seam).

    .. deprecated::
        Used by cli_agent code; will be removed when CLI commands
        are refactored to use ClientSettings.
    """
    from zndraw.server_manager import TokenStore

    return TokenStore()


def resolve_token(
    base_url: str,
    token: str | None = None,
    user: str | None = None,
    password: SecretStr | str | None = None,
) -> str:
    """Resolve an auth token from explicit credentials, stored token, or guest.

    .. deprecated::
        Used by cli_agent code; will be removed when CLI commands
        are refactored to use ClientSettings.
    """
    validate_credentials(token, user, password)

    raw_password: str | None = None
    if isinstance(password, SecretStr):
        raw_password = password.get_secret_value()
    elif isinstance(password, str):
        raw_password = password

    if token is not None:
        return token

    if user is not None and raw_password is not None:
        return login_with_credentials(base_url, user, raw_password)

    store = get_token_store()
    entry = store.get(base_url)
    if entry is not None:
        with httpx.Client(
            base_url=base_url,
            headers={"Authorization": f"Bearer {entry.access_token}"},
            timeout=10.0,
        ) as client:
            resp = client.get("/v1/auth/users/me")
            if resp.status_code == 200:
                return entry.access_token
            if resp.status_code in {401, 403}:
                store.delete(base_url)
            else:
                resp.raise_for_status()

    return guest_login(base_url)
