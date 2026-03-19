"""Shared token resolution for CLI and Python client.

Single source of truth for the authentication fallback chain:
1. Explicit credentials (--token OR --user/--password)
2. Stored token from ~/.zndraw/tokens.json
3. Guest session fallback
"""

from __future__ import annotations

import os
import warnings

import httpx
from pydantic import SecretStr

from zndraw.server_manager import TokenStore


def get_token_store() -> TokenStore:
    """Return the default TokenStore (testable seam)."""
    return TokenStore()


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


def resolve_token(
    base_url: str,
    token: str | None = None,
    user: str | None = None,
    password: SecretStr | str | None = None,
) -> str:
    """Resolve an auth token from explicit credentials, stored token, or guest.

    Resolution order (explicit always wins over implicit):

    1. Explicit: ``token`` returned as-is, OR ``user``+``password``
       login via ``POST /v1/auth/jwt/login`` (token NOT stored).
    2. Stored token from ``~/.zndraw/tokens.json`` (validated, deleted on 401).
    3. Guest fallback via ``POST /v1/auth/guest``.

    Parameters
    ----------
    base_url
        Server URL for network calls.
    token
        Explicit JWT token (from ``--token`` / ``ZNDRAW_TOKEN``).
    user
        User email (from ``--user`` / ``ZNDRAW_USER``).
    password
        Password (from ``--password`` / ``ZNDRAW_PASSWORD``).
        Accepts ``SecretStr`` or plain ``str``.

    Raises
    ------
    ValueError
        On invalid credential combinations.
    httpx.HTTPStatusError
        On failed login or guest session creation.
    """
    # --- Migration: bridge ZNDRAW_EMAIL → user ---
    legacy_user = os.environ.get("ZNDRAW_EMAIL")
    if (
        token is None
        and user is None
        and legacy_user
        and not os.environ.get("ZNDRAW_USER")
    ):
        warnings.warn(
            "ZNDRAW_EMAIL is deprecated, use ZNDRAW_USER instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        user = legacy_user

    # --- Validation (fail fast) ---
    validate_credentials(token, user, password)

    # Unwrap SecretStr
    raw_password: str | None = None
    if isinstance(password, SecretStr):
        raw_password = password.get_secret_value()
    elif isinstance(password, str):
        raw_password = password

    # --- Tier 1: Explicit credentials ---
    if token is not None:
        return token

    if user is not None and raw_password is not None:
        with httpx.Client(base_url=base_url, timeout=10.0) as client:
            resp = client.post(
                "/v1/auth/jwt/login",
                data={"username": user, "password": raw_password},
            )
            resp.raise_for_status()
            return resp.json()["access_token"]

    # --- Tier 2: Stored token ---
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

    # --- Tier 3: Guest fallback ---
    with httpx.Client(base_url=base_url, timeout=10.0) as client:
        resp = client.post("/v1/auth/guest")
        resp.raise_for_status()
        return resp.json()["access_token"]
