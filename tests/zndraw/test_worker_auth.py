"""Tests for internal worker auth security."""

from __future__ import annotations

from types import SimpleNamespace
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from httpx import AsyncClient

    from zndraw.config import Settings


@pytest.mark.anyio
async def test_regular_login_still_works(client: AsyncClient) -> None:
    """A non-worker email should not be blocked by the default login route."""
    resp = await client.post(
        "/v1/auth/register",
        json={
            "email": "test-login@example.com",
            "password": "testpassword123",
        },
    )
    assert resp.status_code == 201

    resp = await client.post(
        "/v1/auth/jwt/login",
        data={
            "username": "test-login@example.com",
            "password": "testpassword123",
        },
    )
    assert resp.status_code == 200
    assert "access_token" in resp.json()


@pytest.mark.anyio
async def test_worker_token_dep_mints_valid_jwt(
    client: AsyncClient, settings: Settings
) -> None:
    """The default get_worker_token dependency must mint a valid JWT.

    After the deadlock fix, get_worker_token reads the cached worker user
    from ``app.state.internal_worker_user`` and takes no session arg.
    """
    from zndraw_joblib.dependencies import get_worker_token

    app = client._transport.app  # type: ignore[union-attr]
    request = SimpleNamespace(app=app)

    token = await get_worker_token(request)  # type: ignore[arg-type]
    assert isinstance(token, str)
    assert len(token) > 0

    # Verify token authenticates as the internal worker
    resp = await client.get(
        "/v1/auth/users/me",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["email"] == settings.internal_worker_email
    assert data["is_superuser"] is True
