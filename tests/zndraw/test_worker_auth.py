"""Tests for internal worker auth security."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from zndraw.config import Settings

if TYPE_CHECKING:
    from httpx import AsyncClient
    from sqlalchemy.ext.asyncio import AsyncSession


@pytest.fixture(name="settings")
def settings_fixture() -> Settings:
    """Return a Settings instance for the current test environment."""
    return Settings()


@pytest.mark.anyio
async def test_worker_login_blocked(client: AsyncClient, settings: Settings) -> None:
    """POST /v1/auth/jwt/login with the worker email must return 403."""
    resp = await client.post(
        "/v1/auth/jwt/login",
        data={
            "username": settings.internal_worker_email,
            "password": "does-not-matter",
        },
    )
    assert resp.status_code == 403


@pytest.mark.anyio
async def test_regular_login_still_works(client: AsyncClient) -> None:
    """A non-worker email should not be blocked by the guard."""
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
    client: AsyncClient, session: AsyncSession, settings: Settings
) -> None:
    """The WorkerTokenDep override (wired by conftest) must mint a valid JWT."""
    from zndraw_joblib.dependencies import get_worker_token

    app = client._transport.app  # type: ignore[union-attr]
    override = app.dependency_overrides.get(get_worker_token)
    assert override is not None, "get_worker_token override not configured"

    # Call override with the test session (mirrors FastAPI DI resolution)
    token = await override(session=session)
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
