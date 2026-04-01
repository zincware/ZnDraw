"""Tests for internal worker auth security."""

from __future__ import annotations

import pytest
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.config import Settings


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
    # Register a user first
    resp = await client.post(
        "/v1/auth/register",
        json={
            "email": "test-login@example.com",
            "password": "testpassword123",
        },
    )
    assert resp.status_code == 201

    # Login should work
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
    client: AsyncClient, settings: Settings, session: AsyncSession
) -> None:
    """The WorkerTokenDep should mint a JWT for the internal worker user."""
    from fastapi_users.authentication import JWTStrategy

    from zndraw.database import ensure_internal_worker, lookup_worker_user
    from zndraw_auth.settings import AuthSettings
    from zndraw_joblib.dependencies import get_worker_token

    # Ensure the worker user exists in the test DB
    await ensure_internal_worker(session, settings.internal_worker_email)

    auth_settings = AuthSettings()

    # Wire the override exactly as the lifespan does
    app = client._transport.app  # type: ignore[union-attr]

    async def _mint_worker_token() -> str:
        worker = await lookup_worker_user(session, settings.internal_worker_email)
        strategy = JWTStrategy(
            secret=auth_settings.secret_key.get_secret_value(),
            lifetime_seconds=auth_settings.token_lifetime_seconds,
        )
        return await strategy.write_token(worker)

    app.dependency_overrides[get_worker_token] = _mint_worker_token

    # Verify the override is configured
    override = app.dependency_overrides.get(get_worker_token)
    assert override is not None, "get_worker_token override not configured"

    token = await override()
    assert isinstance(token, str)
    assert len(token) > 0

    # Verify token is valid by calling /v1/auth/users/me
    resp = await client.get(
        "/v1/auth/users/me",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["email"] == settings.internal_worker_email
    assert data["is_superuser"] is True
