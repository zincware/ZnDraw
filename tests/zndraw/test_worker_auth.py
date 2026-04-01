"""Tests for internal worker auth security."""

import pytest
from httpx import AsyncClient

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
