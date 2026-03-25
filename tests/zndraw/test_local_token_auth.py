"""Tests for local admin token authentication.

Tests cover:
- Local token grants admin access to the shutdown endpoint
- Wrong local token is rejected
- Normal JWT auth still works when no local token is configured
"""

from unittest.mock import patch

import pytest
import pytest_asyncio
from conftest import create_test_token, create_test_user_model
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession
from zndraw_auth import User


@pytest_asyncio.fixture
async def admin_user(session: AsyncSession) -> User:
    """Create an admin user in the database."""
    user = create_test_user_model(email="admin@local.test", is_superuser=True)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user


@pytest_asyncio.fixture
async def regular_user(session: AsyncSession) -> User:
    """Create a regular (non-admin) user in the database."""
    user = create_test_user_model(email="regular@local.test")
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user


@pytest.mark.asyncio
@patch("zndraw.routes.admin.os.kill")
async def test_local_token_grants_admin_access(mock_kill, client: AsyncClient) -> None:
    """Bearer local_token grants admin access to shutdown endpoint."""
    from zndraw.app import app

    app.state.local_token = "test-local-token"
    try:
        response = await client.post(
            "/v1/admin/shutdown",
            headers={"Authorization": "Bearer test-local-token"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "shutting_down"
    finally:
        app.state.local_token = None


@pytest.mark.asyncio
async def test_wrong_local_token_rejected(client: AsyncClient) -> None:
    """Wrong local_token is rejected with 401."""
    from zndraw.app import app

    app.state.local_token = "correct-token"
    try:
        response = await client.post(
            "/v1/admin/shutdown",
            headers={"Authorization": "Bearer wrong-token"},
        )
        assert response.status_code == 401
    finally:
        app.state.local_token = None


@pytest.mark.asyncio
@patch("zndraw.routes.admin.os.kill")
async def test_no_local_token_configured_normal_auth_works(
    mock_kill, client: AsyncClient, admin_user: User
) -> None:
    """When no local_token is configured, normal JWT admin auth still works."""
    from zndraw.app import app

    assert app.state.local_token is None

    token = create_test_token(admin_user)
    response = await client.post(
        "/v1/admin/shutdown",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "shutting_down"


@pytest.mark.asyncio
async def test_no_local_token_configured_regular_user_rejected(
    client: AsyncClient, regular_user: User
) -> None:
    """When no local_token is configured, regular users are rejected."""
    from zndraw.app import app

    assert app.state.local_token is None

    token = create_test_token(regular_user)
    response = await client.post(
        "/v1/admin/shutdown",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_no_auth_header_rejected(client: AsyncClient) -> None:
    """Requests without Authorization header are rejected."""
    from zndraw.app import app

    app.state.local_token = "test-token"
    try:
        response = await client.post("/v1/admin/shutdown")
        assert response.status_code == 401
    finally:
        app.state.local_token = None
