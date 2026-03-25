"""Tests for authentication endpoints (zndraw-auth / fastapi-users)."""

import pytest
from httpx import AsyncClient
from zndraw_auth import User

# =============================================================================
# Login Tests
# =============================================================================


@pytest.mark.asyncio
async def test_login(client: AsyncClient, test_user: User) -> None:
    """Test login with valid credentials returns JWT."""
    response = await client.post(
        "/v1/auth/jwt/login",
        data={"username": "testuser@local.test", "password": "testpassword"},
    )
    assert response.status_code == 200

    body = response.json()
    assert body["token_type"] == "bearer"
    assert body["access_token"]


@pytest.mark.asyncio
async def test_login_invalid_password(client: AsyncClient, test_user: User) -> None:
    """Test that invalid password returns 400."""
    response = await client.post(
        "/v1/auth/jwt/login",
        data={"username": "testuser@local.test", "password": "wrongpassword"},
    )
    assert response.status_code == 400


@pytest.mark.asyncio
async def test_login_nonexistent_user(client: AsyncClient) -> None:
    """Test that nonexistent user returns 400."""
    response = await client.post(
        "/v1/auth/jwt/login",
        data={"username": "nonexistent@local.test", "password": "password"},
    )
    assert response.status_code == 400


# =============================================================================
# Guest Session Tests
# =============================================================================


@pytest.mark.asyncio
async def test_guest_session(client: AsyncClient) -> None:
    """Test creating a guest session returns JWT."""
    response = await client.post("/v1/auth/guest")
    assert response.status_code == 200

    body = response.json()
    assert body["token_type"] == "bearer"
    assert body["access_token"]


@pytest.mark.asyncio
async def test_guest_sessions_unique_tokens(client: AsyncClient) -> None:
    """Test that multiple guest sessions get unique tokens."""
    response1 = await client.post("/v1/auth/guest")
    response2 = await client.post("/v1/auth/guest")

    assert response1.status_code == 200
    assert response2.status_code == 200

    assert response1.json()["access_token"] != response2.json()["access_token"]


# =============================================================================
# Registration Tests
# =============================================================================


@pytest.mark.asyncio
async def test_register(client: AsyncClient) -> None:
    """Test registering a new user."""
    response = await client.post(
        "/v1/auth/register",
        json={"email": "newuser@example.com", "password": "newpassword"},
    )
    assert response.status_code == 201

    body = response.json()
    assert body["email"] == "newuser@example.com"
    assert "id" in body


@pytest.mark.asyncio
async def test_register_duplicate_email(client: AsyncClient) -> None:
    """Test that duplicate email returns 400."""
    # First registration succeeds
    await client.post(
        "/v1/auth/register",
        json={"email": "duplicate@example.com", "password": "password123"},
    )
    # Second registration with same email fails
    response = await client.post(
        "/v1/auth/register",
        json={"email": "duplicate@example.com", "password": "password456"},
    )
    assert response.status_code == 400


@pytest.mark.asyncio
async def test_register_missing_fields(client: AsyncClient) -> None:
    """Test that missing required fields returns 422."""
    response = await client.post(
        "/v1/auth/register",
        json={"email": "newuser@example.com"},
    )
    assert response.status_code == 422
