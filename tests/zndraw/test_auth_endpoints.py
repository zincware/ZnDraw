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
@pytest.mark.parametrize(
    ("username", "password"),
    [
        ("testuser@local.test", "wrongpassword"),
        ("nonexistent@local.test", "password"),
    ],
    ids=["invalid_password", "nonexistent_user"],
)
async def test_login_fails(
    client: AsyncClient, test_user: User, username: str, password: str
) -> None:
    """Login returns 400 for invalid credentials."""
    response = await client.post(
        "/v1/auth/jwt/login",
        data={"username": username, "password": password},
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
@pytest.mark.parametrize(
    ("setup_email", "email", "password", "expected_status"),
    [
        ("duplicate@example.com", "duplicate@example.com", "password456", 400),
        (None, "newuser@example.com", None, 422),
    ],
    ids=["duplicate_email", "missing_password"],
)
async def test_register_fails(
    client: AsyncClient,
    setup_email: str | None,
    email: str,
    password: str | None,
    expected_status: int,
) -> None:
    """Registration returns error for invalid input."""
    if setup_email:
        setup_resp = await client.post(
            "/v1/auth/register",
            json={"email": setup_email, "password": "password123"},
        )
        assert setup_resp.status_code == 201
    body = {"email": email}
    if password is not None:
        body["password"] = password
    response = await client.post("/v1/auth/register", json=body)
    assert response.status_code == expected_status
