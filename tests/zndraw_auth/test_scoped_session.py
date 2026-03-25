"""Tests for current_user_scoped_session dependency."""

import pytest
from httpx import AsyncClient

from zndraw_auth import TokenResponse, UserCreate


async def _register_and_get_token(
    client: AsyncClient, email: str, password: str, login_form_class: type
) -> str:
    """Register a user, login, and return the raw access token."""
    user_data = UserCreate(email=email, password=password)
    await client.post("/auth/register", json=user_data.model_dump())
    form = login_form_class(username=email, password=password)
    response = await client.post("/auth/jwt/login", data=form.model_dump())
    token = TokenResponse.model_validate(response.json())
    return token.access_token


# --- Happy path ---


@pytest.mark.asyncio
async def test_scoped_session_returns_user(
    client: AsyncClient, login_form_class: type
) -> None:
    """Valid token returns user with accessible attributes."""
    token = await _register_and_get_token(
        client, "scoped@example.com", "password123", login_form_class
    )
    headers = {"Authorization": f"Bearer {token}"}

    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["email"] == "scoped@example.com"
    assert "user_id" in data


# --- Missing / invalid token ---


@pytest.mark.asyncio
async def test_scoped_session_no_token(client: AsyncClient) -> None:
    """Request without token returns 401."""
    response = await client.get("/test/scoped-session")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_scoped_session_invalid_token(client: AsyncClient) -> None:
    """Garbage token returns 401."""
    headers = {"Authorization": "Bearer not-a-valid-jwt"}
    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_scoped_session_expired_token(client: AsyncClient) -> None:
    """Expired JWT returns 401."""
    # A structurally valid but expired JWT (HS256, exp=0)
    import jwt as pyjwt

    expired = pyjwt.encode(
        {
            "sub": "00000000-0000-0000-0000-000000000000",
            "aud": "fastapi-users:auth",
            "exp": 0,
        },
        "test-secret-key",
        algorithm="HS256",
    )
    headers = {"Authorization": f"Bearer {expired}"}
    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 401


# --- Malformed sub claim ---


@pytest.mark.asyncio
async def test_scoped_session_malformed_sub(client: AsyncClient) -> None:
    """Token with non-UUID sub claim returns 401."""
    import jwt as pyjwt

    token = pyjwt.encode(
        {"sub": "not-a-uuid", "aud": "fastapi-users:auth"},
        "test-secret-key",
        algorithm="HS256",
    )
    headers = {"Authorization": f"Bearer {token}"}
    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_scoped_session_missing_sub(client: AsyncClient) -> None:
    """Token without sub claim returns 401."""
    import jwt as pyjwt

    token = pyjwt.encode(
        {"aud": "fastapi-users:auth"},
        "test-secret-key",
        algorithm="HS256",
    )
    headers = {"Authorization": f"Bearer {token}"}
    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 401


# --- Non-existent / inactive user ---


@pytest.mark.asyncio
async def test_scoped_session_nonexistent_user(client: AsyncClient) -> None:
    """Valid token for non-existent user returns 401."""
    import uuid

    import jwt as pyjwt

    token = pyjwt.encode(
        {"sub": str(uuid.uuid4()), "aud": "fastapi-users:auth"},
        "test-secret-key",
        algorithm="HS256",
    )
    headers = {"Authorization": f"Bearer {token}"}
    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_scoped_session_inactive_user(
    client: AsyncClient, login_form_class: type
) -> None:
    """Inactive user returns 401 even with valid token."""
    # Register and get token
    token = await _register_and_get_token(
        client, "inactive@example.com", "password123", login_form_class
    )

    # Deactivate the user via admin
    admin_form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=admin_form.model_dump())
    admin_token = TokenResponse.model_validate(response.json())
    admin_headers = {"Authorization": f"Bearer {admin_token.access_token}"}

    # Get user id
    headers = {"Authorization": f"Bearer {token}"}
    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 200
    user_id = response.json()["user_id"]

    # Deactivate via admin endpoint
    response = await client.patch(
        f"/users/{user_id}",
        headers=admin_headers,
        json={"is_active": False},
    )
    assert response.status_code == 200

    # Now the original token should be rejected
    response = await client.get("/test/scoped-session", headers=headers)
    assert response.status_code == 401
