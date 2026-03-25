"""Tests for authentication flows."""

import pytest
from httpx import AsyncClient

from zndraw_auth import TokenResponse, UserCreate, UserRead


@pytest.mark.asyncio
async def test_register_user(client: AsyncClient) -> None:
    """Test user registration."""
    user_data = UserCreate(email="test@example.com", password="testpassword123")
    response = await client.post(
        "/auth/register",
        json=user_data.model_dump(),
    )
    assert response.status_code == 201
    user = UserRead.model_validate(response.json())
    assert user.email == "test@example.com"
    assert user.id is not None
    assert user.is_active is True
    assert user.is_superuser is False


@pytest.mark.asyncio
async def test_login_user(client: AsyncClient, login_form_class: type) -> None:
    """Test user login."""
    # Register first
    user_data = UserCreate(email="login@example.com", password="testpassword123")
    await client.post("/auth/register", json=user_data.model_dump())

    # Login
    form = login_form_class(username=user_data.email, password="testpassword123")
    response = await client.post("/auth/jwt/login", data=form.model_dump())
    assert response.status_code == 200
    token = TokenResponse.model_validate(response.json())
    assert token.access_token
    assert token.token_type == "bearer"


@pytest.mark.asyncio
async def test_login_invalid_credentials(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test login with wrong password."""
    # Register
    user_data = UserCreate(email="wrong@example.com", password="correctpassword")
    await client.post("/auth/register", json=user_data.model_dump())

    # Login with wrong password
    form = login_form_class(username=user_data.email, password="wrongpassword")
    response = await client.post("/auth/jwt/login", data=form.model_dump())
    assert response.status_code == 400


@pytest.mark.asyncio
async def test_register_duplicate_email(client: AsyncClient) -> None:
    """Test registering with existing email."""
    # Register first user
    user_data = UserCreate(email="dupe@example.com", password="password123")
    await client.post("/auth/register", json=user_data.model_dump())

    # Try to register again with same email
    duplicate_data = UserCreate(email="dupe@example.com", password="password456")
    response = await client.post("/auth/register", json=duplicate_data.model_dump())
    assert response.status_code == 400


# --- Tests for exported dependency injections ---


async def _get_auth_header(
    client: AsyncClient, email: str, password: str, login_form_class: type
) -> dict[str, str]:
    """Helper to register, login, and return auth header."""
    user_data = UserCreate(email=email, password=password)
    await client.post("/auth/register", json=user_data.model_dump())
    form = login_form_class(username=email, password=password)
    response = await client.post("/auth/jwt/login", data=form.model_dump())
    token = TokenResponse.model_validate(response.json())
    return {"Authorization": f"Bearer {token.access_token}"}


@pytest.mark.asyncio
async def test_current_active_user_dependency(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test current_active_user dependency injection."""
    headers = await _get_auth_header(
        client, "active@example.com", "password123", login_form_class
    )

    response = await client.get("/test/protected", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["email"] == "active@example.com"
    assert "user_id" in data


@pytest.mark.asyncio
async def test_current_active_user_unauthorized(client: AsyncClient) -> None:
    """Test current_active_user rejects unauthenticated requests."""
    response = await client.get("/test/protected")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_current_superuser_forbidden(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test current_superuser rejects non-superusers."""
    headers = await _get_auth_header(
        client, "regular@example.com", "password123", login_form_class
    )

    response = await client.get("/test/superuser", headers=headers)
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_current_optional_user_authenticated(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test current_optional_user with authenticated user."""
    headers = await _get_auth_header(
        client, "optional@example.com", "password123", login_form_class
    )

    response = await client.get("/test/optional", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["authenticated"] == "true"
    assert data["user_id"] is not None


@pytest.mark.asyncio
async def test_current_optional_user_anonymous(client: AsyncClient) -> None:
    """Test current_optional_user without authentication."""
    response = await client.get("/test/optional")
    assert response.status_code == 200
    data = response.json()
    assert data["authenticated"] == "false"
    assert data["user_id"] is None


@pytest.mark.asyncio
async def test_get_session_dependency(client: AsyncClient) -> None:
    """Test get_session dependency injection."""
    response = await client.get("/test/session")
    assert response.status_code == 200
    data = response.json()
    assert data["db_check"] == "1"


# --- Tests for dev mode and default admin ---


@pytest.mark.asyncio
async def test_default_admin_created_on_startup(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test that default admin is created on startup when configured."""
    # The test_settings fixture has default_admin_email="admin@test.com"
    # This admin should have been created during ensure_default_admin()
    form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=form.model_dump())
    assert response.status_code == 200
    token = TokenResponse.model_validate(response.json())

    # Verify admin can access superuser route
    headers = {"Authorization": f"Bearer {token.access_token}"}
    response = await client.get("/test/superuser", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["is_superuser"] == "True"


@pytest.mark.asyncio
async def test_dev_mode_all_users_superuser(
    client_dev_mode: AsyncClient, login_form_class: type
) -> None:
    """Test that users are superusers in dev mode (no admin configured)."""
    # Register a user
    user_data = UserCreate(email="devuser@example.com", password="password123")
    response = await client_dev_mode.post(
        "/auth/register",
        json=user_data.model_dump(),
    )
    assert response.status_code == 201
    user = UserRead.model_validate(response.json())
    assert user.is_superuser is True

    # Login and verify we can access superuser-protected route
    form = login_form_class(username=user_data.email, password="password123")
    response = await client_dev_mode.post("/auth/jwt/login", data=form.model_dump())
    token = TokenResponse.model_validate(response.json())
    headers = {"Authorization": f"Bearer {token.access_token}"}

    response = await client_dev_mode.get("/test/superuser", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["is_superuser"] == "True"
