"""Tests for the users router endpoints (/users/*)."""

import uuid

import pytest
from httpx import AsyncClient

from zndraw_auth import TokenResponse, UserCreate, UserRead, UserUpdate


async def _register_and_login(
    client: AsyncClient, email: str, password: str, login_form_class: type
) -> tuple[UserRead, dict[str, str]]:
    """Helper to register, login, and return user + auth header."""
    user_data = UserCreate(email=email, password=password)
    response = await client.post("/auth/register", json=user_data.model_dump())
    assert response.status_code == 201
    user = UserRead.model_validate(response.json())

    form = login_form_class(username=email, password=password)
    response = await client.post("/auth/jwt/login", data=form.model_dump())
    token = TokenResponse.model_validate(response.json())
    headers = {"Authorization": f"Bearer {token.access_token}"}

    return user, headers


# --- GET /users/me Tests ---


@pytest.mark.asyncio
async def test_get_current_user_me(client: AsyncClient, login_form_class: type) -> None:
    """Test GET /users/me returns current user profile with email and is_superuser."""
    user, headers = await _register_and_login(
        client, "me@example.com", "password123", login_form_class
    )

    response = await client.get("/users/me", headers=headers)
    assert response.status_code == 200

    me = UserRead.model_validate(response.json())
    assert me.id == user.id
    assert me.email == "me@example.com"
    assert me.is_active is True
    assert me.is_superuser is False  # Production mode: regular users
    assert me.is_verified is False


@pytest.mark.asyncio
async def test_get_me_unauthorized(client: AsyncClient, login_form_class: type) -> None:
    """Test GET /users/me requires authentication."""
    response = await client.get("/users/me")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_get_me_admin_is_superuser(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test GET /users/me shows is_superuser=true for admin."""
    # Login as the pre-created admin
    login_form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    token = TokenResponse.model_validate(response.json())
    headers = {"Authorization": f"Bearer {token.access_token}"}

    response = await client.get("/users/me", headers=headers)
    assert response.status_code == 200

    admin = UserRead.model_validate(response.json())
    assert admin.email == "admin@test.com"
    assert admin.is_superuser is True
    assert admin.is_active is True


# --- PATCH /users/me Tests ---


@pytest.mark.asyncio
async def test_update_current_user_password(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test PATCH /users/me can update password."""
    user, headers = await _register_and_login(
        client, "update@example.com", "oldpassword", login_form_class
    )

    # Update password
    update_data = UserUpdate(password="newpassword123")
    response = await client.patch(
        "/users/me", headers=headers, json=update_data.model_dump(exclude_unset=True)
    )
    assert response.status_code == 200

    updated_user = UserRead.model_validate(response.json())
    assert updated_user.id == user.id
    assert updated_user.email == user.email

    # Verify old password doesn't work
    login_form = login_form_class(username="update@example.com", password="oldpassword")
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    assert response.status_code == 400  # Bad credentials

    # Verify new password works
    login_form = login_form_class(
        username="update@example.com", password="newpassword123"
    )
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_update_me_cannot_change_superuser(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test regular users cannot make themselves superusers."""
    user, headers = await _register_and_login(
        client, "nosuper@example.com", "password123", login_form_class
    )
    assert user.is_superuser is False

    # Try to update is_superuser (should be ignored for non-superusers)
    update_data = {"is_superuser": True}
    response = await client.patch("/users/me", headers=headers, json=update_data)
    assert response.status_code == 200

    updated_user = UserRead.model_validate(response.json())
    assert updated_user.is_superuser is False  # Should remain False


@pytest.mark.asyncio
async def test_update_me_unauthorized(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test PATCH /users/me requires authentication."""
    update_data = UserUpdate(password="newpassword")
    response = await client.patch(
        "/users/me", json=update_data.model_dump(exclude_unset=True)
    )
    assert response.status_code == 401


# --- GET /users/{user_id} Tests (Superuser Only) ---


@pytest.mark.asyncio
async def test_get_user_by_id_as_superuser(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test superuser can get any user by ID."""
    # Create a regular user
    regular_user, _ = await _register_and_login(
        client, "regular@example.com", "password123", login_form_class
    )

    # Login as admin
    login_form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    token = TokenResponse.model_validate(response.json())
    admin_headers = {"Authorization": f"Bearer {token.access_token}"}

    # Get regular user by ID
    response = await client.get(f"/users/{regular_user.id}", headers=admin_headers)
    assert response.status_code == 200

    fetched_user = UserRead.model_validate(response.json())
    assert fetched_user.id == regular_user.id
    assert fetched_user.email == "regular@example.com"
    assert fetched_user.is_superuser is False


@pytest.mark.asyncio
async def test_get_user_by_id_forbidden_for_regular_user(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test regular users cannot get other users by ID."""
    # Create two users
    _, headers1 = await _register_and_login(
        client, "user1@example.com", "password123", login_form_class
    )
    user2, _ = await _register_and_login(
        client, "user2@example.com", "password123", login_form_class
    )

    # User1 tries to get User2's profile
    response = await client.get(f"/users/{user2.id}", headers=headers1)
    assert response.status_code == 403  # Forbidden


@pytest.mark.asyncio
async def test_get_user_by_id_not_found(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test getting non-existent user returns 404."""
    # Login as admin
    login_form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    token = TokenResponse.model_validate(response.json())
    admin_headers = {"Authorization": f"Bearer {token.access_token}"}

    # Try to get non-existent user
    fake_id = uuid.uuid4()
    response = await client.get(f"/users/{fake_id}", headers=admin_headers)
    assert response.status_code == 404


# --- PATCH /users/{user_id} Tests (Superuser Only) ---


@pytest.mark.asyncio
async def test_update_user_by_id_as_superuser(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test superuser can update any user."""
    # Create a regular user
    regular_user, _ = await _register_and_login(
        client, "target@example.com", "password123", login_form_class
    )

    # Login as admin
    login_form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    token = TokenResponse.model_validate(response.json())
    admin_headers = {"Authorization": f"Bearer {token.access_token}"}

    # Update the regular user to be a superuser
    update_data = {"is_superuser": True, "is_active": True}
    response = await client.patch(
        f"/users/{regular_user.id}", headers=admin_headers, json=update_data
    )
    assert response.status_code == 200

    updated_user = UserRead.model_validate(response.json())
    assert updated_user.id == regular_user.id
    assert updated_user.is_superuser is True


@pytest.mark.asyncio
async def test_update_user_by_id_forbidden_for_regular_user(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test regular users cannot update other users."""
    # Create two users
    _, headers1 = await _register_and_login(
        client, "user1@example.com", "password123", login_form_class
    )
    user2, _ = await _register_and_login(
        client, "user2@example.com", "password123", login_form_class
    )

    # User1 tries to update User2
    update_data = {"is_superuser": True}
    response = await client.patch(
        f"/users/{user2.id}", headers=headers1, json=update_data
    )
    assert response.status_code == 403  # Forbidden


# --- DELETE /users/{user_id} Tests (Superuser Only) ---


@pytest.mark.asyncio
async def test_delete_user_by_id_as_superuser(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test superuser can delete users."""
    # Create a regular user
    regular_user, _ = await _register_and_login(
        client, "deleteme@example.com", "password123", login_form_class
    )

    # Login as admin
    login_form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    token = TokenResponse.model_validate(response.json())
    admin_headers = {"Authorization": f"Bearer {token.access_token}"}

    # Delete the user
    response = await client.delete(f"/users/{regular_user.id}", headers=admin_headers)
    assert response.status_code == 204  # No content

    # Verify user is deleted (trying to get returns 404)
    response = await client.get(f"/users/{regular_user.id}", headers=admin_headers)
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_delete_user_by_id_forbidden_for_regular_user(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test regular users cannot delete other users."""
    # Create two users
    _, headers1 = await _register_and_login(
        client, "user1@example.com", "password123", login_form_class
    )
    user2, _ = await _register_and_login(
        client, "user2@example.com", "password123", login_form_class
    )

    # User1 tries to delete User2
    response = await client.delete(f"/users/{user2.id}", headers=headers1)
    assert response.status_code == 403  # Forbidden


@pytest.mark.asyncio
async def test_delete_nonexistent_user(
    client: AsyncClient, login_form_class: type
) -> None:
    """Test deleting non-existent user returns 404."""
    # Login as admin
    login_form = login_form_class(username="admin@test.com", password="admin-password")
    response = await client.post("/auth/jwt/login", data=login_form.model_dump())
    token = TokenResponse.model_validate(response.json())
    admin_headers = {"Authorization": f"Bearer {token.access_token}"}

    # Try to delete non-existent user
    fake_id = uuid.uuid4()
    response = await client.delete(f"/users/{fake_id}", headers=admin_headers)
    assert response.status_code == 404


# --- Dev Mode Tests ---


@pytest.mark.asyncio
async def test_get_me_in_dev_mode_shows_superuser(
    client_dev_mode: AsyncClient, login_form_class: type
) -> None:
    """Test that users are superusers in dev mode."""
    _, headers = await _register_and_login(
        client_dev_mode, "devuser@example.com", "password123", login_form_class
    )

    response = await client_dev_mode.get("/users/me", headers=headers)
    assert response.status_code == 200

    me = UserRead.model_validate(response.json())
    assert me.email == "devuser@example.com"
    assert me.is_superuser is True  # Dev mode: all users are superusers
    assert me.is_active is True
