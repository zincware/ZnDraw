"""Tests for admin functionality.

Tests cover:
- User management endpoints (list, get, update, delete)
- Server settings endpoints (get, set, unset default room)
- Room creation with @-prefixed presets (@empty, @none)
- Auth required protection
"""

import uuid

import pytest
import pytest_asyncio
from conftest import create_test_token, create_test_user_model
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession
from zndraw_auth import User

from zndraw.exceptions import ProblemDetail
from zndraw.schemas import StatusResponse

# =============================================================================
# Fixtures
# =============================================================================


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
    """Create a regular user in the database."""
    user = create_test_user_model(email="regular@local.test")
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user


def auth_header(user: User) -> dict[str, str]:
    """Create Authorization header for a user."""
    token = create_test_token(user)
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# User Management Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_users(
    client: AsyncClient, admin_user: User, regular_user: User
) -> None:
    """Admin can list users."""
    response = await client.get("/v1/admin/users", headers=auth_header(admin_user))
    assert response.status_code == 200

    data = response.json()
    assert "items" in data
    assert "total" in data
    assert data["total"] >= 2  # At least admin_user and regular_user


@pytest.mark.asyncio
async def test_get_user(
    client: AsyncClient, admin_user: User, regular_user: User
) -> None:
    """Admin can get user details."""
    response = await client.get(
        f"/v1/admin/users/{regular_user.id}", headers=auth_header(admin_user)
    )
    assert response.status_code == 200

    data = response.json()
    assert data["id"] == str(regular_user.id)
    assert data["email"] == regular_user.email
    assert data["is_superuser"] is False


@pytest.mark.asyncio
async def test_get_user_not_found(client: AsyncClient, admin_user: User) -> None:
    """Getting a non-existent user returns 404."""
    fake_id = uuid.uuid4()
    response = await client.get(
        f"/v1/admin/users/{fake_id}", headers=auth_header(admin_user)
    )
    assert response.status_code == 404
    assert response.headers["content-type"] == "application/problem+json"

    problem = ProblemDetail.model_validate(response.json())
    assert problem.status == 404


@pytest.mark.asyncio
async def test_update_user_is_superuser(
    client: AsyncClient, admin_user: User, regular_user: User
) -> None:
    """Admin can promote a user to superuser."""
    response = await client.patch(
        f"/v1/admin/users/{regular_user.id}",
        headers=auth_header(admin_user),
        json={"is_superuser": True},
    )
    assert response.status_code == 200

    data = response.json()
    assert data["is_superuser"] is True


@pytest.mark.asyncio
async def test_cannot_demote_self(client: AsyncClient, admin_user: User) -> None:
    """Admin cannot demote themselves."""
    response = await client.patch(
        f"/v1/admin/users/{admin_user.id}",
        headers=auth_header(admin_user),
        json={"is_superuser": False},
    )
    assert response.status_code == 403
    assert response.headers["content-type"] == "application/problem+json"


@pytest.mark.asyncio
async def test_delete_user(
    client: AsyncClient, admin_user: User, regular_user: User
) -> None:
    """Admin can delete a user."""
    response = await client.delete(
        f"/v1/admin/users/{regular_user.id}", headers=auth_header(admin_user)
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify user is deleted
    response = await client.get(
        f"/v1/admin/users/{regular_user.id}", headers=auth_header(admin_user)
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_cannot_delete_self(client: AsyncClient, admin_user: User) -> None:
    """Admin cannot delete themselves."""
    response = await client.delete(
        f"/v1/admin/users/{admin_user.id}", headers=auth_header(admin_user)
    )
    assert response.status_code == 403
    assert response.headers["content-type"] == "application/problem+json"


# =============================================================================
# Server Settings Tests (Integration - requires Redis via server fixture)
# =============================================================================


async def _get_token(http_client: AsyncClient) -> str:
    """Get a guest token for testing."""
    response = await http_client.post("/v1/auth/guest")
    assert response.status_code == 200
    return response.json()["access_token"]


async def _create_room(http_client: AsyncClient, token: str) -> str:
    """Create a room and return its ID."""
    room_id = str(uuid.uuid4())
    response = await http_client.post(
        "/v1/rooms",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id, "description": "Test room"},
    )
    assert response.status_code == 201
    return room_id


@pytest.mark.asyncio
async def test_get_default_room_empty(server: str, http_client: AsyncClient) -> None:
    """Getting default room when none set returns None."""
    token = await _get_token(http_client)
    response = await http_client.get(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 200
    assert response.json()["room_id"] is None


@pytest.mark.asyncio
async def test_set_default_room_dev_mode(server: str, http_client: AsyncClient) -> None:
    """In dev mode, any user can set default room."""
    token = await _get_token(http_client)
    room_id = await _create_room(http_client, token)

    response = await http_client.put(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id},
    )
    assert response.status_code == 200
    assert response.json()["room_id"] == room_id

    # Verify it was set
    response = await http_client.get(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 200
    assert response.json()["room_id"] == room_id


@pytest.mark.asyncio
async def test_set_default_room_nonexistent(
    server: str, http_client: AsyncClient
) -> None:
    """Setting default room to non-existent room returns 404."""
    token = await _get_token(http_client)
    response = await http_client.put(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": "nonexistent-room"},
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_unset_default_room_dev_mode(
    server: str, http_client: AsyncClient
) -> None:
    """In dev mode, any user can unset default room."""
    token = await _get_token(http_client)
    room_id = await _create_room(http_client, token)

    # First set it
    await http_client.put(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id},
    )

    # Then unset it
    response = await http_client.delete(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify it was unset
    response = await http_client.get(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 200
    assert response.json()["room_id"] is None


@pytest.mark.asyncio
async def test_create_room_with_at_empty(server: str, http_client: AsyncClient) -> None:
    """copyFrom=@empty creates a room with one empty frame, bypassing default."""
    token = await _get_token(http_client)

    # Set a default room with frames
    default_id = await _create_room(http_client, token)

    await http_client.put(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": default_id},
    )

    # Create room with @empty — should NOT copy from default
    room_id = str(uuid.uuid4())
    response = await http_client.post(
        "/v1/rooms",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id, "copy_from": "@empty"},
    )
    assert response.status_code == 201
    assert response.json()["frame_count"] == 1


@pytest.mark.asyncio
async def test_create_room_with_at_none(server: str, http_client: AsyncClient) -> None:
    """copyFrom=@none creates a room with zero frames."""
    token = await _get_token(http_client)
    room_id = str(uuid.uuid4())
    response = await http_client.post(
        "/v1/rooms",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id, "copy_from": "@none"},
    )
    assert response.status_code == 201
    assert response.json()["frame_count"] == 0


@pytest.mark.asyncio
async def test_create_room_uses_default_when_no_copy_from(
    server: str, http_client: AsyncClient
) -> None:
    """Room creation without copyFrom uses server default."""
    token = await _get_token(http_client)

    # Create and populate a default room
    default_id = await _create_room(http_client, token)

    await http_client.put(
        "/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": default_id},
    )

    # Create room without copyFrom — should copy from default
    room_id = str(uuid.uuid4())
    response = await http_client.post(
        "/v1/rooms",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id},
    )
    assert response.status_code == 201
    # Default room has 1 frame, so new room should also have 1 frame
    assert response.json()["frame_count"] == 1


@pytest.mark.asyncio
async def test_create_room_with_invalid_preset(
    server: str, http_client: AsyncClient
) -> None:
    """copyFrom with unknown @-prefixed preset returns 422."""
    token = await _get_token(http_client)
    room_id = str(uuid.uuid4())
    response = await http_client.post(
        "/v1/rooms",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id, "copy_from": "@invalid"},
    )
    assert response.status_code == 422


# =============================================================================
# Auth Required Tests
# =============================================================================


@pytest.mark.asyncio
async def test_admin_endpoints_require_auth(client: AsyncClient) -> None:
    """Admin endpoints require authentication."""
    response = await client.get("/v1/admin/users")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_server_settings_mutation_requires_auth(
    server: str, http_client: AsyncClient
) -> None:
    """Server settings mutation endpoints require authentication."""
    response = await http_client.put(
        "/v1/server-settings/default-room",
        json={"room_id": "test-room"},
    )
    assert response.status_code == 401


# =============================================================================
# Pagination Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_users_pagination(
    client: AsyncClient, session: AsyncSession, admin_user: User
) -> None:
    """User list supports pagination."""
    # Create additional users
    for i in range(5):
        user = create_test_user_model(email=f"user{i}@local.test")
        session.add(user)
    await session.commit()

    # Test with limit
    response = await client.get(
        "/v1/admin/users?limit=3", headers=auth_header(admin_user)
    )
    assert response.status_code == 200

    data = response.json()
    assert len(data["items"]) == 3
    assert data["total"] >= 6  # admin + 5 created
    assert data["limit"] == 3
    assert data["offset"] == 0

    # Test with offset
    response = await client.get(
        "/v1/admin/users?offset=2&limit=2", headers=auth_header(admin_user)
    )
    assert response.status_code == 200

    data = response.json()
    assert len(data["items"]) == 2
    assert data["offset"] == 2
