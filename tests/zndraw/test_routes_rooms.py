"""Tests for room REST API endpoints."""

import pytest
from httpx import AsyncClient

from zndraw.schemas import RoomCreateResponse


async def _get_user_token(
    http_client: AsyncClient,
    email: str,
    password: str = "testpassword",
) -> str:
    """Register a user and get their auth token."""
    reg_response = await http_client.post(
        "/v1/auth/register", json={"email": email, "password": password}
    )
    assert reg_response.status_code == 201, f"Register failed: {reg_response.text}"
    login_response = await http_client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": password},
    )
    assert login_response.status_code == 200, f"Login failed: {login_response.text}"
    return login_response.json()["access_token"]


async def _register_and_login(
    http_client: AsyncClient,
    email: str,
    password: str = "testpassword",
) -> str:
    """Alias of ``_get_user_token`` — kept for readability in new tests."""
    return await _get_user_token(http_client, email, password)


async def _get_user_display_name(http_client: AsyncClient, token: str) -> str:
    """Get the authenticated user's display_name."""
    me_resp = await http_client.get(
        "/v1/auth/users/me",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert me_resp.status_code == 200
    return me_resp.json()["display_name"]


@pytest.mark.asyncio
async def test_create_room_with_valid_name(http_client: AsyncClient):
    """Test creating a room with a valid name and display-name owner."""
    token = await _get_user_token(http_client, "validname@example.com")
    display_name = await _get_user_display_name(http_client, token)
    room_name = "abc-123-def-456"
    response = await http_client.post(
        "/v1/rooms",
        json={"owner": display_name, "name": room_name},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 201
    result = RoomCreateResponse.model_validate(response.json())
    assert result.room_id == f"{display_name}/{room_name}"
    assert result.status == "ok"
    assert result.created is True


@pytest.mark.asyncio
async def test_create_room_with_alphanumeric_only(http_client: AsyncClient):
    """Test creating a room with alphanumeric characters only."""
    token = await _get_user_token(http_client, "alphanumeric@example.com")
    display_name = await _get_user_display_name(http_client, token)
    room_name = "abc123def456"
    response = await http_client.post(
        "/v1/rooms",
        json={"owner": display_name, "name": room_name},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 201
    result = RoomCreateResponse.model_validate(response.json())
    assert result.room_id == f"{display_name}/{room_name}"


@pytest.mark.asyncio
async def test_create_room_with_underscores(http_client: AsyncClient):
    """Test creating a room with underscores (allowed)."""
    token = await _get_user_token(http_client, "underscores@example.com")
    display_name = await _get_user_display_name(http_client, token)
    room_name = "test_file_123"
    response = await http_client.post(
        "/v1/rooms",
        json={"owner": display_name, "name": room_name},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 201
    result = RoomCreateResponse.model_validate(response.json())
    assert result.room_id == f"{display_name}/{room_name}"


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "invalid_room_name",
    [
        "room with spaces",  # Spaces not allowed
        "room@example",  # @ not allowed
        "room.with.dots",  # Dots not allowed
        "room/with/slashes",  # Slashes not allowed
        "room\\with\\backslashes",  # Backslashes not allowed
        "room:with:colons",  # Colons not allowed
    ],
)
async def test_create_room_with_invalid_characters(
    http_client: AsyncClient, invalid_room_name: str
):
    """Test that room creation rejects invalid name characters."""
    token = await _get_user_token(
        http_client, f"invalid{abs(hash(invalid_room_name)) % 10**9}@example.com"
    )
    display_name = await _get_user_display_name(http_client, token)
    response = await http_client.post(
        "/v1/rooms",
        json={"owner": display_name, "name": invalid_room_name},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_owner_segment_uuid_returns_422(
    http_client: AsyncClient,
) -> None:
    """A UUID-shaped owner segment must hit the path regex gate (422)."""
    email = "rooms-uuid-422@example.com"
    token = await _register_and_login(http_client, email)
    fake_uuid = "00000000-0000-0000-0000-000000000000"
    resp = await http_client.get(
        f"/v1/rooms/{fake_uuid}/any-room",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 422


@pytest.mark.asyncio
async def test_owner_segment_unknown_display_name_returns_404(
    http_client: AsyncClient,
) -> None:
    """Path-regex-valid but unknown owner display_name resolves to 404."""
    email = "rooms-404@example.com"
    token = await _register_and_login(http_client, email)
    resp = await http_client.get(
        "/v1/rooms/no-such-display-name/any-room",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 404
