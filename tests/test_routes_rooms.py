"""Tests for room REST API endpoints."""

import pytest
from httpx import AsyncClient

from zndraw.schemas import RoomCreateResponse


async def _get_user_token(
    http_client: AsyncClient, email: str, password: str = "testpassword"
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


@pytest.mark.asyncio
async def test_create_room_with_valid_uuid(http_client: AsyncClient):
    """Test creating a room with a valid UUID format."""
    token = await _get_user_token(http_client, "validuuid@example.com")
    room_id = "abc-123-def-456"
    response = await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 201
    result = RoomCreateResponse.model_validate(response.json())
    assert result.room_id == room_id
    assert result.status == "ok"
    assert result.created is True


@pytest.mark.asyncio
async def test_create_room_with_alphanumeric_only(http_client: AsyncClient):
    """Test creating a room with alphanumeric characters only."""
    token = await _get_user_token(http_client, "alphanumeric@example.com")
    room_id = "abc123def456"
    response = await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 201
    result = RoomCreateResponse.model_validate(response.json())
    assert result.room_id == room_id


@pytest.mark.asyncio
async def test_create_room_with_underscores(http_client: AsyncClient):
    """Test creating a room with underscores (now allowed)."""
    token = await _get_user_token(http_client, "underscores@example.com")
    room_id = "test_file_123"
    response = await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 201
    result = RoomCreateResponse.model_validate(response.json())
    assert result.room_id == room_id


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "invalid_room_id",
    [
        "@overview",  # System room prefix not allowed
        "room with spaces",  # Spaces not allowed
        "room@example",  # @ in middle not allowed
        "room.with.dots",  # Dots not allowed
        "room/with/slashes",  # Slashes not allowed
        "room\\with\\backslashes",  # Backslashes not allowed
        "room:with:colons",  # Colons not allowed
    ],
)
async def test_create_room_with_invalid_characters(
    http_client: AsyncClient, invalid_room_id: str
):
    """Test that room creation rejects invalid characters."""
    token = await _get_user_token(
        http_client, f"invalid{hash(invalid_room_id)}@example.com"
    )
    response = await http_client.post(
        "/v1/rooms",
        json={"room_id": invalid_room_id},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 400
    # RFC 9457 Problem JSON response
    problem = response.json()
    assert "type" in problem
    assert "detail" in problem
    assert "alphanumeric" in problem["detail"].lower()
