"""Tests for Progress REST API endpoints."""

import json

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from redis.asyncio import Redis
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.redis import RedisKey


# =============================================================================
# POST /v1/rooms/{room_id}/progress
# =============================================================================


@pytest.mark.asyncio
async def test_create_progress(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
    redis_client: Redis,
) -> None:
    """POST creates tracker, returns 201, emits progress_start, stores in Redis."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    data = response.json()
    assert data["progress_id"] == "task-1"
    assert data["description"] == "Loading data"
    assert data["n"] == 0
    assert data["total"] is None
    assert data["elapsed"] == 0.0
    assert data["unit"] == "it"

    # Verify socket broadcast
    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "progress_start"

    # Verify stored in Redis
    redis_key = RedisKey.room_progress(room.id)
    stored = await redis_client.hget(redis_key, "task-1")
    assert stored is not None
    stored_data = json.loads(stored)
    assert stored_data["progress_id"] == "task-1"
    assert stored_data["description"] == "Loading data"


@pytest.mark.asyncio
async def test_create_progress_with_unit(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """POST with custom unit returns it in the response."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/progress",
        json={
            "progress_id": "task-1",
            "description": "Uploading",
            "unit": "frames",
        },
        headers=auth_header(token),
    )
    assert response.status_code == 201
    assert response.json()["unit"] == "frames"

    # Verify broadcast includes unit
    assert mock_sio.emitted[0]["data"]["unit"] == "frames"


@pytest.mark.asyncio
async def test_create_progress_requires_auth(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """POST without token returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
    )
    assert response.status_code == 401


# =============================================================================
# PATCH /v1/rooms/{room_id}/progress/{progress_id}
# =============================================================================


@pytest.mark.asyncio
async def test_update_progress(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """PATCH updates tqdm fields, returns 200, emits progress_update."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Create a tracker first
    await client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    mock_sio.emitted.clear()

    # Update with tqdm-like fields
    response = await client.patch(
        f"/v1/rooms/{room.id}/progress/task-1",
        json={"n": 42, "total": 100, "elapsed": 5.3, "unit": "frames"},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["progress_id"] == "task-1"
    assert data["n"] == 42
    assert data["total"] == 100
    assert data["elapsed"] == 5.3
    assert data["unit"] == "frames"

    # Verify socket broadcast
    update_events = [e for e in mock_sio.emitted if e["event"] == "progress_update"]
    assert len(update_events) == 1
    assert update_events[0]["data"]["n"] == 42


@pytest.mark.asyncio
async def test_update_progress_not_found(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """PATCH for non-existent tracker returns 404 with progress-not-found type."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.patch(
        f"/v1/rooms/{room.id}/progress/nonexistent",
        json={"n": 10},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "progress-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_update_progress_description(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """PATCH can update description alongside tqdm fields."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Create a tracker first
    await client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    mock_sio.emitted.clear()

    # Update both description and tqdm fields
    response = await client.patch(
        f"/v1/rooms/{room.id}/progress/task-1",
        json={
            "description": "Processing step 2",
            "n": 75,
            "total": 100,
            "elapsed": 8.1,
        },
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["description"] == "Processing step 2"
    assert data["n"] == 75
    assert data["total"] == 100
    assert data["elapsed"] == 8.1

    # Verify socket broadcast
    update_events = [e for e in mock_sio.emitted if e["event"] == "progress_update"]
    assert len(update_events) == 1


# =============================================================================
# DELETE /v1/rooms/{room_id}/progress/{progress_id}
# =============================================================================


@pytest.mark.asyncio
async def test_delete_progress(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
    redis_client: Redis,
) -> None:
    """DELETE removes tracker, returns 204, emits progress_complete."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Create a tracker first
    await client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    mock_sio.emitted.clear()

    # Delete it
    response = await client.delete(
        f"/v1/rooms/{room.id}/progress/task-1",
        headers=auth_header(token),
    )
    assert response.status_code == 204

    # Verify socket broadcast
    complete_events = [e for e in mock_sio.emitted if e["event"] == "progress_complete"]
    assert len(complete_events) == 1

    # Verify removed from Redis
    redis_key = RedisKey.room_progress(room.id)
    stored = await redis_client.hget(redis_key, "task-1")
    assert stored is None


@pytest.mark.asyncio
async def test_delete_progress_not_found(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """DELETE for non-existent tracker returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(
        f"/v1/rooms/{room.id}/progress/nonexistent",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "progress-not-found" in response.json()["type"]


# =============================================================================
# Room not found
# =============================================================================


@pytest.mark.asyncio
async def test_progress_room_not_found(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """POST to non-existent room returns 404 with room-not-found type."""
    _, token = await create_test_user_in_db(session)

    response = await client.post(
        "/v1/rooms/nonexistent/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
