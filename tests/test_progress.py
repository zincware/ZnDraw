"""Tests for Progress REST API endpoints."""

import json
from collections.abc import AsyncIterator
from unittest.mock import AsyncMock

import pytest
import pytest_asyncio
from conftest import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from zndraw_auth.settings import AuthSettings

from zndraw.config import Settings
from zndraw.redis import RedisKey

# =============================================================================
# Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="progress_session")
async def progress_session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session for each test."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)
        factory = async_sessionmaker(
            bind=engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MockSioServer:
    return MockSioServer()


@pytest_asyncio.fixture(name="mock_redis")
async def mock_redis_fixture() -> AsyncMock:
    redis = AsyncMock()
    redis.get = AsyncMock(return_value=None)  # no edit lock
    redis._progress_store: dict[str, dict[str, str]] = {}

    async def hset(key: str, field: str, value: str) -> int:
        if key not in redis._progress_store:
            redis._progress_store[key] = {}
        redis._progress_store[key][field] = value
        return 1

    async def hget(key: str, field: str) -> str | None:
        return redis._progress_store.get(key, {}).get(field)

    async def hdel(key: str, *fields: str) -> int:
        deleted = 0
        if key in redis._progress_store:
            for f in fields:
                if f in redis._progress_store[key]:
                    del redis._progress_store[key][f]
                    deleted += 1
        return deleted

    async def hgetall(key: str) -> dict[str, str]:
        return redis._progress_store.get(key, {})

    redis.hset = AsyncMock(side_effect=hset)
    redis.hget = AsyncMock(side_effect=hget)
    redis.hdel = AsyncMock(side_effect=hdel)
    redis.hgetall = AsyncMock(side_effect=hgetall)
    redis.expire = AsyncMock(return_value=True)
    return redis


@pytest_asyncio.fixture(name="progress_client")
async def progress_client_fixture(
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> AsyncIterator[AsyncClient]:
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield progress_session

    def get_sio_override() -> MockSioServer:
        return mock_sio

    def get_redis_override() -> AsyncMock:
        return mock_redis

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_tsio] = get_sio_override
    app.dependency_overrides[get_redis] = get_redis_override
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        yield client

    app.dependency_overrides.clear()


# =============================================================================
# POST /v1/rooms/{room_id}/progress
# =============================================================================


@pytest.mark.asyncio
async def test_create_progress(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """POST creates tracker, returns 201, emits progress_start, stores in Redis."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.post(
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
    stored = await mock_redis.hget(redis_key, "task-1")
    assert stored is not None
    stored_data = json.loads(stored)
    assert stored_data["progress_id"] == "task-1"
    assert stored_data["description"] == "Loading data"

    # Verify TTL was set
    mock_redis.expire.assert_called()


@pytest.mark.asyncio
async def test_create_progress_with_unit(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """POST with custom unit returns it in the response."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.post(
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
    progress_client: AsyncClient,
    progress_session: AsyncSession,
) -> None:
    """POST without token returns 401."""
    user, _ = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
    )
    assert response.status_code == 401


# =============================================================================
# PATCH /v1/rooms/{room_id}/progress/{progress_id}
# =============================================================================


@pytest.mark.asyncio
async def test_update_progress(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """PATCH updates tqdm fields, returns 200, emits progress_update."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    # Create a tracker first
    await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    mock_sio.emitted.clear()

    # Update with tqdm-like fields
    response = await progress_client.patch(
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
    progress_client: AsyncClient,
    progress_session: AsyncSession,
) -> None:
    """PATCH for non-existent tracker returns 404 with progress-not-found type."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.patch(
        f"/v1/rooms/{room.id}/progress/nonexistent",
        json={"n": 10},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "progress-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_update_progress_description(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """PATCH can update description alongside tqdm fields."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    # Create a tracker first
    await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    mock_sio.emitted.clear()

    # Update both description and tqdm fields
    response = await progress_client.patch(
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
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """DELETE removes tracker, returns 204, emits progress_complete, removed from Redis."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    # Create a tracker first
    await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    mock_sio.emitted.clear()

    # Delete it
    response = await progress_client.delete(
        f"/v1/rooms/{room.id}/progress/task-1",
        headers=auth_header(token),
    )
    assert response.status_code == 204

    # Verify socket broadcast
    complete_events = [e for e in mock_sio.emitted if e["event"] == "progress_complete"]
    assert len(complete_events) == 1

    # Verify removed from Redis
    redis_key = RedisKey.room_progress(room.id)
    stored = await mock_redis.hget(redis_key, "task-1")
    assert stored is None


@pytest.mark.asyncio
async def test_delete_progress_not_found(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
) -> None:
    """DELETE for non-existent tracker returns 404."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.delete(
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
    progress_client: AsyncClient,
    progress_session: AsyncSession,
) -> None:
    """POST to non-existent room returns 404 with room-not-found type."""
    _, token = await create_test_user_in_db(progress_session)

    response = await progress_client.post(
        "/v1/rooms/nonexistent/progress",
        json={"progress_id": "task-1", "description": "Loading data"},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
