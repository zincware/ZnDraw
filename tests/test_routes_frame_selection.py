"""Tests for frame selection REST API endpoints."""

import json
from collections.abc import AsyncIterator
from unittest.mock import AsyncMock, MagicMock

import pytest
import pytest_asyncio
from conftest import create_test_token, create_test_user_model
from httpx import ASGITransport, AsyncClient
from redis.asyncio import Redis
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from zndraw_auth import User
from zndraw_auth.settings import AuthSettings

from zndraw.config import Settings
from zndraw.models import MemberRole, Room, RoomMembership
from zndraw.socket_events import FrameSelectionUpdate

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="fs_session")
async def fs_session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session for each test."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )

    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)

        async_session_factory = async_sessionmaker(
            bind=engine,
            class_=AsyncSession,
            expire_on_commit=False,
        )

        async with async_session_factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="fs_redis")
async def fs_redis_fixture() -> AsyncIterator[Redis]:
    """Create a Redis client for frame selection testing."""
    redis: Redis = Redis.from_url("redis://localhost", decode_responses=True)
    await redis.flushdb()
    yield redis
    await redis.flushdb()
    await redis.aclose()


@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MagicMock:
    """Create a mock Socket.IO server for testing."""
    sio_mock = MagicMock()
    sio_mock.emit = AsyncMock()
    return sio_mock


@pytest_asyncio.fixture(name="fs_client")
async def fs_client_fixture(
    fs_session: AsyncSession,
    fs_redis: Redis,
    mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with session, Redis, and sio overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield fs_session

    def get_redis_override() -> Redis:
        return fs_redis

    def get_sio_override() -> MagicMock:
        return mock_sio

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_redis] = get_redis_override
    app.dependency_overrides[get_tsio] = get_sio_override
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client

    app.dependency_overrides.clear()


# =============================================================================
# Helpers
# =============================================================================


async def _create_user(
    session: AsyncSession,
    email: str = "testuser@local.test",
    is_superuser: bool = False,
) -> tuple[User, str]:
    """Create a user and return the user and access token."""
    user = create_test_user_model(email=email, is_superuser=is_superuser)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    token = create_test_token(user)
    return user, token


async def _create_room(
    session: AsyncSession, user: User, description: str = "Test Room"
) -> Room:
    """Create a room with user as owner."""
    room = Room(
        description=description,
        created_by_id=user.id,  # type: ignore
        is_public=True,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore
        user_id=user.id,  # type: ignore
        role=MemberRole.OWNER,
    )
    session.add(membership)
    await session.commit()
    return room


def _auth(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# GET /v1/rooms/{room_id}/frame-selection
# =============================================================================


@pytest.mark.asyncio
async def test_get_returns_null_when_empty(
    fs_client: AsyncClient, fs_session: AsyncSession
) -> None:
    """GET returns null frameSelection for a new room."""
    user, token = await _create_user(fs_session)
    room = await _create_room(fs_session, user)

    response = await fs_client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=_auth(token)
    )
    assert response.status_code == 200
    assert response.json()["frame_selection"] is None


@pytest.mark.asyncio
async def test_get_returns_stored_indices(
    fs_client: AsyncClient, fs_session: AsyncSession
) -> None:
    """GET returns stored indices when frame_selection column is set."""
    user, token = await _create_user(fs_session)
    room = await _create_room(fs_session, user)

    # Set column directly
    room.frame_selection = json.dumps([2, 5, 10])
    fs_session.add(room)
    await fs_session.commit()

    response = await fs_client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=_auth(token)
    )
    assert response.status_code == 200
    assert response.json()["frame_selection"] == [2, 5, 10]


@pytest.mark.asyncio
async def test_get_returns_404_for_nonexistent_room(
    fs_client: AsyncClient, fs_session: AsyncSession
) -> None:
    """GET for non-existent room returns 404."""
    _, token = await _create_user(fs_session)

    response = await fs_client.get(
        "/v1/rooms/nonexistent/frame-selection", headers=_auth(token)
    )
    assert response.status_code == 404


# =============================================================================
# PUT /v1/rooms/{room_id}/frame-selection
# =============================================================================


@pytest.mark.asyncio
async def test_put_stores_indices(
    fs_client: AsyncClient, fs_session: AsyncSession
) -> None:
    """PUT stores indices and GET returns them."""
    user, token = await _create_user(fs_session)
    room = await _create_room(fs_session, user)

    put_resp = await fs_client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [1, 3, 7]},
        headers=_auth(token),
    )
    assert put_resp.status_code == 200
    assert put_resp.json()["success"] is True

    get_resp = await fs_client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=_auth(token)
    )
    assert get_resp.json()["frame_selection"] == [1, 3, 7]


@pytest.mark.asyncio
async def test_put_broadcasts_socket_event(
    fs_client: AsyncClient,
    fs_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """PUT broadcasts FrameSelectionUpdate socket event."""
    user, token = await _create_user(fs_session)
    room = await _create_room(fs_session, user)

    await fs_client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [0, 4]},
        headers=_auth(token),
    )

    mock_sio.emit.assert_called()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, FrameSelectionUpdate)
    assert model.indices == [0, 4]


@pytest.mark.asyncio
async def test_put_empty_list_clears_selection(
    fs_client: AsyncClient, fs_session: AsyncSession
) -> None:
    """PUT with empty list clears selection (GET returns null)."""
    user, token = await _create_user(fs_session)
    room = await _create_room(fs_session, user)

    # Set some indices first
    await fs_client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [1, 2]},
        headers=_auth(token),
    )

    # Clear with empty list
    await fs_client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": []},
        headers=_auth(token),
    )

    get_resp = await fs_client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=_auth(token)
    )
    assert get_resp.json()["frame_selection"] is None


@pytest.mark.asyncio
async def test_put_rejects_negative_indices(
    fs_client: AsyncClient, fs_session: AsyncSession
) -> None:
    """PUT with negative indices returns 400."""
    user, token = await _create_user(fs_session)
    room = await _create_room(fs_session, user)

    response = await fs_client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [1, -2, 3]},
        headers=_auth(token),
    )
    assert response.status_code == 400


@pytest.mark.asyncio
async def test_put_returns_404_for_nonexistent_room(
    fs_client: AsyncClient, fs_session: AsyncSession
) -> None:
    """PUT for non-existent room returns 404."""
    _, token = await _create_user(fs_session)

    response = await fs_client.put(
        "/v1/rooms/nonexistent/frame-selection",
        json={"indices": [0]},
        headers=_auth(token),
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_roundtrip(fs_client: AsyncClient, fs_session: AsyncSession) -> None:
    """PUT then GET returns consistent data."""
    user, token = await _create_user(fs_session)
    room = await _create_room(fs_session, user)
    indices = [0, 2, 4, 6, 8]

    await fs_client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": indices},
        headers=_auth(token),
    )

    get_resp = await fs_client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=_auth(token)
    )
    assert get_resp.status_code == 200
    assert get_resp.json()["frame_selection"] == indices
