"""Tests for Edit Lock REST API endpoints and WritableRoomDep enforcement."""

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
from zndraw.redis import RedisKey
from zndraw.schemas import StatusResponse
from zndraw.socket_events import LockUpdate

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="el_session")
async def el_session_fixture() -> AsyncIterator[AsyncSession]:
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


@pytest_asyncio.fixture(name="el_redis")
async def el_redis_fixture() -> AsyncIterator[Redis]:
    """Create a Redis client for edit lock testing."""
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


@pytest_asyncio.fixture(name="el_client")
async def el_client_fixture(
    el_session: AsyncSession,
    el_redis: Redis,
    mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with session, Redis, and sio overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield el_session

    def get_redis_override() -> Redis:
        return el_redis

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
# GET /v1/rooms/{room_id}/edit-lock
# =============================================================================


@pytest.mark.asyncio
async def test_get_edit_lock_returns_unlocked_when_no_lock(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test GET returns locked=False when no edit lock exists."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    response = await el_client.get(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["locked"] is False
    assert data["user_id"] is None


@pytest.mark.asyncio
async def test_get_edit_lock_returns_locked_when_lock_exists(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test GET returns lock info when an edit lock is held."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Set lock directly in Redis
    lock_data = json.dumps(
        {"user_id": str(user.id), "msg": "testing", "acquired_at": 1000.0}
    )
    await el_redis.set(RedisKey.edit_lock(room.id), lock_data, ex=10)

    response = await el_client.get(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["locked"] is True
    assert data["user_id"] == str(user.id)
    assert data["msg"] == "testing"
    assert data["ttl"] is not None
    assert data["ttl"] > 0


@pytest.mark.asyncio
async def test_get_edit_lock_returns_404_for_nonexistent_room(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await _create_user(el_session)

    response = await el_client.get(
        "/v1/rooms/nonexistent/edit-lock", headers=_auth(token)
    )
    assert response.status_code == 404


# =============================================================================
# PUT /v1/rooms/{room_id}/edit-lock (Acquire / Refresh)
# =============================================================================


@pytest.mark.asyncio
async def test_acquire_edit_lock(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test PUT acquires the edit lock."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    response = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing geometries"},
        headers=_auth(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["locked"] is True
    assert data["user_id"] == str(user.id)
    assert data["msg"] == "editing geometries"
    assert data["acquired_at"] is not None
    assert data["ttl"] is not None

    # Verify Redis entry
    raw = await el_redis.get(RedisKey.edit_lock(room.id))
    assert raw is not None
    lock = json.loads(raw)
    assert lock["user_id"] == str(user.id)


@pytest.mark.asyncio
async def test_acquire_edit_lock_broadcasts_lock_update(
    el_client: AsyncClient,
    el_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT broadcasts LockUpdate socket event."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "drawing"},
        headers=_auth(token),
    )

    mock_sio.emit.assert_called()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, LockUpdate)
    assert model.action == "acquired"
    assert model.user_id == str(user.id)
    assert model.msg == "drawing"


@pytest.mark.asyncio
async def test_refresh_edit_lock_idempotent(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test PUT refreshes when called by the same holder."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    resp1 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "first"},
        headers=_auth(token),
    )
    assert resp1.status_code == 200

    # Refresh (same user, different message)
    resp2 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "updated"},
        headers=_auth(token),
    )
    assert resp2.status_code == 200
    data = resp2.json()
    assert data["locked"] is True
    assert data["user_id"] == str(user.id)


@pytest.mark.asyncio
async def test_acquire_edit_lock_conflict_with_other_user(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test PUT returns 423 when another user holds the lock."""
    user1, token1 = await _create_user(el_session, "user1@test")
    user2, token2 = await _create_user(el_session, "user2@test")
    room = await _create_room(el_session, user1)

    # Give user2 room membership
    membership = RoomMembership(
        room_id=room.id,  # type: ignore
        user_id=user2.id,  # type: ignore
        role=MemberRole.MEMBER,
    )
    el_session.add(membership)
    await el_session.commit()

    # User1 acquires lock
    resp1 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "user1 editing"},
        headers=_auth(token1),
    )
    assert resp1.status_code == 200

    # User2 tries to acquire → 423
    resp2 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "user2 editing"},
        headers=_auth(token2),
    )
    assert resp2.status_code == 423
    assert "locked" in resp2.json()["type"]


@pytest.mark.asyncio
async def test_acquire_edit_lock_blocked_by_admin_lock(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test PUT returns 423 when room is admin-locked."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)
    room.locked = True
    el_session.add(room)
    await el_session.commit()

    response = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token),
    )
    assert response.status_code == 423


# =============================================================================
# DELETE /v1/rooms/{room_id}/edit-lock (Release)
# =============================================================================


@pytest.mark.asyncio
async def test_release_edit_lock(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test DELETE releases the lock."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "temp"},
        headers=_auth(token),
    )

    # Release
    response = await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token)
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify Redis entry removed
    raw = await el_redis.get(RedisKey.edit_lock(room.id))
    assert raw is None


@pytest.mark.asyncio
async def test_release_edit_lock_broadcasts_lock_update(
    el_client: AsyncClient,
    el_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE broadcasts LockUpdate with action=released."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire then release
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=_auth(token),
    )
    mock_sio.emit.reset_mock()

    await el_client.delete(f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token))

    mock_sio.emit.assert_called()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, LockUpdate)
    assert model.action == "released"


@pytest.mark.asyncio
async def test_release_edit_lock_idempotent_when_no_lock(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test DELETE succeeds even if no lock exists."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    response = await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token)
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())


@pytest.mark.asyncio
async def test_release_edit_lock_forbidden_for_non_holder(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test DELETE returns 403 when non-holder tries to release."""
    user1, token1 = await _create_user(el_session, "user1@test")
    user2, token2 = await _create_user(el_session, "user2@test")
    room = await _create_room(el_session, user1)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore
        user_id=user2.id,  # type: ignore
        role=MemberRole.MEMBER,
    )
    el_session.add(membership)
    await el_session.commit()

    # User1 acquires
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=_auth(token1),
    )

    # User2 tries to release → 403
    response = await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token2)
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_admin_can_release_any_lock(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test admin can release another user's lock."""
    user, token = await _create_user(el_session, "user@test")
    admin, admin_token = await _create_user(el_session, "admin@test", is_superuser=True)
    room = await _create_room(el_session, user)

    # Give admin membership
    membership = RoomMembership(
        room_id=room.id,  # type: ignore
        user_id=admin.id,  # type: ignore
        role=MemberRole.MEMBER,
    )
    el_session.add(membership)
    await el_session.commit()

    # User acquires
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=_auth(token),
    )

    # Admin releases
    response = await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(admin_token)
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())


# =============================================================================
# WritableRoomDep: Admin Lock Enforcement
# =============================================================================


@pytest.mark.asyncio
async def test_writable_room_blocks_non_admin_on_admin_locked_room(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test mutation returns 423 when room is admin-locked and user is not admin."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)
    room.locked = True
    el_session.add(room)
    await el_session.commit()

    # Try to set a bookmark (mutation endpoint using WritableRoomDep)
    response = await el_client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "test"},
        headers=_auth(token),
    )
    assert response.status_code == 423
    assert "locked" in response.json()["type"]


@pytest.mark.asyncio
async def test_writable_room_allows_admin_on_admin_locked_room(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test admin can mutate even when room is admin-locked."""
    admin, token = await _create_user(el_session, is_superuser=True)
    room = await _create_room(el_session, admin)
    room.locked = True
    el_session.add(room)
    await el_session.commit()

    response = await el_client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "admin edit"},
        headers=_auth(token),
    )
    assert response.status_code == 200


# =============================================================================
# WritableRoomDep: Edit Lock Enforcement
# =============================================================================


@pytest.mark.asyncio
async def test_writable_room_blocks_non_holder(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test mutation returns 423 when another user holds the edit lock."""
    user1, token1 = await _create_user(el_session, "user1@test")
    user2, token2 = await _create_user(el_session, "user2@test")
    room = await _create_room(el_session, user1)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore
        user_id=user2.id,  # type: ignore
        role=MemberRole.MEMBER,
    )
    el_session.add(membership)
    await el_session.commit()

    # User1 acquires lock
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token1),
    )

    # User2 tries mutation → 423
    response = await el_client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "blocked"},
        headers=_auth(token2),
    )
    assert response.status_code == 423
    assert "locked" in response.json()["type"]


@pytest.mark.asyncio
async def test_writable_room_allows_lock_holder(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test lock holder can mutate normally."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire lock
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token),
    )

    # Holder can still mutate
    response = await el_client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "allowed"},
        headers=_auth(token),
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_writable_room_allows_get_when_locked(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test GET endpoints still work when room has edit lock."""
    user1, token1 = await _create_user(el_session, "user1@test")
    user2, token2 = await _create_user(el_session, "user2@test")
    room = await _create_room(el_session, user1)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore
        user_id=user2.id,  # type: ignore
        role=MemberRole.MEMBER,
    )
    el_session.add(membership)
    await el_session.commit()

    # User1 acquires lock
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token1),
    )

    # User2 can still read bookmarks
    response = await el_client.get(
        f"/v1/rooms/{room.id}/bookmarks",
        headers=_auth(token2),
    )
    assert response.status_code == 200
