"""Tests for Edit Lock REST API endpoints and WritableRoomDep enforcement."""

import asyncio
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
    assert data["lock_token"] is None


@pytest.mark.asyncio
async def test_get_edit_lock_returns_locked_when_lock_exists(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test GET returns lock info when an edit lock is held."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Set lock directly in Redis (new format with lock_token)
    lock_data = json.dumps(
        {
            "lock_token": "test-token-123",
            "user_id": str(user.id),
            "sid": None,
            "msg": "testing",
            "acquired_at": 1000.0,
        }
    )
    await el_redis.set(RedisKey.edit_lock(room.id), lock_data, ex=10)

    response = await el_client.get(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["locked"] is True
    assert data["user_id"] == str(user.id)
    assert data["lock_token"] == "test-token-123"
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
    """Test PUT acquires the edit lock and returns a lock_token."""
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
    assert data["lock_token"] is not None
    assert data["msg"] == "editing geometries"
    assert data["acquired_at"] is not None
    assert data["ttl"] is not None

    # Verify Redis entry has lock_token
    raw = await el_redis.get(RedisKey.edit_lock(room.id))
    assert raw is not None
    lock = json.loads(raw)
    assert lock["user_id"] == str(user.id)
    assert lock["lock_token"] == data["lock_token"]


@pytest.mark.asyncio
async def test_acquire_edit_lock_stores_session_id(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test PUT stores X-Session-ID header as sid in Redis."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    response = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers={**_auth(token), "X-Session-ID": "my-sid-123"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["sid"] == "my-sid-123"

    # Verify Redis
    raw = await el_redis.get(RedisKey.edit_lock(room.id))
    lock = json.loads(raw)
    assert lock["sid"] == "my-sid-123"


@pytest.mark.asyncio
async def test_acquire_edit_lock_broadcasts_lock_update(
    el_client: AsyncClient,
    el_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT broadcasts LockUpdate socket event with ttl."""
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
    assert model.ttl is not None
    assert model.ttl > 0


@pytest.mark.asyncio
async def test_refresh_edit_lock_with_token(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test PUT with Lock-Token header refreshes the lock."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    resp1 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "first"},
        headers=_auth(token),
    )
    assert resp1.status_code == 200
    lock_token = resp1.json()["lock_token"]

    # Refresh with Lock-Token
    resp2 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "updated"},
        headers={**_auth(token), "Lock-Token": lock_token},
    )
    assert resp2.status_code == 200
    data = resp2.json()
    assert data["locked"] is True
    assert data["lock_token"] == lock_token
    assert data["user_id"] == str(user.id)


@pytest.mark.asyncio
async def test_refresh_with_expired_lock_returns_409(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test refresh returns 409 when lock has expired."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    resp = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token),
    )
    lock_token = resp.json()["lock_token"]

    # Simulate expiry by deleting the key
    await el_redis.delete(RedisKey.edit_lock(room.id))

    # Try to refresh → 409
    resp2 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers={**_auth(token), "Lock-Token": lock_token},
    )
    assert resp2.status_code == 409


@pytest.mark.asyncio
async def test_acquire_without_token_when_lock_exists_returns_423(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test PUT without Lock-Token when lock exists returns 423 (no silent re-create)."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    resp1 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "first"},
        headers=_auth(token),
    )
    assert resp1.status_code == 200

    # Same user, no Lock-Token → 423 (lock already exists)
    resp2 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "second attempt"},
        headers=_auth(token),
    )
    assert resp2.status_code == 423


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


@pytest.mark.asyncio
async def test_lock_auto_expires_after_ttl(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Lock disappears from Redis after edit_lock_ttl without refresh."""
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
    resp = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "no refresh"},
        headers=_auth(token1),
    )
    assert resp.status_code == 200

    # Wait for TTL + buffer (default 10s + 1s)
    await asyncio.sleep(11)

    # Lock should be gone
    status = await el_client.get(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token1)
    )
    assert status.json()["locked"] is False

    # User2 can now acquire
    resp2 = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "after expiry"},
        headers=_auth(token2),
    )
    assert resp2.status_code == 200
    assert resp2.json()["user_id"] == str(user2.id)


# =============================================================================
# DELETE /v1/rooms/{room_id}/edit-lock (Release)
# =============================================================================


@pytest.mark.asyncio
async def test_release_edit_lock_with_token(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test DELETE with Lock-Token releases the lock."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    resp = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "temp"},
        headers=_auth(token),
    )
    lock_token = resp.json()["lock_token"]

    # Release with Lock-Token
    response = await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock",
        headers={**_auth(token), "Lock-Token": lock_token},
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify Redis entry removed
    raw = await el_redis.get(RedisKey.edit_lock(room.id))
    assert raw is None


@pytest.mark.asyncio
async def test_release_edit_lock_by_user_id_fallback(
    el_client: AsyncClient, el_session: AsyncSession, el_redis: Redis
) -> None:
    """Test DELETE without Lock-Token falls back to user_id check."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "temp"},
        headers=_auth(token),
    )

    # Release without Lock-Token (user_id fallback)
    response = await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=_auth(token)
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

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
    resp = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=_auth(token),
    )
    lock_token = resp.json()["lock_token"]
    mock_sio.emit.reset_mock()

    await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock",
        headers={**_auth(token), "Lock-Token": lock_token},
    )

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
async def test_release_edit_lock_wrong_token_returns_403(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test DELETE with wrong Lock-Token returns 403."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=_auth(token),
    )

    # Release with wrong Lock-Token
    response = await el_client.delete(
        f"/v1/rooms/{room.id}/edit-lock",
        headers={**_auth(token), "Lock-Token": "wrong-token"},
    )
    assert response.status_code == 403


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

    # Admin releases (no Lock-Token needed for admin)
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
# WritableRoomDep: Edit Lock Enforcement with Lock-Token
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
async def test_writable_room_allows_lock_holder_with_token(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test lock holder can mutate when sending Lock-Token header."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire lock
    resp = await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token),
    )
    lock_token = resp.json()["lock_token"]

    # Holder can mutate with Lock-Token
    response = await el_client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "allowed"},
        headers={**_auth(token), "Lock-Token": lock_token},
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_writable_room_allows_lock_holder_by_user_id(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test lock holder can mutate by user_id fallback (no Lock-Token)."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire lock
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token),
    )

    # Holder can mutate without Lock-Token (user_id fallback)
    response = await el_client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "allowed"},
        headers=_auth(token),
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_writable_room_blocks_wrong_lock_token(
    el_client: AsyncClient, el_session: AsyncSession
) -> None:
    """Test mutation returns 423 when wrong Lock-Token is sent."""
    user, token = await _create_user(el_session)
    room = await _create_room(el_session, user)

    # Acquire lock
    await el_client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=_auth(token),
    )

    # Wrong Lock-Token → 423
    response = await el_client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "blocked"},
        headers={**_auth(token), "Lock-Token": "wrong-token"},
    )
    assert response.status_code == 423


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
