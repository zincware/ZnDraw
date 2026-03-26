"""Tests for Edit Lock REST API endpoints and WritableRoomDep enforcement."""

import asyncio
import json

import pytest
from helpers import auth_header, create_test_room, create_test_user_in_db
from httpx import AsyncClient
from redis.asyncio import Redis
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.models import MemberRole, RoomMembership
from zndraw.redis import RedisKey
from zndraw.schemas import StatusResponse

# =============================================================================
# GET /v1/rooms/{room_id}/edit-lock
# =============================================================================


@pytest.mark.asyncio
async def test_get_edit_lock_returns_unlocked_when_no_lock(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET returns locked=False when no edit lock exists."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/edit-lock", headers=auth_header(token)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["locked"] is False
    assert data["user_id"] is None
    assert data["lock_token"] is None


@pytest.mark.asyncio
async def test_get_edit_lock_returns_locked_when_lock_exists(
    client: AsyncClient, session: AsyncSession, redis_client: Redis
) -> None:
    """Test GET returns lock info when an edit lock is held."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

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
    await redis_client.set(RedisKey.edit_lock(room.id), lock_data, ex=10)

    response = await client.get(
        f"/v1/rooms/{room.id}/edit-lock", headers=auth_header(token)
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
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/nonexistent/edit-lock", headers=auth_header(token)
    )
    assert response.status_code == 404


# =============================================================================
# PUT /v1/rooms/{room_id}/edit-lock (Acquire / Refresh)
# =============================================================================


@pytest.mark.asyncio
async def test_acquire_edit_lock(
    client: AsyncClient, session: AsyncSession, redis_client: Redis
) -> None:
    """Test PUT acquires the edit lock and returns a lock_token."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing geometries"},
        headers=auth_header(token),
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
    raw = await redis_client.get(RedisKey.edit_lock(room.id))
    assert raw is not None
    lock = json.loads(raw)
    assert lock["user_id"] == str(user.id)
    assert lock["lock_token"] == data["lock_token"]


@pytest.mark.asyncio
async def test_acquire_edit_lock_stores_session_id(
    client: AsyncClient, session: AsyncSession, redis_client: Redis
) -> None:
    """Test PUT stores X-Session-ID header as sid in Redis."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers={**auth_header(token), "X-Session-ID": "my-sid-123"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["sid"] == "my-sid-123"

    # Verify Redis
    raw = await redis_client.get(RedisKey.edit_lock(room.id))
    lock = json.loads(raw)
    assert lock["sid"] == "my-sid-123"


@pytest.mark.asyncio
async def test_acquire_edit_lock_broadcasts_lock_update(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio,
) -> None:
    """Test PUT broadcasts LockUpdate socket event with ttl."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "drawing"},
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) >= 1
    evt = mock_sio.emitted[-1]
    assert evt["event"] == "lock_update"
    assert evt["data"]["action"] == "acquired"
    assert evt["data"]["user_id"] == str(user.id)
    assert evt["data"]["msg"] == "drawing"
    assert evt["data"]["ttl"] is not None
    assert evt["data"]["ttl"] > 0


@pytest.mark.asyncio
async def test_refresh_edit_lock_with_token(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT with Lock-Token header refreshes the lock."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire
    resp1 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "first"},
        headers=auth_header(token),
    )
    assert resp1.status_code == 200
    lock_token = resp1.json()["lock_token"]

    # Refresh with Lock-Token
    resp2 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "updated"},
        headers={**auth_header(token), "Lock-Token": lock_token},
    )
    assert resp2.status_code == 200
    data = resp2.json()
    assert data["locked"] is True
    assert data["lock_token"] == lock_token
    assert data["user_id"] == str(user.id)


@pytest.mark.asyncio
async def test_refresh_with_expired_lock_returns_409(
    client: AsyncClient, session: AsyncSession, redis_client: Redis
) -> None:
    """Test refresh returns 409 when lock has expired."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire
    resp = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=auth_header(token),
    )
    lock_token = resp.json()["lock_token"]

    # Simulate expiry by deleting the key
    await redis_client.delete(RedisKey.edit_lock(room.id))

    # Try to refresh → 409
    resp2 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers={**auth_header(token), "Lock-Token": lock_token},
    )
    assert resp2.status_code == 409


@pytest.mark.asyncio
async def test_acquire_without_token_when_lock_exists_returns_423(
    client: AsyncClient, session: AsyncSession
) -> None:
    """PUT without Lock-Token when lock exists returns 423."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire
    resp1 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "first"},
        headers=auth_header(token),
    )
    assert resp1.status_code == 200

    # Same user, no Lock-Token → 423 (lock already exists)
    resp2 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "second attempt"},
        headers=auth_header(token),
    )
    assert resp2.status_code == 423


@pytest.mark.asyncio
async def test_acquire_edit_lock_conflict_with_other_user(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT returns 423 when another user holds the lock."""
    user1, token1 = await create_test_user_in_db(session, "user1@test")
    user2, token2 = await create_test_user_in_db(session, "user2@test")
    room = await create_test_room(session, user1)

    # Give user2 room membership
    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user2.id,  # type: ignore[arg-type]
        role=MemberRole.MEMBER,
    )
    session.add(membership)
    await session.commit()

    # User1 acquires lock
    resp1 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "user1 editing"},
        headers=auth_header(token1),
    )
    assert resp1.status_code == 200

    # User2 tries to acquire → 423
    resp2 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "user2 editing"},
        headers=auth_header(token2),
    )
    assert resp2.status_code == 423
    assert "locked" in resp2.json()["type"]


@pytest.mark.asyncio
async def test_acquire_edit_lock_blocked_by_admin_lock(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT returns 423 when room is admin-locked."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    room.locked = True
    session.add(room)
    await session.commit()

    response = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=auth_header(token),
    )
    assert response.status_code == 423


@pytest.mark.asyncio
async def test_lock_auto_expires_after_ttl(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Lock disappears from Redis after edit_lock_ttl without refresh."""
    user1, token1 = await create_test_user_in_db(session, "user1@test")
    user2, token2 = await create_test_user_in_db(session, "user2@test")
    room = await create_test_room(session, user1)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user2.id,  # type: ignore[arg-type]
        role=MemberRole.MEMBER,
    )
    session.add(membership)
    await session.commit()

    # User1 acquires
    resp = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "no refresh"},
        headers=auth_header(token1),
    )
    assert resp.status_code == 200

    # Wait for TTL + buffer (default 10s + 1s)
    await asyncio.sleep(11)

    # Lock should be gone
    status = await client.get(
        f"/v1/rooms/{room.id}/edit-lock", headers=auth_header(token1)
    )
    assert status.json()["locked"] is False

    # User2 can now acquire
    resp2 = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "after expiry"},
        headers=auth_header(token2),
    )
    assert resp2.status_code == 200
    assert resp2.json()["user_id"] == str(user2.id)


# =============================================================================
# DELETE /v1/rooms/{room_id}/edit-lock (Release)
# =============================================================================


@pytest.mark.asyncio
async def test_release_edit_lock_with_token(
    client: AsyncClient, session: AsyncSession, redis_client: Redis
) -> None:
    """Test DELETE with Lock-Token releases the lock."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire
    resp = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "temp"},
        headers=auth_header(token),
    )
    lock_token = resp.json()["lock_token"]

    # Release with Lock-Token
    response = await client.delete(
        f"/v1/rooms/{room.id}/edit-lock",
        headers={**auth_header(token), "Lock-Token": lock_token},
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify Redis entry removed
    raw = await redis_client.get(RedisKey.edit_lock(room.id))
    assert raw is None


@pytest.mark.asyncio
async def test_release_edit_lock_broadcasts_lock_update(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio,
) -> None:
    """Test DELETE broadcasts LockUpdate with action=released."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire then release
    resp = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=auth_header(token),
    )
    lock_token = resp.json()["lock_token"]
    mock_sio.emitted.clear()

    await client.delete(
        f"/v1/rooms/{room.id}/edit-lock",
        headers={**auth_header(token), "Lock-Token": lock_token},
    )

    assert len(mock_sio.emitted) >= 1
    evt = mock_sio.emitted[-1]
    assert evt["event"] == "lock_update"
    assert evt["data"]["action"] == "released"


@pytest.mark.asyncio
async def test_release_edit_lock_idempotent_when_no_lock(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test DELETE succeeds even if no lock exists."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=auth_header(token)
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())


@pytest.mark.asyncio
async def test_release_edit_lock_wrong_token_returns_403(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test DELETE with wrong Lock-Token returns 403."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire
    await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=auth_header(token),
    )

    # Release with wrong Lock-Token
    response = await client.delete(
        f"/v1/rooms/{room.id}/edit-lock",
        headers={**auth_header(token), "Lock-Token": "wrong-token"},
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_release_edit_lock_forbidden_for_non_holder(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test DELETE returns 403 when non-holder tries to release."""
    user1, token1 = await create_test_user_in_db(session, "user1@test")
    user2, token2 = await create_test_user_in_db(session, "user2@test")
    room = await create_test_room(session, user1)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user2.id,  # type: ignore[arg-type]
        role=MemberRole.MEMBER,
    )
    session.add(membership)
    await session.commit()

    # User1 acquires
    await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=auth_header(token1),
    )

    # User2 tries to release → 403
    response = await client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=auth_header(token2)
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_admin_can_release_any_lock(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test admin can release another user's lock."""
    user, token = await create_test_user_in_db(session, "user@test")
    admin, admin_token = await create_test_user_in_db(
        session, "admin@test", is_superuser=True
    )
    room = await create_test_room(session, user)

    # Give admin membership
    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=admin.id,  # type: ignore[arg-type]
        role=MemberRole.MEMBER,
    )
    session.add(membership)
    await session.commit()

    # User acquires
    await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={},
        headers=auth_header(token),
    )

    # Admin releases (no Lock-Token needed for admin)
    response = await client.delete(
        f"/v1/rooms/{room.id}/edit-lock", headers=auth_header(admin_token)
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())


# =============================================================================
# WritableRoomDep: Admin Lock Enforcement
# =============================================================================


@pytest.mark.asyncio
async def test_writable_room_blocks_non_admin_on_admin_locked_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test mutation returns 423 when room is admin-locked and user is not admin."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    room.locked = True
    session.add(room)
    await session.commit()

    # Try to set a bookmark (mutation endpoint using WritableRoomDep)
    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "test"},
        headers=auth_header(token),
    )
    assert response.status_code == 423
    assert "locked" in response.json()["type"]


@pytest.mark.asyncio
async def test_writable_room_allows_admin_on_admin_locked_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test admin can mutate even when room is admin-locked."""
    admin, token = await create_test_user_in_db(session, is_superuser=True)
    room = await create_test_room(session, admin)
    room.locked = True
    session.add(room)
    await session.commit()

    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "admin edit"},
        headers=auth_header(token),
    )
    assert response.status_code == 200


# =============================================================================
# WritableRoomDep: Edit Lock Enforcement with Lock-Token
# =============================================================================


@pytest.mark.asyncio
async def test_writable_room_blocks_non_holder(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test mutation returns 423 when another user holds the edit lock."""
    user1, token1 = await create_test_user_in_db(session, "user1@test")
    user2, token2 = await create_test_user_in_db(session, "user2@test")
    room = await create_test_room(session, user1)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user2.id,  # type: ignore[arg-type]
        role=MemberRole.MEMBER,
    )
    session.add(membership)
    await session.commit()

    # User1 acquires lock
    await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=auth_header(token1),
    )

    # User2 tries mutation → 423
    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "blocked"},
        headers=auth_header(token2),
    )
    assert response.status_code == 423
    assert "locked" in response.json()["type"]


@pytest.mark.asyncio
async def test_writable_room_allows_lock_holder_with_token(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test lock holder can mutate when sending Lock-Token header."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire lock
    resp = await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=auth_header(token),
    )
    lock_token = resp.json()["lock_token"]

    # Holder can mutate with Lock-Token
    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "allowed"},
        headers={**auth_header(token), "Lock-Token": lock_token},
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_writable_room_blocks_wrong_lock_token(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test mutation returns 423 when wrong Lock-Token is sent."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Acquire lock
    await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=auth_header(token),
    )

    # Wrong Lock-Token → 423
    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/0",
        json={"label": "blocked"},
        headers={**auth_header(token), "Lock-Token": "wrong-token"},
    )
    assert response.status_code == 423


@pytest.mark.asyncio
async def test_writable_room_allows_get_when_locked(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET endpoints still work when room has edit lock."""
    user1, token1 = await create_test_user_in_db(session, "user1@test")
    user2, token2 = await create_test_user_in_db(session, "user2@test")
    room = await create_test_room(session, user1)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user2.id,  # type: ignore[arg-type]
        role=MemberRole.MEMBER,
    )
    session.add(membership)
    await session.commit()

    # User1 acquires lock
    await client.put(
        f"/v1/rooms/{room.id}/edit-lock",
        json={"msg": "editing"},
        headers=auth_header(token1),
    )

    # User2 can still read bookmarks
    response = await client.get(
        f"/v1/rooms/{room.id}/bookmarks",
        headers=auth_header(token2),
    )
    assert response.status_code == 200
