"""Edit lock REST API endpoints for room-level exclusive editing."""

import json
import time

from fastapi import APIRouter, status

from zndraw.dependencies import (
    CurrentUserDep,
    RedisDep,
    SessionDep,
    SioDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    Forbidden,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    problem_responses,
)
from zndraw.redis import RedisKey
from zndraw.schemas import (
    EDIT_LOCK_TTL,
    EditLockRequest,
    EditLockResponse,
    StatusResponse,
)
from zndraw.socket_events import LockUpdate

router = APIRouter(prefix="/v1/rooms/{room_id}/edit-lock", tags=["edit-lock"])


async def _read_lock(redis: RedisDep, room_id: str) -> EditLockResponse:
    """Read the current edit lock state from Redis."""
    raw = await redis.get(RedisKey.edit_lock(room_id))
    if raw is None:
        return EditLockResponse(locked=False)
    data = json.loads(raw)
    ttl = await redis.ttl(RedisKey.edit_lock(room_id))
    return EditLockResponse(
        locked=True,
        user_id=data["user_id"],
        msg=data.get("msg"),
        acquired_at=data["acquired_at"],
        ttl=max(ttl, 0),
    )


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def get_edit_lock(
    session: SessionDep,
    redis: RedisDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> EditLockResponse:
    """Get current edit lock status for a room."""
    await verify_room(session, room_id)
    return await _read_lock(redis, room_id)


@router.put(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound, RoomLocked),
)
async def acquire_edit_lock(
    session: SessionDep,
    redis: RedisDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    room_id: str,
    request: EditLockRequest,
) -> EditLockResponse:
    """Acquire or refresh the room edit lock (idempotent for the holder)."""
    room = await verify_room(session, room_id)

    # Admin lock blocks non-admins
    if room.locked and not current_user.is_superuser:
        raise RoomLocked.exception("Room is locked by an administrator")

    user_id = str(current_user.id)
    key = RedisKey.edit_lock(room_id)

    # Check existing lock
    raw = await redis.get(key)
    if raw is not None:
        holder = json.loads(raw)
        if holder["user_id"] != user_id:
            raise RoomLocked.exception("Room is being edited by another user")
        # Refresh: keep original acquired_at
        acquired_at = holder["acquired_at"]
        action = "refreshed"
    else:
        acquired_at = time.time()
        action = "acquired"

    # Set/refresh lock with TTL
    lock_data = json.dumps(
        {"user_id": user_id, "msg": request.msg, "acquired_at": acquired_at}
    )
    await redis.set(key, lock_data, ex=EDIT_LOCK_TTL)

    await sio.emit(
        LockUpdate(
            room_id=room_id,
            action=action,
            user_id=user_id,
            msg=request.msg,
        ),
        room=room_channel(room_id),
    )

    ttl = await redis.ttl(key)
    return EditLockResponse(
        locked=True,
        user_id=user_id,
        msg=request.msg,
        acquired_at=acquired_at,
        ttl=max(ttl, 0),
    )


@router.delete(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound, Forbidden),
)
async def release_edit_lock(
    session: SessionDep,
    redis: RedisDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> StatusResponse:
    """Release the room edit lock."""
    await verify_room(session, room_id)

    key = RedisKey.edit_lock(room_id)
    raw = await redis.get(key)

    # No lock held — idempotent success
    if raw is None:
        return StatusResponse()

    holder = json.loads(raw)
    user_id = str(current_user.id)

    # Only the holder or an admin can release
    if holder["user_id"] != user_id and not current_user.is_superuser:
        raise Forbidden.exception("Only the lock holder or an admin can release")

    await redis.delete(key)

    await sio.emit(
        LockUpdate(
            room_id=room_id,
            action="released",
            user_id=user_id,
        ),
        room=room_channel(room_id),
    )

    return StatusResponse()
