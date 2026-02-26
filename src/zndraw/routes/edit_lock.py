"""Edit lock REST API endpoints for room-level exclusive editing."""

import json
import time
import uuid
from typing import Annotated

from fastapi import APIRouter, Header

from zndraw.config import SettingsDep
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
    LockExpired,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    problem_responses,
)
from zndraw.redis import RedisKey
from zndraw.schemas import (
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
        lock_token=data["lock_token"],
        user_id=data["user_id"],
        sid=data.get("sid"),
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
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, LockExpired
    ),
)
async def acquire_edit_lock(
    session: SessionDep,
    redis: RedisDep,
    sio: SioDep,
    settings: SettingsDep,
    current_user: CurrentUserDep,
    room_id: str,
    request: EditLockRequest,
    lock_token: Annotated[str | None, Header(alias="Lock-Token")] = None,
    x_session_id: Annotated[str | None, Header(alias="X-Session-ID")] = None,
) -> EditLockResponse:
    """Acquire or refresh the room edit lock.

    - Without Lock-Token header: acquire a new lock (generates token).
    - With Lock-Token header: refresh an existing lock (resets TTL).
    - Returns 409 if Lock-Token is provided but the lock has expired.
    - Returns 423 if another session holds the lock.
    """
    room = await verify_room(session, room_id)

    if room.locked and not current_user.is_superuser:
        raise RoomLocked.exception("Room is locked by an administrator")

    user_id = str(current_user.id)
    key = RedisKey.edit_lock(room_id)
    raw = await redis.get(key)

    if lock_token is not None:
        # --- Refresh path: client sends Lock-Token header ---
        if raw is None:
            raise LockExpired.exception("Lock expired, re-acquire required")
        holder = json.loads(raw)
        if holder["lock_token"] != lock_token:
            raise RoomLocked.exception("Room is being edited by another session")
        # Refresh: keep all original fields, reset TTL
        await redis.set(key, raw, ex=settings.edit_lock_ttl)
        ttl = await redis.ttl(key)
        return EditLockResponse(
            locked=True,
            lock_token=holder["lock_token"],
            user_id=holder["user_id"],
            sid=holder.get("sid"),
            msg=request.msg if request.msg is not None else holder.get("msg"),
            acquired_at=holder["acquired_at"],
            ttl=max(ttl, 0),
        )

    # --- Acquire path: no Lock-Token header ---
    new_token = str(uuid.uuid4())
    acquired_at = time.time()
    lock_data = json.dumps(
        {
            "lock_token": new_token,
            "user_id": user_id,
            "sid": x_session_id,
            "msg": request.msg,
            "acquired_at": acquired_at,
        }
    )
    # Atomic set-if-not-exists — prevents TOCTOU race between concurrent acquires
    was_set = await redis.set(key, lock_data, ex=settings.edit_lock_ttl, nx=True)
    if not was_set:
        raise RoomLocked.exception("Room is being edited by another session")

    ttl = await redis.ttl(key)
    await sio.emit(
        LockUpdate(
            room_id=room_id,
            action="acquired",
            user_id=user_id,
            sid=x_session_id,
            msg=request.msg,
            ttl=max(ttl, 0),
        ),
        room=room_channel(room_id),
    )

    return EditLockResponse(
        locked=True,
        lock_token=new_token,
        user_id=user_id,
        sid=x_session_id,
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
    lock_token: Annotated[str | None, Header(alias="Lock-Token")] = None,
) -> StatusResponse:
    """Release the room edit lock."""
    await verify_room(session, room_id)

    key = RedisKey.edit_lock(room_id)
    raw = await redis.get(key)

    if raw is None:
        return StatusResponse()

    holder = json.loads(raw)

    # Admin can always release
    if not current_user.is_superuser:
        if lock_token is not None:
            if holder["lock_token"] != lock_token:
                raise Forbidden.exception("Lock token does not match")
        else:
            # Fallback: user_id check for backwards compat during migration
            if holder["user_id"] != str(current_user.id):
                raise Forbidden.exception(
                    "Only the lock holder or an admin can release"
                )

    await redis.delete(key)

    await sio.emit(
        LockUpdate(
            room_id=room_id,
            action="released",
            user_id=holder["user_id"],
            sid=holder.get("sid"),
        ),
        room=room_channel(room_id),
    )

    return StatusResponse()
