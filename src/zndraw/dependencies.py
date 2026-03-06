"""FastAPI dependencies for database sessions, Redis, and authentication.

All resources are accessed from request.app.state following Flask/Celery pattern.
Authentication uses zndraw-auth package.
"""

import json
from pathlib import Path as FilePath
from typing import Annotated, NamedTuple

from fastapi import Depends, Path, Request
from redis.asyncio import Redis as AsyncRedis
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker
from zndraw_auth import (
    SessionDep,
    User,
    current_active_user,
    current_optional_user,
    current_superuser,
    current_user_scoped_session,
)
from zndraw_auth.db import get_session_maker
from zndraw_joblib.dependencies import ResultBackend, validate_room_id
from zndraw_joblib.settings import JobLibSettings
from zndraw_socketio import AsyncServerWrapper

from zndraw.exceptions import Forbidden, RoomLocked, RoomNotFound, SessionNotFound
from zndraw.geometries import geometries as geometry_models
from zndraw.geometries.camera import Camera
from zndraw.models import Room, RoomGeometry
from zndraw.redis import RedisKey
from zndraw.storage.base import StorageBackend
from zndraw.storage.router import StorageRouter

# Re-export auth dependencies for convenience
CurrentUserDep = Annotated[User, Depends(current_active_user)]
AdminUserDep = Annotated[User, Depends(current_superuser)]
OptionalUserDep = Annotated[User | None, Depends(current_optional_user)]

# Scoped-session variants — session closed before endpoint body runs,
# so the SQLite asyncio.Lock is NOT held during long-polling.
SessionMakerDep = Annotated[
    async_sessionmaker[AsyncSession], Depends(get_session_maker)
]
CurrentUserFactoryDep = Annotated[User, Depends(current_user_scoped_session)]


def get_redis(request: Request) -> AsyncRedis:  # type: ignore[type-arg]
    """Get the async Redis client from app.state."""
    return request.app.state.redis


RedisDep = Annotated[AsyncRedis, Depends(get_redis)]  # type: ignore[type-arg]


def get_storage(request: Request) -> StorageBackend:
    """Get frame storage backend from app.state."""
    return request.app.state.frame_storage


StorageDep = Annotated[StorageBackend, Depends(get_storage)]


def get_storage_router(request: Request) -> StorageRouter:
    """Get StorageRouter for mount management."""
    storage = request.app.state.frame_storage
    if not isinstance(storage, StorageRouter):
        raise RuntimeError(f"Expected StorageRouter, got {type(storage).__name__}")
    return storage


StorageRouterDep = Annotated[StorageRouter, Depends(get_storage_router)]


def get_tsio(request: Request) -> AsyncServerWrapper:
    """Get the zndraw-socketio typed wrapper from app.state."""
    return request.app.state.tsio


SioDep = Annotated[AsyncServerWrapper, Depends(get_tsio)]


def get_result_backend(request: Request) -> ResultBackend:
    """Get ResultBackend from app.state."""
    return request.app.state.result_backend


ResultBackendDep = Annotated[ResultBackend, Depends(get_result_backend)]


def get_joblib_settings(request: Request) -> JobLibSettings:
    """Get JobLibSettings from app.state."""
    return request.app.state.joblib_settings


JobLibSettingsDep = Annotated[JobLibSettings, Depends(get_joblib_settings)]


def get_media_path(request: Request) -> FilePath:
    """Get media path from app.state.settings."""
    return request.app.state.settings.media_path


MediaPathDep = Annotated[FilePath, Depends(get_media_path)]


async def verify_room(session: AsyncSession, room_id: str) -> Room:
    """Verify room exists and return it, or raise RoomNotFound."""
    room = await session.get(Room, room_id)
    if room is None:
        raise RoomNotFound.exception(f"Room with id {room_id} not found")
    return room


def room_channel(room_id: str) -> str:
    """Get Socket.IO room channel name."""
    return f"room:{room_id}"


async def _check_locks(
    redis: AsyncRedis,
    room_id: str,
    room: Room,
    user: User,  # type: ignore[type-arg]
    lock_token: str | None = None,
) -> None:
    """Check admin lock and edit lock; raise RoomLocked (423) if blocked."""
    if room.locked and not user.is_superuser:
        raise RoomLocked.exception("Room is locked by an administrator")

    raw = await redis.get(RedisKey.edit_lock(room_id))
    if raw is not None:
        holder = json.loads(raw)
        if lock_token is not None:
            if holder["lock_token"] != lock_token:
                raise RoomLocked.exception("Room is being edited by another session")
        else:
            if holder["user_id"] != str(user.id):
                raise RoomLocked.exception("Room is being edited by another user")


async def get_writable_room(
    request: Request,
    session: SessionDep,
    current_user: CurrentUserDep,
    redis: RedisDep,
    room_id: str = Path(),
) -> Room:
    """Verify room exists and is writable by the current user.

    Validates room_id format, checks admin lock (SQL) and edit lock (Redis).
    Raises RoomLocked (423) if the room cannot be modified.
    """
    validate_room_id(room_id)
    room = await verify_room(session, room_id)
    lock_token = request.headers.get("Lock-Token")
    await _check_locks(redis, room_id, room, current_user, lock_token)
    return room


WritableRoomDep = Annotated[Room, Depends(get_writable_room)]


async def get_verified_session_id(
    current_user: CurrentUserDep,
    redis: RedisDep,
    room_id: str = Path(),
    session_id: str = Path(),
) -> str:
    """Verify session belongs to the current user, or raise 404.

    Checks the presence key to confirm the session exists and is owned
    by the requesting user. Returns 404 for both missing and non-owned
    sessions (prevents enumeration).
    """
    user_id = await redis.get(  # type: ignore[misc]
        RedisKey.presence_sid(room_id, session_id)
    )
    if user_id is None or user_id != str(current_user.id):
        raise SessionNotFound.exception("Session not found")
    return session_id


VerifiedSessionDep = Annotated[str, Depends(get_verified_session_id)]


async def get_owner_from_geometry(
    redis: AsyncRedis,  # type: ignore[type-arg]
    session: AsyncSession,
    room_id: str,
    key: str,
) -> str | None:
    """Read owner from geometry config via Pydantic validation.

    Tries Redis hash first (session cameras), then SQL.
    Returns None if geometry doesn't exist or has no owner.
    """
    raw = await redis.hget(RedisKey.room_cameras(room_id), key)  # type: ignore[misc]
    if raw is not None:
        entry = json.loads(raw)
        return Camera(**entry["data"]).owner

    row = await session.get(RoomGeometry, (room_id, key))
    if row is not None:
        model_cls = geometry_models.get(row.type)
        if model_cls is not None:
            return model_cls(**json.loads(row.config)).owner

    return None


class WritableGeometryInfo(NamedTuple):
    """Resolved room and current owner for a writable geometry."""

    room: Room
    current_owner: str | None


async def check_geometry_write_access(
    session: AsyncSession,
    redis: AsyncRedis,  # type: ignore[type-arg]
    room_id: str,
    geometry_key: str,
    current_user: User,
    lock_token: str | None = None,
) -> WritableGeometryInfo:
    """Check that the current user can write to a geometry.

    Combines lock checks with ownership: owners can edit their own geometries
    even when the room is admin-locked.

    Lock precedence
    ---------------
    1. Edit lock: blocks non-holders (even owners)
    2. Admin lock: blocks non-owners (owners bypass)
    3. Ownership: blocks non-owners always
    """
    validate_room_id(room_id)
    room = await verify_room(session, room_id)

    # Edit lock: always blocks non-holders (even owners)
    raw = await redis.get(RedisKey.edit_lock(room_id))
    if raw is not None:
        holder = json.loads(raw)
        if lock_token is not None:
            if holder["lock_token"] != lock_token:
                raise RoomLocked.exception("Room is being edited by another session")
        else:
            if holder["user_id"] != str(current_user.id):
                raise RoomLocked.exception("Room is being edited by another user")

    # Resolve owner once — reused for lock bypass and ownership check
    current_owner = await get_owner_from_geometry(redis, session, room_id, geometry_key)
    user_id_str = str(current_user.id)
    is_unowned = current_owner is None
    is_mine = current_owner == user_id_str

    # Admin lock: only owners of their OWN geometry bypass
    if room.locked and not current_user.is_superuser and not is_mine:
        raise RoomLocked.exception("Room is locked by an administrator")

    # Ownership check: non-owners blocked even without locks
    if not current_user.is_superuser and not is_unowned and not is_mine:
        raise Forbidden.exception("Not the geometry owner")

    return WritableGeometryInfo(room=room, current_owner=current_owner)


async def get_writable_geometry(
    request: Request,
    session: SessionDep,
    current_user: CurrentUserDep,
    redis: RedisDep,
    room_id: str = Path(),
    key: str = Path(),
) -> WritableGeometryInfo:
    """FastAPI dependency wrapping check_geometry_write_access."""
    lock_token = request.headers.get("Lock-Token")
    return await check_geometry_write_access(
        session, redis, room_id, key, current_user, lock_token
    )


WritableGeometryDep = Annotated[WritableGeometryInfo, Depends(get_writable_geometry)]


async def get_writable_room_id(
    request: Request,
    session: SessionDep,
    current_user: CurrentUserDep,
    redis: RedisDep,
    room_id: str = Path(),
) -> str:
    """Verify room is writable and return its id.

    Override for zndraw-joblib's ``verify_writable_room``.
    Same checks as ``get_writable_room`` but returns ``str``.
    Virtual rooms (@global, @internal) skip DB/lock checks since they
    are not stored as Room records — access is controlled by joblib's
    own admin check in the registration endpoint.
    """
    validate_room_id(room_id)
    if room_id not in ("@global", "@internal"):
        room = await verify_room(session, room_id)
        lock_token = request.headers.get("Lock-Token")
        await _check_locks(redis, room_id, room, current_user, lock_token)
    return room_id
