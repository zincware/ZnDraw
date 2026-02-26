"""Socket.IO server with room support and presence tracking."""

import json
from uuid import UUID

import socketio
from fastapi import Depends
from fastapi_users.jwt import decode_jwt
from jwt import InvalidTokenError
from sqlalchemy import select as sa_select
from sqlmodel import and_, select
from zndraw_auth import AuthSettings, SessionDep, User, get_auth_settings
from zndraw_joblib import Worker, cleanup_worker
from zndraw_joblib.events import emit as joblib_emit
from zndraw_joblib.models import ProviderRecord
from zndraw_socketio import EventContext, wrap

from zndraw.config import Settings, get_zndraw_settings
from zndraw.dependencies import RedisDep, StorageDep, room_channel
from zndraw.exceptions import (
    NotInRoom,
    NotRoomMember,
    ProblemException,
    RoomNotFound,
    UserNotFound,
)
from zndraw.geometries.camera import Camera
from zndraw.models import Room, RoomMembership
from zndraw.redis import RedisKey
from zndraw.schemas import ProgressResponse
from zndraw.socket_events import (
    FramesInvalidate,
    GeometryInvalidate,
    Heartbeat,
    HeartbeatResponse,
    LockUpdate,
    RoomJoin,
    RoomJoinResponse,
    RoomLeave,
    RoomLeaveResponse,
    SessionJoined,
    SessionLeft,
    Typing,
    TypingResponse,
    TypingStart,
    TypingStop,
    UserGet,
    UserGetResponse,
)

# Module-level Socket.IO server wrapped with zndraw-socketio (drop-in replacement)
# tsio.app is set in database.py lifespan to enable DI (resolved at event time)
tsio = wrap(socketio.AsyncServer(async_mode="asgi", cors_allowed_origins="*"))


async def _cleanup_session(
    redis: RedisDep, sid: str, sio_session: dict, room_id: str
) -> None:
    """Remove this SID's camera and settings from the room and broadcast deletion."""
    # Remove active camera tracking
    await redis.hdel(RedisKey.active_cameras(room_id), sid)  # type: ignore[misc]

    # Remove session settings
    await redis.delete(RedisKey.session_settings(room_id, sid))

    camera_key: str | None = sio_session.get("camera_key")
    if camera_key is None:
        return
    hash_key = RedisKey.room_cameras(room_id)
    deleted = await redis.hdel(hash_key, camera_key)  # type: ignore[misc]
    if deleted:
        await tsio.emit(
            GeometryInvalidate(room_id=room_id, operation="delete", key=camera_key),
            room=room_channel(room_id),
        )


# =============================================================================
# Exception Handler for RFC 9457 Errors
# =============================================================================


@tsio.exception_handler(ProblemException)
async def handle_problem(ctx: EventContext, exc: ProblemException) -> dict:
    """Convert ProblemException to RFC 9457 response dict."""
    return exc.problem.model_dump(exclude_none=True)


# =============================================================================
# Connection Lifecycle Handlers
# =============================================================================


@tsio.on("connect")
async def on_connect(
    sid: str,
    environ: dict,
    auth: dict | None = None,
    auth_settings: AuthSettings = Depends(get_auth_settings),
) -> bool:
    """Handle Socket.IO connection with JWT validation."""
    token = auth.get("token") if isinstance(auth, dict) else None

    if token is None:
        raise ConnectionRefusedError("No authentication token provided")

    try:
        payload = decode_jwt(
            token,
            auth_settings.secret_key.get_secret_value(),
            audience=["fastapi-users:auth"],
        )
    except InvalidTokenError as e:
        raise ConnectionRefusedError(f"Invalid token: {e}") from None

    user_id = UUID(payload["sub"])
    await tsio.save_session(sid, {"user_id": user_id, "current_room_id": None})
    await tsio.enter_room(sid, f"user:{user_id}")
    return True


@tsio.on(
    "disconnect",
    emits=[GeometryInvalidate, LockUpdate, SessionLeft, FramesInvalidate],
)
async def on_disconnect(
    sid: str,
    reason: str,
    redis: RedisDep,
    db_session: SessionDep,
) -> None:
    """Handle disconnect - clean up presence, camera, locks, and workers."""
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")

    if current_room_id is not None:
        await redis.delete(RedisKey.presence_sid(current_room_id, sid))

        # Delete session camera from room hash
        await _cleanup_session(redis, sid, sio_session, current_room_id)

        await tsio.emit(
            SessionLeft(room_id=current_room_id, user_id=user_id, sid=sid),
            room=room_channel(current_room_id),
        )

        # Release edit lock if this disconnecting session holds it
        lock_key = RedisKey.edit_lock(current_room_id)
        raw_lock = await redis.get(lock_key)
        if raw_lock is not None:
            holder = json.loads(raw_lock)
            if holder.get("sid") == sid:
                await redis.delete(lock_key)
                await tsio.emit(
                    LockUpdate(
                        room_id=current_room_id,
                        action="released",
                        user_id=str(user_id),
                        sid=sid,
                    ),
                    room=room_channel(current_room_id),
                )

    # Worker cleanup — fail tasks, remove links, soft-delete orphan jobs
    worker_id = sio_session.get("worker_id")
    if worker_id:
        worker = await db_session.get(Worker, UUID(worker_id))
        if worker:
            # Find rooms with frames providers BEFORE cleanup deletes them
            result = await db_session.execute(
                sa_select(ProviderRecord.room_id).where(
                    ProviderRecord.worker_id == worker.id,
                    ProviderRecord.category == "frames",
                )
            )
            frame_provider_rooms = result.scalars().all()

            emissions = await cleanup_worker(db_session, worker)
            await db_session.commit()
            await joblib_emit(tsio, emissions)

            # Clear provider frame counts and notify frontends
            for rid in frame_provider_rooms:
                await redis.delete(RedisKey.provider_frame_count(rid))  # type: ignore[misc]
                await tsio.emit(
                    FramesInvalidate(
                        room_id=rid,
                        action="clear",
                        count=0,
                        reason="provider_disconnected",
                    ),
                    room=room_channel(rid),
                )


# =============================================================================
# Event Handlers - Use Depends pattern (same deps as FastAPI routes!)
# =============================================================================


@tsio.on(UserGet)
async def user_get(sid: str, data: UserGet, session: SessionDep) -> UserGetResponse:
    """Return the authenticated user's information."""
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]

    user = await session.get(User, user_id)
    if user is None:
        raise UserNotFound.exception("User not found")

    return UserGetResponse(id=user.id, email=user.email, is_superuser=user.is_superuser)


@tsio.on(RoomJoin, emits=[SessionLeft, SessionJoined, GeometryInvalidate])
async def room_join(
    sid: str,
    data: RoomJoin,
    redis: RedisDep,
    storage: StorageDep,
    session: SessionDep,
    settings: Settings = Depends(get_zndraw_settings),
) -> RoomJoinResponse:
    """Join a Socket.IO room for real-time updates.

    Supports special system rooms with '@' prefix (e.g., @overview) that skip
    database validation, presence tracking, and camera creation.
    """
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]

    # Check if this is a special system room (@ prefix)
    is_system_room = data.room_id.startswith("@")

    # Validate and load room from database (skip for system rooms)
    room_locked = False
    email = None
    room: Room | None = None
    if not is_system_room:
        room = await session.get(Room, data.room_id)
        if room is None:
            raise RoomNotFound.exception(f"Room with id {data.room_id} not found")

        room_locked = room.locked

        result = await session.exec(  # type: ignore[attr-defined]
            select(RoomMembership).where(
                and_(
                    RoomMembership.room_id == data.room_id,
                    RoomMembership.user_id == user_id,
                )
            )
        )
        membership = result.first()

        if membership is None and not room.is_public:
            raise NotRoomMember.exception("Not a member of this private room")

    user = await session.get(User, user_id)
    email = user.email if user else None

    # Leave previous room if any
    old_room_id: str | None = sio_session.get("current_room_id")
    if old_room_id is not None:
        await tsio.leave_room(sid, room_channel(old_room_id))
        # Clean up presence and camera for non-system rooms
        if not old_room_id.startswith("@"):
            await redis.delete(RedisKey.presence_sid(old_room_id, sid))
            await _cleanup_session(redis, sid, sio_session, old_room_id)
        await tsio.emit(
            SessionLeft(room_id=old_room_id, user_id=user_id, sid=sid),
            room=room_channel(old_room_id),
        )

    # Join new room
    await tsio.enter_room(sid, room_channel(data.room_id))
    sio_session["current_room_id"] = data.room_id
    sio_session["client_type"] = data.client_type
    await tsio.save_session(sid, sio_session)

    if data.client_type == "frontend":
        await tsio.enter_room(sid, "frontend")
        await tsio.enter_room(sid, room_channel("@global"))

    # System rooms (@-prefixed) skip presence, camera, and return minimal data
    if is_system_room:
        return RoomJoinResponse(
            room_id=data.room_id,
            session_id=sid,
            step=0,
            frame_count=0,
            locked=False,
        )

    assert room is not None  # is_system_room is False → room was assigned above

    # Regular rooms: set presence, create camera, broadcast join
    await redis.set(
        RedisKey.presence_sid(data.room_id, sid),
        str(user_id),
        ex=settings.presence_ttl,
    )

    await tsio.emit(
        SessionJoined(room_id=data.room_id, user_id=user_id, sid=sid, email=email),
        room=room_channel(data.room_id),
        skip_sid=sid,
    )

    # Capture room state before potential commit
    room_step = room.step

    # Create session camera in Redis hash (frontend only — pyclients don't need a viewport)
    camera_key: str | None = None
    if data.client_type == "frontend":
        camera_key = f"cam:{email}:{sid[:8]}"
        camera = Camera(owner=str(user_id))
        camera_value = json.dumps(
            {
                "sid": sid,
                "email": email,
                "data": camera.model_dump(),
            }
        )
        hash_key = RedisKey.room_cameras(data.room_id)
        await redis.hset(hash_key, camera_key, camera_value)  # type: ignore[misc]

        # Store camera_key in session for cleanup on leave/disconnect
        sio_session["camera_key"] = camera_key
        await tsio.save_session(sid, sio_session)

        # Track as active camera for deletion prevention
        await redis.hset(RedisKey.active_cameras(data.room_id), sid, camera_key)  # type: ignore[misc]

        await tsio.emit(
            GeometryInvalidate(room_id=data.room_id, operation="set", key=camera_key),
            room=room_channel(data.room_id),
        )

    frame_count = await storage.get_length(data.room_id)

    # Load active progress trackers from Redis
    progress_raw = await redis.hgetall(RedisKey.room_progress(data.room_id))  # type: ignore[misc]
    progress_trackers = {
        pid: ProgressResponse(**json.loads(pdata))
        for pid, pdata in progress_raw.items()
    }

    return RoomJoinResponse(
        room_id=data.room_id,
        session_id=sid,
        step=room_step,
        frame_count=frame_count,
        locked=room_locked,
        camera_key=camera_key,
        progress_trackers=progress_trackers,
    )


@tsio.on(RoomLeave, emits=[SessionLeft, GeometryInvalidate])
async def room_leave(sid: str, data: RoomLeave, redis: RedisDep) -> RoomLeaveResponse:
    """Leave a Socket.IO room.

    Idempotent: if already in a different room (room switch race), succeeds silently.
    Raises NotInRoom if not in any room (helps catch bugs).
    """
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")

    # Not in any room - this is likely a bug
    if current_room_id is None:
        raise NotInRoom.exception("Not currently in a room")

    # Already in a different room - room_join already handled the leave
    # This happens during room switching when cleanup runs after join
    if current_room_id != data.room_id:
        return RoomLeaveResponse(room_id=data.room_id)

    await tsio.leave_room(sid, room_channel(data.room_id))
    await redis.delete(RedisKey.presence_sid(data.room_id, sid))

    # Clean up camera
    await _cleanup_session(redis, sid, sio_session, data.room_id)

    sio_session["current_room_id"] = None
    await tsio.save_session(sid, sio_session)

    await tsio.emit(
        SessionLeft(room_id=data.room_id, user_id=user_id, sid=sid),
        room=room_channel(data.room_id),
    )
    return RoomLeaveResponse(room_id=data.room_id)


@tsio.on(Heartbeat)
async def heartbeat(
    sid: str,
    data: Heartbeat,
    redis: RedisDep,
    settings: Settings = Depends(get_zndraw_settings),
) -> HeartbeatResponse:
    """Keep presence alive (should be called every 30s).

    Camera cleanup is explicit (on leave/disconnect), so no TTL refresh needed.
    """
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")

    if current_room_id != data.room_id:
        raise NotInRoom.exception("Not in this room")

    await redis.set(
        RedisKey.presence_sid(data.room_id, sid),
        str(user_id),
        ex=settings.presence_ttl,
    )

    # Refresh settings TTL alongside presence
    await redis.expire(  # type: ignore[misc]
        RedisKey.session_settings(data.room_id, sid),
        settings.presence_ttl,
    )

    return HeartbeatResponse()


@tsio.on(TypingStart, emits=[Typing])
async def typing_start(
    sid: str, data: TypingStart, session: SessionDep
) -> TypingResponse:
    """Broadcast that user started typing."""
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")

    if current_room_id != data.room_id:
        raise NotInRoom.exception("Not in this room")

    user = await session.get(User, user_id)
    email = user.email if user else None

    await tsio.emit(
        Typing(room_id=data.room_id, user_id=user_id, email=email, is_typing=True),
        room=room_channel(data.room_id),
        skip_sid=sid,
    )
    return TypingResponse()


@tsio.on(TypingStop, emits=[Typing])
async def typing_stop(
    sid: str, data: TypingStop, session: SessionDep
) -> TypingResponse:
    """Broadcast that user stopped typing."""
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")

    if current_room_id != data.room_id:
        raise NotInRoom.exception("Not in this room")

    user = await session.get(User, user_id)
    email = user.email if user else None

    await tsio.emit(
        Typing(room_id=data.room_id, user_id=user_id, email=email, is_typing=False),
        room=room_channel(data.room_id),
        skip_sid=sid,
    )
    return TypingResponse()
