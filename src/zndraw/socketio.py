"""Socket.IO server with room support."""

import json
from uuid import UUID

import socketio
from fastapi import Depends
from fastapi_users.jwt import decode_jwt
from jwt import InvalidTokenError
from zndraw_socketio import EventContext, wrap

from zndraw.broadcast import room_channel
from zndraw.dependencies import FrameStorageDep, RedisDep
from zndraw.exceptions import (
    NotInRoom,
    ProblemError,
    RoomNotFound,
    UserNotFound,
)
from zndraw.geometries.camera import Camera
from zndraw.models import RoomGeometry
from zndraw.redis import RedisKey
from zndraw.schemas import ProgressResponse
from zndraw.socket_events import (
    GeometryInvalidate,
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
from zndraw_auth import AuthSettings, SessionDep, User, get_auth_settings

# Module-level Socket.IO server wrapped with zndraw-socketio (drop-in replacement)
# tsio.app is set in database.py lifespan to enable DI (resolved at event time)
tsio = wrap(socketio.AsyncServer(async_mode="asgi", cors_allowed_origins="*"))


async def _cleanup_session(
    redis: RedisDep, sid: str, sio_session: dict, room_id: str
) -> None:
    """Remove this SID's camera from the room and broadcast deletion."""
    # Remove active camera tracking
    await redis.hdel(RedisKey.active_cameras(room_id), sid)  # type: ignore[misc]

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


@tsio.exception_handler(ProblemError)
async def handle_problem(_ctx: EventContext, exc: ProblemError) -> dict:
    """Convert ProblemError to RFC 9457 response dict."""
    return exc.problem.model_dump(exclude_none=True)


# =============================================================================
# Connection Lifecycle Handlers
# =============================================================================


@tsio.on("connect")
async def on_connect(
    sid: str,
    _environ: dict,
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
    share_token = auth.get("share_token") if isinstance(auth, dict) else None
    await tsio.save_session(
        sid,
        {
            "user_id": user_id,
            "current_room_id": None,
            "current_room_address": None,
            "share_token": share_token,
        },
    )
    await tsio.enter_room(sid, f"user:{user_id}")
    await tsio.enter_room(sid, "rooms:feed")
    return True


@tsio.on(
    "disconnect",
    emits=[GeometryInvalidate, LockUpdate, SessionLeft],
)
async def on_disconnect(
    sid: str,
    _reason: str,
    redis: RedisDep,
) -> None:
    """Handle disconnect - clean up camera and locks."""
    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")
    current_room_address: str | None = sio_session.get("current_room_address")

    if current_room_id is not None:
        # Delete session camera from room hash
        await _cleanup_session(redis, sid, sio_session, current_room_id)

        await tsio.emit(
            SessionLeft(
                room_id=current_room_address or current_room_id,
                user_id=user_id,
                sid=sid,
            ),
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


# =============================================================================
# Event Handlers - Use Depends pattern (same deps as FastAPI routes!)
# =============================================================================


@tsio.on(UserGet)
async def user_get(sid: str, _data: UserGet, session: SessionDep) -> UserGetResponse:
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
    storage: FrameStorageDep,
    session: SessionDep,
) -> RoomJoinResponse:
    """Join a Socket.IO room for real-time updates."""
    from zndraw.access import can_read
    from zndraw.dependencies import (
        _load_room_by_address,
        fetch_group_role,
        resolve_share_token,
    )

    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    share_token: str | None = sio_session.get("share_token")

    composed = f"{data.owner_id}/{data.room_name}"

    room = await _load_room_by_address(session, data.owner_id, data.room_name)
    if room is None:
        raise RoomNotFound.exception(f"Room {composed} not found")

    user = await session.get(User, user_id)
    if user is None:
        raise UserNotFound.exception("User not found")

    share = await resolve_share_token(session, share_token, room.id)
    group_role = None
    if room.owner_group_id is not None:
        group_role = await fetch_group_role(session, user_id, room.owner_group_id)
    if not can_read(user, room, share, group_role=group_role):
        raise RoomNotFound.exception(f"Room {composed} not found")

    email = user.email

    # Leave previous room if any.
    old_room_id: str | None = sio_session.get("current_room_id")
    old_room_address: str | None = sio_session.get("current_room_address")
    if old_room_id is not None:
        await tsio.leave_room(sid, room_channel(old_room_id))
        if not old_room_id.startswith("@"):
            await _cleanup_session(redis, sid, sio_session, old_room_id)
        await tsio.emit(
            SessionLeft(
                room_id=old_room_address or old_room_id,
                user_id=user_id,
                sid=sid,
            ),
            room=room_channel(old_room_id),
        )

    # Join new room.
    await tsio.enter_room(sid, room_channel(room.id))
    sio_session["current_room_id"] = room.id
    sio_session["current_room_address"] = room.public_address
    sio_session["client_type"] = data.client_type
    await tsio.save_session(sid, sio_session)

    if data.client_type == "frontend":
        await tsio.enter_room(sid, "frontend")
        await tsio.enter_room(sid, room_channel("@global"))

    await tsio.emit(
        SessionJoined(
            room_id=room.public_address, user_id=user_id, sid=sid, email=email
        ),
        room=room_channel(room.id),
        skip_sid=sid,
    )

    room_step = room.step

    camera_key: str | None = None
    if data.client_type == "frontend":
        camera = Camera(owner=str(user_id))
        if room.default_camera:
            default_row = await session.get(
                RoomGeometry, (room.id, room.default_camera)
            )
            if default_row and default_row.type == "Camera":
                default_data = json.loads(default_row.config)
                updates = {
                    field: default_data[field]
                    for field in (
                        "position",
                        "target",
                        "up",
                        "fov",
                        "near",
                        "far",
                        "zoom",
                        "camera_type",
                    )
                    if field in default_data
                }
                if updates:
                    camera = camera.model_copy(update=updates)

        camera_key = f"cam:{email}:{sid[:8]}"
        camera_value = json.dumps(
            {"sid": sid, "email": email, "data": camera.model_dump()}
        )
        hash_key = RedisKey.room_cameras(room.id)
        await redis.hset(hash_key, camera_key, camera_value)  # type: ignore[misc]

        sio_session["camera_key"] = camera_key
        await tsio.save_session(sid, sio_session)

        await redis.hset(RedisKey.active_cameras(room.id), sid, camera_key)  # type: ignore[misc]

        await tsio.emit(
            GeometryInvalidate(room_id=room.id, operation="set", key=camera_key),
            room=room_channel(room.id),
        )

    frame_count = await storage.get_length(room.id)

    progress_raw = await redis.hgetall(RedisKey.room_progress(room.id))  # type: ignore[misc]
    progress_trackers = {
        pid: ProgressResponse(**json.loads(pdata))
        for pid, pdata in progress_raw.items()
    }

    return RoomJoinResponse(
        room_id=room.public_address,
        session_id=sid,
        step=room_step,
        frame_count=frame_count,
        camera_key=camera_key,
        default_camera=room.default_camera,
        progress_trackers=progress_trackers,
    )


@tsio.on(RoomLeave, emits=[SessionLeft, GeometryInvalidate])
async def room_leave(
    sid: str, data: RoomLeave, redis: RedisDep, session: SessionDep
) -> RoomLeaveResponse:
    """Leave a Socket.IO room."""
    from zndraw.dependencies import _load_room_by_address

    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")

    if current_room_id is None:
        raise NotInRoom.exception("Not currently in a room")

    room = await _load_room_by_address(session, data.owner_id, data.room_name)
    if room is None or current_room_id != room.id:
        return RoomLeaveResponse(room_id=f"{data.owner_id}/{data.room_name}")

    await tsio.leave_room(sid, room_channel(room.id))
    await _cleanup_session(redis, sid, sio_session, room.id)

    sio_session["current_room_id"] = None
    await tsio.save_session(sid, sio_session)

    await tsio.emit(
        SessionLeft(room_id=room.public_address, user_id=user_id, sid=sid),
        room=room_channel(room.id),
    )
    return RoomLeaveResponse(room_id=room.public_address)


async def _handle_typing(
    sid: str,
    owner_id: UUID,
    room_name: str,
    session: SessionDep,
    *,
    is_typing: bool,
) -> TypingResponse:
    """Broadcast typing status change to room."""
    from zndraw.dependencies import _load_room_by_address

    sio_session = await tsio.get_session(sid)
    user_id: UUID = sio_session["user_id"]
    current_room_id: str | None = sio_session.get("current_room_id")

    room = await _load_room_by_address(session, owner_id, room_name)
    if room is None or current_room_id != room.id:
        raise NotInRoom.exception("Not in this room")

    user = await session.get(User, user_id)
    email = user.email if user else None

    await tsio.emit(
        Typing(room_id=room.id, user_id=user_id, email=email, is_typing=is_typing),
        room=room_channel(room.id),
        skip_sid=sid,
    )
    return TypingResponse()


@tsio.on(TypingStart, emits=[Typing])
async def typing_start(
    sid: str, data: TypingStart, session: SessionDep
) -> TypingResponse:
    return await _handle_typing(
        sid, data.owner_id, data.room_name, session, is_typing=True
    )


@tsio.on(TypingStop, emits=[Typing])
async def typing_stop(
    sid: str, data: TypingStop, session: SessionDep
) -> TypingResponse:
    return await _handle_typing(
        sid, data.owner_id, data.room_name, session, is_typing=False
    )
