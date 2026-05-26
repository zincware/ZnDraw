"""FastAPI dependencies for database sessions, Redis, and authentication.

All resources are accessed from request.app.state.
Authentication uses zndraw-auth package.
"""

import json
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path as FilePath
from typing import Annotated, NamedTuple
from uuid import UUID

from fastapi import Depends, Header, Path, Request
from fastapi_users.authentication import JWTStrategy
from redis.asyncio import Redis as AsyncRedis
from sqlalchemy import or_
from sqlalchemy.ext.asyncio import async_sessionmaker
from sqlmodel import col, select
from sqlmodel.ext.asyncio.session import AsyncSession
from zndraw_socketio import AsyncServerWrapper

from zndraw.access import (
    GroupRole,
    ShareContext,
    can_edit,
    can_manage,
    can_read,
)
from zndraw.broadcast import room_channel  # re-export for legacy importers (Task 6 removes)
from zndraw.exceptions import (
    Forbidden,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    RoomReadOnly,
    SessionNotFound,
)
from zndraw.geometries import geometries as geometry_models
from zndraw.geometries.camera import Camera
from zndraw.models import Group, GroupMembership, Room, RoomGeometry, RoomShareLink
from zndraw.redis import RedisKey
from zndraw.storage import FrameStorage
from zndraw_auth import (
    SessionDep,
    User,
    current_active_user,
    current_superuser,
    current_user_scoped_session,
    get_auth_settings,
)
from zndraw_auth.db import get_session_maker
from zndraw_auth.settings import AuthSettings
from zndraw_joblib.dependencies import ResultBackend, validate_room_id
from zndraw_joblib.settings import JobLibSettings

# Re-export auth dependencies for convenience
CurrentUserDep = Annotated[User, Depends(current_active_user)]
AdminUserDep = Annotated[User, Depends(current_superuser)]


async def get_local_token_or_admin(
    request: Request,
    session: SessionDep,
    auth_settings: Annotated[AuthSettings, Depends(get_auth_settings)],
) -> User:
    """Admin-or-local-token auth. Tries local token first, then JWT superuser.

    - Local token matches ``app.state.local_token`` → synthetic superuser.
    - Otherwise, extract JWT from Authorization header and require superuser.
    """
    auth_header = request.headers.get("Authorization", "")

    # Path 1: Local admin token (fastest — no DB hit)
    local_token: str | None = getattr(request.app.state, "local_token", None)
    if local_token is not None and auth_header == f"Bearer {local_token}":
        return User(
            email="local-admin@localhost", hashed_password="", is_superuser=True
        )

    # Path 2: JWT → must resolve to an active superuser
    if auth_header.startswith("Bearer "):
        token = auth_header.removeprefix("Bearer ")
        strategy = JWTStrategy(
            secret=auth_settings.secret_key.get_secret_value(),
            lifetime_seconds=auth_settings.token_lifetime_seconds,
        )
        from fastapi_users.jwt import decode_jwt
        from jwt import InvalidTokenError

        user: User | None = None
        try:
            data = decode_jwt(
                token,
                secret=strategy.decode_key,
                audience=strategy.token_audience,
                algorithms=[strategy.algorithm],
            )
            user_id_raw = data.get("sub")
            if user_id_raw:
                user = await session.get(User, UUID(user_id_raw))
        except (InvalidTokenError, ValueError):
            user = None

        if user is not None and user.is_active:
            if not user.is_superuser:
                raise Forbidden.exception("Not a superuser")
            return user

    raise NotAuthenticated.exception("Not authenticated")


LocalTokenOrAdminDep = Annotated[User, Depends(get_local_token_or_admin)]

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


def get_frame_storage(request: Request) -> FrameStorage:
    """Get frame storage registry from app.state."""
    return request.app.state.frame_storage


FrameStorageDep = Annotated[FrameStorage, Depends(get_frame_storage)]


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
    """Verify room exists and return it, or raise RoomNotFound.

    Accepts either the surrogate UUID primary key or a composed
    ``<owner_uuid>/<room_name>`` address.
    """
    if "/" in room_id:
        owner_str, _, name_part = room_id.partition("/")
        try:
            owner_uuid = UUID(owner_str)
        except ValueError as exc:
            raise RoomNotFound.exception(f"Room with id {room_id} not found") from exc
        room = await _load_room_by_address(session, owner_uuid, name_part)
    else:
        room = await session.get(Room, room_id)
    if room is None:
        raise RoomNotFound.exception(f"Room with id {room_id} not found")
    return room


# =============================================================================
# Group membership helpers
# =============================================================================


async def fetch_my_group_ids(session: AsyncSession, user_id: UUID) -> list[UUID]:
    """Return all group ids the user is a member of (any role).

    Parameters
    ----------
    session
        Async database session.
    user_id
        The user whose memberships to query.
    """
    result = await session.exec(
        select(GroupMembership.group_id).where(GroupMembership.user_id == user_id)
    )
    return list(result.all())


async def fetch_group_role(
    session: AsyncSession, user_id: UUID, group_id: UUID
) -> GroupRole | None:
    """Return the user's role in the given group, or None if not a member."""
    result = await session.exec(
        select(GroupMembership).where(
            GroupMembership.user_id == user_id,
            GroupMembership.group_id == group_id,
        )
    )
    membership = result.first()
    return membership.role if membership is not None else None


class OwnerKind(StrEnum):
    USER = "user"
    GROUP = "group"


async def resolve_owner(
    session: AsyncSession, owner_id: UUID
) -> tuple[OwnerKind, str] | None:
    """Look up ``owner_id`` as a user (returns email) or group (returns name)."""
    user = await session.get(User, owner_id)
    if user is not None:
        return OwnerKind.USER, user.email
    group = await session.get(Group, owner_id)
    if group is not None:
        return OwnerKind.GROUP, group.name
    return None


async def get_my_group_ids(
    session: SessionDep, current_user: CurrentUserDep
) -> list[UUID]:
    """FastAPI dependency returning group ids for the current user."""
    return await fetch_my_group_ids(session, current_user.id)


MyGroupIdsDep = Annotated[list[UUID], Depends(get_my_group_ids)]


# =============================================================================
# Share-token resolver
# =============================================================================


async def resolve_share_token(
    session: AsyncSession, token: str | None, room_id: str
) -> ShareContext | None:
    """Resolve an ``X-Room-Share-Token`` header to a validated ShareContext.

    Returns ``None`` — never raises — when the token is missing, unknown,
    targets a different room, revoked, or expired.

    Parameters
    ----------
    session
        Async database session.
    token
        The raw token string from the request header, or ``None`` when the
        header is absent.
    room_id
        The room the token must be scoped to.
    """
    if token is None:
        return None
    result = await session.exec(
        select(RoomShareLink).where(RoomShareLink.token == token)
    )
    link = result.first()
    if link is None or link.room_id != room_id:
        return None
    if link.revoked_at is not None:
        return None
    if link.expires_at is not None and link.expires_at <= datetime.now(UTC):
        return None
    return ShareContext(room_id=link.room_id, access=link.access)


async def get_share_context(
    session: SessionDep,
    room_id: str = Path(),
    x_room_share_token: str | None = Header(default=None, alias="X-Room-Share-Token"),
) -> ShareContext | None:
    """FastAPI dependency resolving the share-token header to a ShareContext."""
    return await resolve_share_token(session, x_room_share_token, room_id)


ShareTokenDep = Annotated[ShareContext | None, Depends(get_share_context)]


# =============================================================================
# Edit lock helper
# =============================================================================


async def _check_edit_lock(
    redis: AsyncRedis,  # type: ignore[type-arg]
    room_id: str,
    lock_token: str | None = None,
) -> None:
    """Redis edit-lock check only — no admin lock, no permission logic."""
    raw = await redis.get(RedisKey.edit_lock(room_id))
    if raw is None:
        return
    holder = json.loads(raw)
    if lock_token is None:
        raise RoomLocked.exception("Room is being edited; Lock-Token required")
    if holder["lock_token"] != lock_token:
        raise RoomLocked.exception("Room is being edited by another session")


# =============================================================================
# Session camera helpers
# =============================================================================


async def get_verified_session_id(
    current_user: CurrentUserDep,
    redis: RedisDep,
    session: SessionDep,
    owner_id: UUID = Path(),
    room_name: str = Path(),
    session_id: str = Path(),
) -> str:
    """Verify session belongs to the current user, or raise 404.

    Finds the session's own camera in room_cameras (by SID match) to
    verify ownership. The active camera may point to another user's
    camera (sessions can view through any camera in the room), so we
    cannot use the active_cameras chain for ownership verification.
    """
    room = await _load_room_by_address(session, owner_id, room_name)
    if room is None:
        raise RoomNotFound.exception(f"Room {owner_id}/{room_name} not found")
    room_id = room.id
    if not await redis.hexists(RedisKey.active_cameras(room_id), session_id):  # type: ignore[misc]
        raise SessionNotFound.exception("Session not found")
    all_cameras: dict[str, str] = await redis.hgetall(  # type: ignore[misc]
        RedisKey.room_cameras(room_id)
    )
    uid = str(current_user.id)
    for raw in all_cameras.values():
        entry = json.loads(raw)
        if entry.get("sid") == session_id and Camera(**entry["data"]).owner == uid:
            return session_id
    raise SessionNotFound.exception("Session not found")


VerifiedSessionDep = Annotated[str, Depends(get_verified_session_id)]


async def get_active_session_cam_id(
    redis: RedisDep,
    session: SessionDep,
    owner_id: UUID = Path(),
    room_name: str = Path(),
    session_id: str = Path(),
) -> str:
    """Verify session exists in active-cameras (no ownership check).

    Use for read-only access where any room participant may view
    session state. For mutations, use ``VerifiedSessionDep``.
    """
    room = await _load_room_by_address(session, owner_id, room_name)
    if room is None:
        raise RoomNotFound.exception(f"Room {owner_id}/{room_name} not found")
    if not await redis.hexists(RedisKey.active_cameras(room.id), session_id):  # type: ignore[misc]
        raise SessionNotFound.exception("Session not found")
    return session_id


ActiveSessionCamDep = Annotated[str, Depends(get_active_session_cam_id)]


# =============================================================================
# Geometry write access
# =============================================================================


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
    share: ShareContext | None,
    lock_token: str | None = None,
) -> WritableGeometryInfo:
    """Verify: (1) Redis edit lock, (2) can_edit permission, (3) geometry ownership."""
    validate_room_id(room_id)
    room = await verify_room(session, room_id)

    # 1. Redis edit-lock (serialization, not auth)
    raw = await redis.get(RedisKey.edit_lock(room_id))
    if raw is not None:
        holder = json.loads(raw)
        if lock_token is not None:
            if holder["lock_token"] != lock_token:
                raise RoomLocked.exception("Room is being edited by another session")
        elif holder["user_id"] != str(current_user.id):
            raise RoomLocked.exception("Room is being edited by another user")

    # 2. can_edit gate
    group_role: GroupRole | None = None
    if room.owner_group_id is not None:
        group_role = await fetch_group_role(
            session, current_user.id, room.owner_group_id
        )
    if not can_edit(current_user, room, share, group_role=group_role):
        raise Forbidden.exception("You may not edit this room")

    # 3. Per-geometry ownership
    current_owner = await get_owner_from_geometry(redis, session, room_id, geometry_key)
    user_id_str = str(current_user.id)
    is_unowned = current_owner is None
    is_mine = current_owner == user_id_str
    if not current_user.is_superuser and not is_unowned and not is_mine:
        raise Forbidden.exception("Not the geometry owner")

    return WritableGeometryInfo(room=room, current_owner=current_owner)


async def get_writable_room_id(
    request: Request,
    session: SessionDep,
    current_user: CurrentUserDep,
    redis: RedisDep,
    room_id: str = Path(),
    x_room_share_token: str | None = Header(default=None, alias="X-Room-Share-Token"),
) -> str:
    """Verify a room is writable and return the surrogate UUID string."""
    validate_room_id(room_id)
    if room_id in ("@global", "@internal"):
        return room_id
    owner_part, _, name_part = room_id.partition("/")
    room = await _load_room_by_address(session, UUID(owner_part), name_part)
    if room is None:
        raise RoomNotFound.exception(f"Room {room_id} not found")
    share = await resolve_share_token(session, x_room_share_token, room.id)
    group_role: GroupRole | None = None
    if room.owner_group_id is not None:
        group_role = await fetch_group_role(
            session, current_user.id, room.owner_group_id
        )
    if not can_edit(current_user, room, share, group_role=group_role):
        raise Forbidden.exception("You may not edit this room")
    lock_token = request.headers.get("Lock-Token")
    await _check_edit_lock(redis, room.id, lock_token)
    return room.id


# =============================================================================
# Access composites — room + auth context bundled together
# =============================================================================


class AccessContext(NamedTuple):
    """Bundle of room + resolved auth context, reused across access deps."""

    room: Room
    share: ShareContext | None
    group_role: GroupRole | None


async def _load_room_by_address(
    session: AsyncSession, owner_id: UUID, room_name: str
) -> Room | None:
    """Look up ``(owner_id, room_name)`` via the unique index."""
    result = await session.exec(
        select(Room).where(
            or_(
                col(Room.owner_user_id) == owner_id,
                col(Room.owner_group_id) == owner_id,
            ),
            col(Room.room_name) == room_name,
        )
    )
    return result.first()


async def _load_access_context(
    session: AsyncSession,
    owner_id: UUID,
    room_name: str,
    current_user: User,
    share: ShareContext | None,
) -> AccessContext:
    """Load room from DB and resolve group role for the current user."""
    room = await _load_room_by_address(session, owner_id, room_name)
    if room is None:
        raise RoomNotFound.exception(f"Room {owner_id}/{room_name} not found")
    group_role: GroupRole | None = None
    if room.owner_group_id is not None:
        group_role = await fetch_group_role(
            session, current_user.id, room.owner_group_id
        )
    return AccessContext(room=room, share=share, group_role=group_role)


async def get_share_context_two_segment(
    session: SessionDep,
    owner_id: UUID = Path(),
    room_name: str = Path(),
    x_room_share_token: str | None = Header(default=None, alias="X-Room-Share-Token"),
) -> ShareContext | None:
    """Resolve the share-token header for the two-segment path."""
    if x_room_share_token is None:
        return None
    room = await _load_room_by_address(session, owner_id, room_name)
    if room is None:
        return None
    return await resolve_share_token(session, x_room_share_token, room.id)


TwoSegmentShareTokenDep = Annotated[
    ShareContext | None, Depends(get_share_context_two_segment)
]


async def get_readable_room(
    session: SessionDep,
    current_user: CurrentUserDep,
    share: TwoSegmentShareTokenDep,
    owner_id: UUID = Path(),
    room_name: str = Path(),
) -> AccessContext:
    """Load room + auth context; raise 404 if caller cannot read."""
    ctx = await _load_access_context(session, owner_id, room_name, current_user, share)
    if not can_read(current_user, ctx.room, ctx.share, group_role=ctx.group_role):
        raise RoomNotFound.exception(f"Room {owner_id}/{room_name} not found")
    return ctx


async def get_editable_room(
    ctx: Annotated[AccessContext, Depends(get_readable_room)],
    current_user: CurrentUserDep,
) -> AccessContext:
    if not can_edit(current_user, ctx.room, ctx.share, group_role=ctx.group_role):
        raise Forbidden.exception("You may not edit this room")
    return ctx


async def get_manageable_room(
    ctx: Annotated[AccessContext, Depends(get_readable_room)],
    current_user: CurrentUserDep,
) -> AccessContext:
    if not can_manage(current_user, ctx.room, group_role=ctx.group_role):
        raise Forbidden.exception("You may not manage this room")
    return ctx


AccessReadDep = Annotated[AccessContext, Depends(get_readable_room)]
AccessEditDep = Annotated[AccessContext, Depends(get_editable_room)]
AccessManageDep = Annotated[AccessContext, Depends(get_manageable_room)]


async def get_writable_geometry(
    request: Request,
    access: AccessEditDep,
    session: SessionDep,
    current_user: CurrentUserDep,
    redis: RedisDep,
    key: str = Path(),
) -> WritableGeometryInfo:
    """Edit-gated room + Redis edit-lock + per-geometry ownership check."""
    lock_token = request.headers.get("Lock-Token")
    await _check_edit_lock(redis, access.room.id, lock_token)
    current_owner = await get_owner_from_geometry(redis, session, access.room.id, key)
    user_id_str = str(current_user.id)
    if (
        not current_user.is_superuser
        and current_owner is not None
        and current_owner != user_id_str
    ):
        raise Forbidden.exception("Not the geometry owner")
    return WritableGeometryInfo(room=access.room, current_owner=current_owner)


WritableGeometryDep = Annotated[WritableGeometryInfo, Depends(get_writable_geometry)]


async def require_writable_room(
    storage: FrameStorageDep,
    access: AccessEditDep,
) -> None:
    """Raise RoomReadOnly if the room has a provider mount."""
    if await storage.has_mount(access.room.id):
        raise RoomReadOnly.exception("Room is provider-backed (read-only)")


RequireWritableDep = Annotated[None, Depends(require_writable_room)]


async def get_writable_room(
    request: Request,
    access: AccessEditDep,
    redis: RedisDep,
) -> Room:
    """Edit-gated room + Redis edit-lock coordination."""
    lock_token = request.headers.get("Lock-Token")
    await _check_edit_lock(redis, access.room.id, lock_token)
    return access.room


WritableRoomDep = Annotated[Room, Depends(get_writable_room)]
