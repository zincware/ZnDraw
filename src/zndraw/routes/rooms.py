"""Room REST API endpoints - simplified for frontend compatibility.

Handles room creation, listing, and basic metadata operations.
Uses string UUIDs for room IDs to match frontend expectations.
"""

import json
import logging
import re
from typing import Annotated, Any

logger = logging.getLogger(__name__)

from fastapi import APIRouter, Query, status
from sqlalchemy.ext.asyncio import AsyncSession
from sqlmodel import select

from zndraw.dependencies import (
    CurrentUserDep,
    OptionalUserDep,
    RedisDep,
    SessionDep,
    SioDep,
    StorageDep,
    StorageRouterDep,
    WritableRoomDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    Forbidden,
    InvalidPayload,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    RoomReadOnly,
    UnprocessableContent,
    problem_responses,
)
from zndraw.geometries import geometries as geometry_models
from zndraw.materials import MeshBasicMaterial
from zndraw.transformations import InArrayTransform
from zndraw.models import (
    Room,
    RoomBookmark,
    RoomFigure,
    RoomGeometry,
    SelectionGroup,
    ServerSettings,
)
from zndraw.redis import RedisKey
from zndraw.schemas import (
    CollectionResponse,
    PresenceResponse,
    PresenceSessionResponse,
    RoomCreate,
    RoomCreateResponse,
    RoomPatchRequest,
    RoomPatchResponse,
    RoomResponse,
    SessionsListResponse,
)
from zndraw.socket_events import FramesInvalidate, RoomUpdate

router = APIRouter(prefix="/v1/rooms", tags=["rooms"])


def _initialize_default_geometries(session: AsyncSession, room_id: str) -> None:
    """Initialize default geometries for a new room.

    Creates the standard geometry set: particles, bonds, curve, cell, floor.
    These geometries use property references (e.g., "arrays.positions") that
    are resolved at render time from frame data.
    Rows are added to the session but NOT committed — caller must commit.
    """
    defaults: dict[str, tuple[str, dict[str, Any]]] = {
        "particles": (
            "Sphere",
            {
                "active": True,
                "position": "arrays.positions",
                "color": "arrays.colors",
                "radius": "arrays.radii",
                "scale": [[0.7, 0.7, 0.7]],
            },
        ),
        "bonds": (
            "Bond",
            {
                "active": True,
                "position": "arrays.positions",
                "color": "arrays.colors",
                "scale": 0.15,
            },
        ),
        "curve": ("Curve", {"active": True}),
        "cell": ("Cell", {"active": True}),
        "floor": ("Floor", {"active": False}),
        "constraints-fixed-atoms": (
            "Sphere",
            {
                "active": True,
                "position": InArrayTransform(
                    source="constraints",
                    path="0.kwargs.indices",
                    filter="arrays.positions",
                ),
                "radius": InArrayTransform(
                    source="constraints",
                    path="0.kwargs.indices",
                    filter="arrays.radii",
                ),
                "color": ["#FF0000"],
                "material": MeshBasicMaterial(wireframe=True),
                "scale": [(0.71, 0.71, 0.71)],
                "selecting": {"enabled": False},
                "hovering": {"enabled": False},
            },
        ),
    }

    for key, (type_name, data) in defaults.items():
        model_cls = geometry_models[type_name]
        config_json = model_cls(**data).model_dump_json()
        session.add(
            RoomGeometry(
                room_id=room_id,
                key=key,
                type=type_name,
                config=config_json,
            )
        )


async def _copy_room_state(
    session: AsyncSession, source_room_id: str, target_room_id: str
) -> None:
    """Copy geometries, bookmarks, figures, and selection groups.

    Copies all room state except frames (handled separately) and
    owned geometries (ephemeral session cameras).
    Rows are added to the session but NOT committed — caller must commit.
    """
    # Copy geometries (skip owned — ephemeral session cameras)
    result = await session.execute(
        select(RoomGeometry).where(RoomGeometry.room_id == source_room_id)
    )
    for row in result.scalars().all():
        config = json.loads(row.config)
        if config.get("owner") is not None:
            continue
        session.add(
            RoomGeometry(
                room_id=target_room_id,
                key=row.key,
                type=row.type,
                config=row.config,
                selection=row.selection,
            )
        )

    # Copy bookmarks
    result = await session.execute(
        select(RoomBookmark).where(RoomBookmark.room_id == source_room_id)
    )
    for row in result.scalars().all():
        session.add(
            RoomBookmark(
                room_id=target_room_id,
                frame_index=row.frame_index,
                label=row.label,
            )
        )

    # Copy figures
    result = await session.execute(
        select(RoomFigure).where(RoomFigure.room_id == source_room_id)
    )
    for row in result.scalars().all():
        session.add(
            RoomFigure(
                room_id=target_room_id,
                key=row.key,
                type=row.type,
                data=row.data,
            )
        )

    # Copy selection groups
    result = await session.execute(
        select(SelectionGroup).where(SelectionGroup.room_id == source_room_id)
    )
    for row in result.scalars().all():
        session.add(
            SelectionGroup(
                room_id=target_room_id,
                name=row.name,
                selections=row.selections,
            )
        )


async def _get_default_room_id(session: AsyncSession) -> str | None:
    """Get the default room ID from ServerSettings."""
    settings = await session.get(ServerSettings, 1)
    return settings.default_room_id if settings else None


async def build_room_update(
    session: AsyncSession,
    storage: Any,
    room: Room,
) -> RoomUpdate:
    """Build a full RoomUpdate snapshot from DB + storage."""
    default_room_id = await _get_default_room_id(session)
    frame_count = await storage.get_length(room.id)
    return RoomUpdate(
        id=room.id,
        description=room.description,
        frame_count=frame_count,
        locked=room.locked,
        is_default=(room.id == default_room_id),
    )


# =============================================================================
# Room CRUD
# =============================================================================


@router.post(
    "",
    response_model=RoomCreateResponse,
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomReadOnly),
)
async def create_room(
    session: SessionDep,
    storage: StorageDep,
    storage_router: StorageRouterDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    request: RoomCreate,
) -> RoomCreateResponse:
    """Create a new room.

    The room ID must contain only alphanumeric characters, hyphens, and
    underscores.  System rooms (prefixed with ``@``) cannot be created via
    the REST API.

    The ``copyFrom`` field controls frame initialization:

    - ``"<room-id>"`` - deep copy from that room (frames, geometries,
      bookmarks, figures, selection groups, step)
    - ``"@empty"`` - one empty frame (bypasses server default)
    - ``"@none"`` - zero frames (bypasses server default)
    - *omitted* - use server default room, fallback to one empty frame
    """
    room_id = request.room_id

    # Validate room ID format: alphanumeric, hyphens, and underscores only
    if not re.match(r"^[a-zA-Z0-9\-_]+$", room_id):
        raise InvalidPayload.exception(
            "Room ID must contain only alphanumeric characters, "
            "hyphens, and underscores"
        )

    # Check if room already exists
    existing = await session.get(Room, room_id)
    if existing is not None:
        frame_count = await storage.get_length(room_id)
        return RoomCreateResponse(
            status="ok",
            room_id=room_id,
            frame_count=frame_count,
            created=False,
        )

    # Resolve copyFrom: @-prefixed presets, room IDs, or server default
    copy_from = request.copy_from
    if copy_from is None:
        # No explicit copyFrom — check server default
        default_room_id = await _get_default_room_id(session)
        copy_from = default_room_id if default_room_id else "@empty"

    # Reject unknown @-prefixed presets
    _PRESETS = {"@empty", "@none"}
    if copy_from.startswith("@") and copy_from not in _PRESETS:
        raise UnprocessableContent.exception(
            f"Unknown preset '{copy_from}'. Valid presets: {', '.join(sorted(_PRESETS))}"
        )

    # Try to resolve source room for room-id copyFrom
    source_room: Room | None = None
    if not copy_from.startswith("@"):
        source_room = await session.get(Room, copy_from)
        if source_room is not None and await storage_router.has_mount(copy_from):
            raise RoomReadOnly.exception(
                "Cannot copy from a room with a mounted source"
            )

    # Create room in database
    room = Room(
        id=room_id,
        description=request.description,
        created_by_id=current_user.id,
        step=source_room.step if source_room else 0,
    )
    session.add(room)

    frame_count = 0
    if copy_from == "@none":
        pass  # zero frames
    elif copy_from == "@empty":
        await storage.extend(room_id, [{}])
        frame_count = 1
    elif source_room is not None:
        # Deep copy from existing room: frames + all state
        source_frames_or_none = await storage.get_range(copy_from, 0, None)
        source_frames = [f for f in source_frames_or_none if f is not None]
        if source_frames:
            await storage.extend(room_id, source_frames)
            frame_count = len(source_frames)
        await _copy_room_state(session, copy_from, room_id)
    else:
        # copyFrom refers to a non-existent room — fall back to empty
        await storage.extend(room_id, [{}])
        frame_count = 1

    # Only create default geometries when NOT copying from an existing room
    if source_room is None:
        _initialize_default_geometries(session, room_id)

    await session.commit()

    # Broadcast room creation to @overview system room
    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")

    return RoomCreateResponse(
        status="ok",
        room_id=room_id,
        frame_count=event.frame_count,
        created=True,
    )


@router.get("")
async def list_rooms(
    session: SessionDep,
    storage: StorageDep,
    _current_user: OptionalUserDep,
    search: Annotated[str | None, Query(description="Search pattern")] = None,
) -> CollectionResponse[RoomResponse]:
    """List available rooms.

    Returns a flat array of room objects matching frontend expectations.
    Authentication is optional - unauthenticated users see public rooms only.
    """
    # Get the default room ID for isDefault computation
    default_room_id = await _get_default_room_id(session)

    # Get all rooms (for PoC, show all public rooms)
    statement = select(Room).where(Room.is_public == True)  # noqa: E712
    result = await session.execute(statement)
    rooms = list(result.scalars().all())

    # Build response with frame counts
    room_responses = []
    for room in rooms:
        frame_count = await storage.get_length(room.id)

        # Apply search filter if provided
        if search:
            search_lower = search.lower()
            if search_lower not in room.id.lower() and (
                room.description is None or search_lower not in room.description.lower()
            ):
                continue

        room_responses.append(
            RoomResponse(
                id=room.id,
                description=room.description,
                frame_count=frame_count,
                locked=room.locked,
                is_default=(room.id == default_room_id),
            )
        )

    return CollectionResponse(items=room_responses)


@router.get(
    "/{room_id}",
    response_model=RoomResponse,
    responses=problem_responses(RoomNotFound),
)
async def get_room(
    session: SessionDep,
    storage: StorageDep,
    _current_user: CurrentUserDep,
    room_id: str,
) -> RoomResponse:
    """Get details of a specific room."""
    room = await verify_room(session, room_id)
    frame_count = await storage.get_length(room_id)
    default_room_id = await _get_default_room_id(session)

    return RoomResponse(
        id=room.id,
        description=room.description,
        frame_count=frame_count,
        locked=room.locked,
        is_default=(room.id == default_room_id),
    )


@router.get(
    "/{room_id}/presence",
    response_model=PresenceResponse,
    responses=problem_responses(RoomNotFound),
)
async def get_room_presence(
    session: SessionDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
) -> PresenceResponse:
    """Get presence (online users) for a room.

    Returns all currently connected sessions in the room.
    Each session (tab/connection) is listed individually.
    """
    # Verify room exists
    await verify_room(session, room_id)

    # Scan for all presence keys for this room
    pattern = RedisKey.presence_sid_pattern(room_id)
    sessions: list[PresenceSessionResponse] = []

    async for key in redis.scan_iter(match=pattern):  # type: ignore[misc]
        # Parse SID from key
        sid = RedisKey.parse_presence_sid(key)
        if sid is None:
            continue

        # Get user_id stored in the key
        user_id_raw = await redis.get(key)  # type: ignore[misc]
        if user_id_raw is None:
            continue

        from uuid import UUID as _UUID

        user_id = _UUID(user_id_raw)

        # Get user info from database
        from zndraw_auth import User

        user = await session.get(User, user_id)
        if user is None:
            continue

        sessions.append(
            PresenceSessionResponse(
                sid=sid,
                user_id=user_id,
                email=user.email,
            )
        )

    return PresenceResponse(items=sessions)


@router.get(
    "/{room_id}/sessions",
    response_model=SessionsListResponse,
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_user_sessions(
    session: SessionDep,
    redis: RedisDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> SessionsListResponse:
    """List the current user's active frontend sessions in this room.

    Frontend sessions are identified by having an entry in the
    active-cameras hash (pyclients don't get cameras).
    """
    await verify_room(session, room_id)

    # Get all frontend SIDs (those with active cameras)
    all_active: dict[str, str] = await redis.hgetall(  # type: ignore[misc]
        RedisKey.active_cameras(room_id)
    )

    # Filter to current user's sessions via presence keys (single round-trip)
    sids = list(all_active.keys())
    if sids:
        presence_keys = [RedisKey.presence_sid(room_id, sid) for sid in sids]
        user_ids: list[str | None] = await redis.mget(presence_keys)  # type: ignore[misc]
        uid = str(current_user.id)
        user_sids = [sid for sid, owner in zip(sids, user_ids) if owner == uid]
    else:
        user_sids = []

    return SessionsListResponse(items=user_sids)


@router.patch(
    "/{room_id}",
    responses=problem_responses(RoomNotFound, RoomLocked, Forbidden),
)
async def update_room(
    session: SessionDep,
    storage: StorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    updates: RoomPatchRequest,
) -> RoomPatchResponse:
    """Update room metadata (description, locked, frame_count).

    Requires writable access (room must not be locked by another user).
    Setting ``frame_count`` stores an external frame count in the storage
    backend (only supported by ``StorageRouter`` with Redis).  Use 0 to
    clear the external count.
    """
    changed = False
    if updates.description is not None:
        room.description = updates.description
        changed = True

    if updates.locked is not None:
        room.locked = updates.locked
        changed = True

    if updates.frame_count is not None:
        count = updates.frame_count
        try:
            if count > 0:
                await storage.set_frame_count(room.id, count)
            else:
                await storage.clear_frame_count(room.id)
        except NotImplementedError:
            raise Forbidden.exception(
                "Storage backend does not support external frame counts"
            )
        await sio.emit(
            FramesInvalidate(room_id=room.id, action="clear", count=count),
            room=room_channel(room.id),
        )
        changed = True

    await session.commit()

    if changed:
        event = await build_room_update(session, storage, room)
        logger.debug("Broadcasting RoomUpdate: %s", event.model_dump())
        await sio.emit(event, room=f"room:{room.id}")
        await sio.emit(event, room="room:@overview")

    return RoomPatchResponse()
