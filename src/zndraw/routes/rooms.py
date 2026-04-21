"""Room REST API endpoints - simplified for frontend compatibility.

Handles room creation, listing, and basic metadata operations.
Uses string UUIDs for room IDs to match frontend expectations.
"""

import json
import re
from typing import Annotated, Any

from fastapi import APIRouter, Query, status
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.dependencies import (
    CurrentUserDep,
    FrameStorageDep,
    OptionalUserDep,
    RedisDep,
    SessionDep,
    SioDep,
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
from zndraw.geometries.camera import Camera
from zndraw.geometries.fog import Fog
from zndraw.geometries.lights import (
    AmbientLight,
    DirectionalLight,
    HemisphereLight,
    LightPosition,
)
from zndraw.geometries.pathtracing import PathTracing
from zndraw.geometries.property_inspector import PropertyInspector
from zndraw.materials import MeshBasicMaterial
from zndraw.models import (
    Room,
    RoomBookmark,
    RoomFigure,
    RoomGeometry,
    RoomMembership,
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
    SessionItem,
    SessionsListResponse,
)
from zndraw.socket_events import FramesInvalidate, RoomUpdate
from zndraw.storage import FrameStorage
from zndraw.transformations import InArrayTransform

router = APIRouter(prefix="/v1/rooms", tags=["rooms"])


def _initialize_default_geometries(session: AsyncSession, room_id: str) -> None:
    """Initialize default geometries for a new room.

    Creates the standard geometry set: particles, bonds, curve, cell, floor.
    These geometries use property references (e.g., "arrays.positions") that
    are resolved at render time from frame data.

    Also creates default scene objects for lighting, fog, pathtracing, and
    property inspector.

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

    # Scene objects: lights, fog, pathtracing, property inspector
    scene_objects: dict[str, tuple[str, Any]] = {
        "key-light": (
            "DirectionalLight",
            DirectionalLight(
                position=LightPosition(x=5.0, y=2.0, z=8.0),
                intensity=0.7,
            ),
        ),
        "fill-light": (
            "DirectionalLight",
            DirectionalLight(
                position=LightPosition(x=-4.0, y=-1.0, z=6.0),
                intensity=0.4,
                color="#a0c4ff",
            ),
        ),
        "rim-light": (
            "DirectionalLight",
            DirectionalLight(
                position=LightPosition(x=0.0, y=0.0, z=-50.0),
                intensity=0.5,
                color="#fff0f5",
            ),
        ),
        "ambient-light": ("AmbientLight", AmbientLight(intensity=0.35)),
        "hemisphere-light": ("HemisphereLight", HemisphereLight(intensity=0.3)),
        "fog": ("Fog", Fog(active=True, near=180.0, far=300.0)),
        "pathtracing": ("PathTracing", PathTracing(active=False)),
        "property-inspector": ("PropertyInspector", PropertyInspector(active=False)),
    }

    for key, (type_name, model) in scene_objects.items():
        session.add(
            RoomGeometry(
                room_id=room_id,
                key=key,
                type=type_name,
                config=model.model_dump_json(),
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
    result = await session.exec(
        select(RoomGeometry).where(RoomGeometry.room_id == source_room_id)
    )
    for row in result.all():
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
    result = await session.exec(
        select(RoomBookmark).where(RoomBookmark.room_id == source_room_id)
    )
    for row in result.all():
        session.add(
            RoomBookmark(
                room_id=target_room_id,
                frame_index=row.frame_index,
                label=row.label,
            )
        )

    # Copy figures
    result = await session.exec(
        select(RoomFigure).where(RoomFigure.room_id == source_room_id)
    )
    for row in result.all():
        session.add(
            RoomFigure(
                room_id=target_room_id,
                key=row.key,
                type=row.type,
                data=row.data,
            )
        )

    # Copy selection groups
    result = await session.exec(
        select(SelectionGroup).where(SelectionGroup.room_id == source_room_id)
    )
    for row in result.all():
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
    storage: FrameStorage,
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


async def broadcast_room_update(
    sio,
    session: AsyncSession,
    storage: FrameStorage,
    room: Room,
) -> None:
    """Broadcast a full RoomUpdate to every authorized viewer.

    Public rooms go to the shared ``rooms:feed`` channel — every
    authenticated socket is a member, including any client currently
    joined to ``room:{id}``, so this single emit reaches in-room
    viewers too.

    Private rooms fan out to each member's ``user:{uid}`` channel,
    which similarly covers both in-room and out-of-room members.
    """
    event = await build_room_update(session, storage, room)
    if room.is_public:
        await sio.emit(event, room="rooms:feed")
        return
    result = await session.exec(
        select(RoomMembership.user_id).where(
            RoomMembership.room_id == room.id
        )
    )
    for uid in result.all():
        await sio.emit(event, room=f"user:{uid}")


# =============================================================================
# Room CRUD
# =============================================================================


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomReadOnly),
)
async def create_room(
    session: SessionDep,
    storage: FrameStorageDep,
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
    presets = {"@empty", "@none"}
    if copy_from.startswith("@") and copy_from not in presets:
        raise UnprocessableContent.exception(
            f"Unknown preset '{copy_from}'. Valid presets: {', '.join(sorted(presets))}"
        )

    # Try to resolve source room for room-id copyFrom
    source_room: Room | None = None
    if not copy_from.startswith("@"):
        source_room = await session.get(Room, copy_from)
        if source_room is not None and await storage.has_mount(copy_from):
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
        await storage[room_id].extend([{}])
        frame_count = 1
    elif source_room is not None:
        # Deep copy from existing room: frames + all state
        source_frames_or_none = await storage[copy_from][0:].to_list()
        source_frames = [f for f in source_frames_or_none if f is not None]
        if source_frames:
            await storage[room_id].extend(source_frames)
            frame_count = len(source_frames)
        await _copy_room_state(session, copy_from, room_id)
    else:
        # copyFrom refers to a non-existent room — fall back to empty
        await storage[room_id].extend([{}])
        frame_count = 1

    # Only create default geometries when NOT copying from an existing room
    if source_room is None:
        _initialize_default_geometries(session, room_id)

    await session.commit()

    await broadcast_room_update(sio, session, storage, room)

    return RoomCreateResponse(
        status="ok",
        room_id=room_id,
        frame_count=frame_count,
        created=True,
    )


@router.get("")
async def list_rooms(
    session: SessionDep,
    storage: FrameStorageDep,
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
    statement = select(Room).where(Room.is_public.is_(True))
    result = await session.exec(statement)
    rooms = list(result.all())

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
    responses=problem_responses(RoomNotFound),
)
async def get_room(
    session: SessionDep,
    storage: FrameStorageDep,
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
    responses=problem_responses(RoomNotFound),
)
async def get_room_presence(
    session: SessionDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
) -> PresenceResponse:
    """Get presence (online users) for a room.

    Derives presence from the camera hash — each frontend session with an
    active camera is considered present. Pyclients do not have cameras and
    are not included.
    """
    from uuid import UUID as _UUID

    await verify_room(session, room_id)

    cameras_raw: dict[str, str] = await redis.hgetall(  # type: ignore[misc]
        RedisKey.room_cameras(room_id)
    )
    sessions_list: list[PresenceSessionResponse] = []
    for raw_value in cameras_raw.values():
        entry = json.loads(raw_value)
        camera = Camera(**entry["data"])
        if camera.owner is None:
            continue
        sessions_list.append(
            PresenceSessionResponse(
                sid=entry["sid"],
                user_id=_UUID(camera.owner),
                email=entry.get("email"),
            )
        )

    return PresenceResponse(items=sessions_list)


@router.get(
    "/{room_id}/sessions",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_sessions(
    session: SessionDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    email: Annotated[str | None, Query(description="Filter by user email")] = None,
) -> SessionsListResponse:
    """List all active frontend sessions in this room.

    Returns every frontend session (identified by having an entry in the
    active-cameras hash). Optional ``email`` query param filters by user.
    """
    await verify_room(session, room_id)

    all_active: dict[str, str] = await redis.hgetall(  # type: ignore[misc]
        RedisKey.active_cameras(room_id)
    )
    if not all_active:
        return SessionsListResponse(items=[])

    sids = list(all_active.keys())
    camera_keys = list(all_active.values())
    raw_cameras: list[str | None] = await redis.hmget(  # type: ignore[assignment]
        RedisKey.room_cameras(room_id), camera_keys
    )

    items: list[SessionItem] = []
    for sid, cam_key, raw in zip(sids, camera_keys, raw_cameras, strict=False):
        if raw is None:
            continue
        entry = json.loads(raw)
        entry_email = entry.get("email", "")
        if email is not None and entry_email != email:
            continue
        items.append(SessionItem(sid=sid, email=entry_email, camera_key=cam_key))

    return SessionsListResponse(items=items)


@router.patch(
    "/{room_id}",
    responses=problem_responses(RoomNotFound, RoomLocked, Forbidden),
)
async def update_room(
    session: SessionDep,
    storage: FrameStorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    updates: RoomPatchRequest,
    room_id: str,  # noqa: ARG001
) -> RoomPatchResponse:
    """Update room metadata (description, locked, frame_count).

    Requires writable access (room must not be locked by another user).
    Setting ``frame_count`` stores an external frame count in Redis
    (used by provider-backed rooms).  Use 0 to clear the external count.
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
        if count > 0:
            await storage.set_frame_count(room.id, count)
        else:
            await storage.clear_frame_count(room.id)
        await sio.emit(
            FramesInvalidate(room_id=room.id, action="clear", count=count),
            room=room_channel(room.id),
        )
        changed = True

    await session.commit()

    if changed:
        await broadcast_room_update(sio, session, storage, room)

    return RoomPatchResponse()
