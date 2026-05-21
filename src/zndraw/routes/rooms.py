"""Room REST API endpoints."""

import json
from typing import Annotated, Any
from uuid import UUID

from fastapi import APIRouter, Query, Response, status
from sqlalchemy import or_
from sqlalchemy.exc import IntegrityError
from sqlmodel import col, select
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import Visibility
from zndraw.config import SettingsDep
from zndraw.dependencies import (
    AccessManageDep,
    AccessReadDep,
    CurrentUserDep,
    FrameStorageDep,
    MyGroupIdsDep,
    OwnerKind,
    RedisDep,
    SessionDep,
    SioDep,
    _load_room_by_address,
    fetch_group_role,
    resolve_owner,
    room_channel,
)
from zndraw.exceptions import (
    Forbidden,
    InvalidPayload,
    NotAuthenticated,
    RoomNotFound,
    RoomReadOnly,
    ShareLinkInvalid,
    TransferTargetInvalid,
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
    owner_id = room.owner_user_id or room.owner_group_id
    assert owner_id is not None, "Room must have exactly one owner"
    resolved = await resolve_owner(session, owner_id)
    owner_kind: str = "user"
    owner_label: str = ""
    if resolved is not None:
        owner_kind = resolved[0].value
        owner_label = resolved[1]
    return RoomUpdate(
        room_id=room.public_address,
        description=room.description,
        frame_count=frame_count,
        visibility=room.visibility,
        owner_id=owner_id,
        owner_kind=owner_kind,  # type: ignore[arg-type]
        owner_label=owner_label,
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
    if room.visibility == Visibility.PUBLIC:
        await sio.emit(event, room="rooms:feed")
        return
    if room.owner_group_id is not None:
        from zndraw.models import GroupMembership

        result = await session.exec(
            select(GroupMembership.user_id).where(
                GroupMembership.group_id == room.owner_group_id
            )
        )
        for uid in result.all():
            await sio.emit(event, room=f"user:{uid}")
        return
    if room.owner_user_id is not None:
        await sio.emit(event, room=f"user:{room.owner_user_id}")


# =============================================================================
# Room CRUD
# =============================================================================


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(
        NotAuthenticated, Forbidden, InvalidPayload, RoomReadOnly,
        TransferTargetInvalid, UnprocessableContent,
    ),
)
async def create_room(
    session: SessionDep,
    storage: FrameStorageDep,
    sio: SioDep,
    settings: SettingsDep,
    current_user: CurrentUserDep,
    request: RoomCreate,
    response: Response,
) -> RoomCreateResponse:
    """Create or idempotently reuse a room in ``owner_id``'s namespace."""
    name = request.name
    target_visibility = request.visibility or settings.default_room_visibility

    # Step 2: permission gate. Resolve target namespace BEFORE existence check.
    owner_user_id: UUID | None = None
    owner_group_id: UUID | None = None

    if request.owner_id == current_user.id or current_user.is_superuser:
        resolved = await resolve_owner(session, request.owner_id)
        if resolved is None:
            raise Forbidden.exception("Not permitted to create in this namespace")
        kind, _ = resolved
        if kind == OwnerKind.USER:
            owner_user_id = request.owner_id
            if target_visibility == Visibility.GROUP:
                raise InvalidPayload.exception(
                    "GROUP visibility requires a group owner"
                )
        else:
            owner_group_id = request.owner_id
            if target_visibility == Visibility.PRIVATE:
                raise InvalidPayload.exception(
                    "PRIVATE visibility requires a user owner"
                )
    else:
        resolved = await resolve_owner(session, request.owner_id)
        if resolved is None:
            raise Forbidden.exception("Not permitted to create in this namespace")
        kind, _ = resolved
        if kind != OwnerKind.GROUP:
            raise Forbidden.exception("Not permitted to create in this namespace")
        role = await fetch_group_role(session, current_user.id, request.owner_id)
        if role is None:
            raise Forbidden.exception("Not permitted to create in this namespace")
        owner_group_id = request.owner_id
        if target_visibility == Visibility.PRIVATE:
            raise InvalidPayload.exception(
                "PRIVATE visibility requires a user owner"
            )

    # Step 4: insert-or-fetch via the unique index.
    existing = await _load_room_by_address(session, request.owner_id, name)
    if existing is not None:
        frame_count = await storage.get_length(existing.id)
        response.status_code = status.HTTP_200_OK
        return RoomCreateResponse(
            room_id=existing.public_address,
            frame_count=frame_count,
            created=False,
        )

    # Resolve copy_from: @-prefixed presets, room IDs, or server default.
    copy_from = request.copy_from
    if copy_from is None:
        default_room_id = await _get_default_room_id(session)
        copy_from = default_room_id if default_room_id else "@empty"

    presets = {"@empty", "@none"}
    if copy_from.startswith("@") and copy_from not in presets:
        raise UnprocessableContent.exception(
            f"Unknown preset '{copy_from}'. "
            f"Valid presets: {', '.join(sorted(presets))}"
        )

    source_room: Room | None = None
    if not copy_from.startswith("@"):
        source_room = await session.get(Room, copy_from)
        if source_room is not None and await storage.has_mount(copy_from):
            raise RoomReadOnly.exception(
                "Cannot copy from a room with a mounted source"
            )

    room = Room(
        room_name=name,
        description=request.description,
        created_by_id=current_user.id,
        owner_user_id=owner_user_id,
        owner_group_id=owner_group_id,
        visibility=target_visibility,
        step=source_room.step if source_room else 0,
    )
    session.add(room)
    await session.flush()  # populate room.id

    frame_count = 0
    if copy_from == "@none":
        pass
    elif copy_from == "@empty":
        await storage[room.id].extend([{}])
        frame_count = 1
    elif source_room is not None:
        source_frames_or_none = await storage[source_room.id][0:].to_list()
        source_frames = [f for f in source_frames_or_none if f is not None]
        if source_frames:
            await storage[room.id].extend(source_frames)
            frame_count = len(source_frames)
        await _copy_room_state(session, source_room.id, room.id)
    else:
        await storage[room.id].extend([{}])
        frame_count = 1

    if source_room is None:
        _initialize_default_geometries(session, room.id)

    await session.commit()

    await broadcast_room_update(sio, session, storage, room)

    return RoomCreateResponse(
        room_id=room.public_address,
        frame_count=frame_count,
        created=True,
    )


@router.get("")
async def list_rooms(
    session: SessionDep,
    storage: FrameStorageDep,
    current_user: CurrentUserDep,
    my_group_ids: MyGroupIdsDep,
    search: Annotated[str | None, Query(description="Search pattern")] = None,
) -> CollectionResponse[RoomResponse]:
    """List rooms visible to the caller."""
    default_room_id = await _get_default_room_id(session)

    conditions: list[Any] = [
        col(Room.visibility) == Visibility.PUBLIC,
        col(Room.owner_user_id) == current_user.id,
    ]
    if my_group_ids:
        conditions.append(col(Room.owner_group_id).in_(my_group_ids))

    stmt = select(Room).where(or_(*conditions))
    result = await session.exec(stmt)
    rooms = list(result.all())

    room_responses: list[RoomResponse] = []
    for room in rooms:
        if search:
            sl = search.lower()
            if sl not in room.room_name.lower() and (
                room.description is None or sl not in room.description.lower()
            ):
                continue
        frame_count = await storage.get_length(room.id)
        owner_id = room.owner_user_id or room.owner_group_id
        assert owner_id is not None
        resolved = await resolve_owner(session, owner_id)
        assert resolved is not None
        kind, label = resolved
        room_responses.append(
            RoomResponse(
                room_id=room.public_address,
                description=room.description,
                frame_count=frame_count,
                visibility=room.visibility,
                owner_id=owner_id,
                owner_kind=kind.value,
                owner_label=label,
                is_default=(room.id == default_room_id),
            )
        )
    return CollectionResponse(items=room_responses)


@router.get(
    "/{owner_id}/{room_name}",
    responses=problem_responses(RoomNotFound, ShareLinkInvalid),
)
async def get_room(
    session: SessionDep,
    storage: FrameStorageDep,
    access: AccessReadDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
) -> RoomResponse:
    """Get details of a specific room (read-gated)."""
    room = access.room
    frame_count = await storage.get_length(room.id)
    default_room_id = await _get_default_room_id(session)
    own_id = room.owner_user_id or room.owner_group_id
    assert own_id is not None
    resolved = await resolve_owner(session, own_id)
    assert resolved is not None
    kind, label = resolved
    return RoomResponse(
        room_id=room.public_address,
        description=room.description,
        frame_count=frame_count,
        visibility=room.visibility,
        owner_id=own_id,
        owner_kind=kind.value,
        owner_label=label,
        is_default=(room.id == default_room_id),
    )


@router.get(
    "/{owner_id}/{room_name}/presence",
    responses=problem_responses(RoomNotFound),
)
async def get_room_presence(
    redis: RedisDep,
    access: AccessReadDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
) -> PresenceResponse:
    """Get presence (online users) for a room."""
    from uuid import UUID as _UUID

    cameras_raw: dict[str, str] = await redis.hgetall(  # type: ignore[misc]
        RedisKey.room_cameras(access.room.id)
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
    "/{owner_id}/{room_name}/sessions",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_sessions(
    redis: RedisDep,
    access: AccessReadDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
    email: Annotated[str | None, Query(description="Filter by user email")] = None,
) -> SessionsListResponse:
    """List all active frontend sessions in this room."""
    room_id = access.room.id
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
    "/{owner_id}/{room_name}",
    responses=problem_responses(
        RoomNotFound, Forbidden, TransferTargetInvalid, InvalidPayload
    ),
)
async def update_room(
    session: SessionDep,
    storage: FrameStorageDep,
    sio: SioDep,
    access: AccessManageDep,
    current_user: CurrentUserDep,
    updates: RoomPatchRequest,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
) -> RoomPatchResponse:
    """Update room metadata, ownership, or visibility (manage-gated).

    A transfer is expressed as a PATCH that sets exactly one of
    ``owner_user_id`` or ``owner_group_id`` (the other is wiped to
    preserve the XOR invariant). Transferring into a group requires
    the caller to be a member of that group.
    """
    room = access.room
    changed = False

    if updates.description is not None:
        room.description = updates.description
        changed = True

    if updates.owner_user_id is not None and updates.owner_group_id is not None:
        raise InvalidPayload.exception(
            "Set exactly one of owner_user_id or owner_group_id, not both"
        )

    if updates.owner_user_id is not None:
        room.owner_user_id = updates.owner_user_id
        room.owner_group_id = None
        changed = True
    elif updates.owner_group_id is not None:
        if not current_user.is_superuser:
            role = await fetch_group_role(
                session, current_user.id, updates.owner_group_id
            )
            if role is None:
                raise TransferTargetInvalid.exception(
                    "You are not a member of the target group"
                )
        room.owner_group_id = updates.owner_group_id
        room.owner_user_id = None
        changed = True

    if updates.visibility is not None:
        room.visibility = updates.visibility
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

    try:
        await session.commit()
    except IntegrityError as exc:
        await session.rollback()
        raise InvalidPayload.exception(
            "Room update violates visibility/owner invariants"
        ) from exc

    if changed:
        await broadcast_room_update(sio, session, storage, room)

    return RoomPatchResponse()
