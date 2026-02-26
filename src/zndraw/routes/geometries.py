"""Geometry REST API endpoints for room geometry management."""

import json
from functools import lru_cache
from typing import Any

from fastapi import APIRouter, status
from pydantic import ValidationError
from sqlmodel import select

from zndraw.dependencies import (
    CurrentUserDep,
    OptionalUserDep,
    RedisDep,
    SessionDep,
    SioDep,
    WritableGeometryDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    Forbidden,
    GeometryNotFound,
    InvalidPayload,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    problem_responses,
)
from zndraw.geometries import geometries as geometry_models
from zndraw.models import Room, RoomGeometry
from zndraw.redis import RedisKey
from zndraw.schemas import (
    DefaultCameraRequest,
    DefaultCameraResponse,
    GeometriesResponse,
    GeometryCreateRequest,
    GeometryData,
    GeometryResponse,
    GeometrySelectionResponse,
    GeometryTypesInfo,
    SelectionUpdateRequest,
    StatusResponse,
)
from zndraw.socket_events import (
    DefaultCameraInvalidate,
    GeometryInvalidate,
    SelectionInvalidate,
)


@lru_cache(maxsize=1)
def _get_geometry_types_info() -> GeometryTypesInfo:
    """Get cached geometry type schemas and defaults from Pydantic models."""
    schemas: dict[str, Any] = {}
    defaults: dict[str, Any] = {}

    for name, model in geometry_models.items():
        schemas[name] = model.model_json_schema()
        defaults[name] = model().model_dump()

    return GeometryTypesInfo(schemas=schemas, defaults=defaults)


def _row_to_geometry_data(row: RoomGeometry) -> GeometryData:
    """Convert a RoomGeometry row to a GeometryData response."""
    data = json.loads(row.config)
    selection = json.loads(row.selection) if row.selection else []
    return GeometryData(type=row.type, data=data, selection=selection)


async def _get_session_cameras_for_room(
    redis: RedisDep, room_id: str
) -> dict[str, GeometryData]:
    """Read all session cameras from the room's Redis hash."""
    hash_key = RedisKey.room_cameras(room_id)
    raw_map: dict[str, str] = await redis.hgetall(hash_key)  # type: ignore[misc]
    cameras: dict[str, GeometryData] = {}
    for cam_key, raw in raw_map.items():
        entry = json.loads(raw)
        cameras[cam_key] = GeometryData(type="Camera", data=entry["data"])
    return cameras


router = APIRouter(prefix="/v1/rooms/{room_id}/geometries", tags=["geometries"])


@router.get(
    "",
    responses=problem_responses(RoomNotFound),
)
async def list_geometries(
    session: SessionDep,
    redis: RedisDep,
    _current_user: OptionalUserDep,
    room_id: str,
) -> GeometriesResponse:
    """List all geometries in a room.

    Merges SQL geometries with Redis session cameras.
    """
    await verify_room(session, room_id)

    # SQL geometries
    result = await session.execute(
        select(RoomGeometry).where(RoomGeometry.room_id == room_id)
    )
    rows = result.scalars().all()
    geometries = {row.key: _row_to_geometry_data(row) for row in rows}

    # Merge Redis session cameras
    session_cameras = await _get_session_cameras_for_room(redis, room_id)
    geometries.update(session_cameras)

    return GeometriesResponse(
        items=geometries,
        types=_get_geometry_types_info(),
    )


@router.get(
    "/{key:path}/selection",
    responses=problem_responses(RoomNotFound, GeometryNotFound),
)
async def get_geometry_selection(
    session: SessionDep,
    _current_user: OptionalUserDep,
    room_id: str,
    key: str,
) -> GeometrySelectionResponse:
    """Get selection for a specific geometry."""
    await verify_room(session, room_id)
    row = await session.get(RoomGeometry, (room_id, key))
    if row is None:
        raise GeometryNotFound.exception(f"Geometry '{key}' not found")
    indices = json.loads(row.selection) if row.selection else []
    return GeometrySelectionResponse(key=key, selection=indices)


@router.get(
    "/{key:path}",
    responses=problem_responses(RoomNotFound, GeometryNotFound),
)
async def get_geometry(
    session: SessionDep,
    redis: RedisDep,
    _current_user: OptionalUserDep,
    room_id: str,
    key: str,
) -> GeometryResponse:
    """Get a single geometry by key."""
    await verify_room(session, room_id)

    # Try Redis hash first (session cameras), then SQL
    hash_key = RedisKey.room_cameras(room_id)
    raw = await redis.hget(hash_key, key)  # type: ignore[misc]
    if raw is not None:
        entry = json.loads(raw)
        return GeometryResponse(
            key=key,
            geometry=GeometryData(type="Camera", data=entry["data"]),
        )

    row = await session.get(RoomGeometry, (room_id, key))
    if row is None:
        raise GeometryNotFound.exception(f"Geometry '{key}' not found")
    return GeometryResponse(key=key, geometry=_row_to_geometry_data(row))


@router.put(
    "/{key:path}/selection",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, GeometryNotFound, RoomLocked, Forbidden
    ),
)
async def update_geometry_selection(
    session: SessionDep,
    sio: SioDep,
    _room: WritableGeometryDep,
    room_id: str,
    key: str,
    request: SelectionUpdateRequest,
) -> StatusResponse:
    """Update selection for a specific geometry.

    Updates only the selection column — no read-modify-write of config.
    """
    row = await session.get(RoomGeometry, (room_id, key))
    if row is None:
        raise GeometryNotFound.exception(f"Geometry '{key}' not found")

    row.selection = json.dumps(request.indices)
    await session.commit()

    await sio.emit(SelectionInvalidate(room_id=room_id), room=room_channel(room_id))
    return StatusResponse()


@router.put(
    "/{key:path}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, RoomLocked, Forbidden),
)
async def upsert_geometry(
    session: SessionDep,
    redis: RedisDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    _geo_info: WritableGeometryDep,
    room_id: str,
    key: str,
    request: GeometryCreateRequest,
) -> StatusResponse:
    """Create or update a geometry.

    Ownership is enforced via DI (WritableGeometryDep); the new_owner
    payload check prevents non-superusers from claiming for others.
    """
    new_owner = request.data.get("owner")
    if (
        new_owner is not None
        and new_owner != str(current_user.id)
        and not current_user.is_superuser
    ):
        raise InvalidPayload.exception("owner must be your own user ID or null")

    # Validate via Pydantic model if type is known
    model_cls = geometry_models.get(request.type)
    if model_cls is not None:
        try:
            config_json = model_cls(**request.data).model_dump_json()
        except ValidationError as exc:
            raise InvalidPayload.exception(str(exc)) from exc
    else:
        config_json = json.dumps(request.data)

    # Try Redis hash first (session cameras), then SQL
    hash_key = RedisKey.room_cameras(room_id)
    raw = await redis.hget(hash_key, key)  # type: ignore[misc]
    if raw is not None:
        entry = json.loads(raw)
        entry["data"] = json.loads(config_json)
        await redis.hset(hash_key, key, json.dumps(entry))  # type: ignore[misc]
    else:
        row = await session.get(RoomGeometry, (room_id, key))
        if row is None:
            row = RoomGeometry(
                room_id=room_id,
                key=key,
                type=request.type,
                config=config_json,
            )
            session.add(row)
        else:
            row.type = request.type
            row.config = config_json
        await session.commit()

    await sio.emit(
        GeometryInvalidate(room_id=room_id, operation="set", key=key),
        room=room_channel(room_id),
    )

    return StatusResponse()


@router.delete(
    "/{key:path}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, RoomLocked, Forbidden),
)
async def delete_geometry(
    session: SessionDep,
    redis: RedisDep,
    sio: SioDep,
    _room: WritableGeometryDep,
    room_id: str,
    key: str,
) -> StatusResponse:
    """Delete a geometry."""
    # Reject deletion if any session is attached to this camera
    all_active: dict = await redis.hgetall(RedisKey.active_cameras(room_id))  # type: ignore[misc]
    for attached_key in all_active.values():
        if attached_key == key:
            raise Forbidden.exception("Cannot delete a camera that is in use")

    # Try Redis hash first (session cameras), then SQL
    hash_key = RedisKey.room_cameras(room_id)
    deleted_from_redis = await redis.hdel(hash_key, key)  # type: ignore[misc]

    if not deleted_from_redis:
        row = await session.get(RoomGeometry, (room_id, key))
        if row is not None:
            await session.delete(row)
            await session.commit()

    # Clear default camera if this geometry was the default
    room = await session.get(Room, room_id)
    if room is not None and room.default_camera == key:
        room.default_camera = None
        session.add(room)
        await session.commit()
        await sio.emit(
            DefaultCameraInvalidate(room_id=room_id, default_camera=None),
            room=room_channel(room_id),
        )

    await sio.emit(
        GeometryInvalidate(room_id=room_id, operation="delete", key=key),
        room=room_channel(room_id),
    )

    return StatusResponse()


# =============================================================================
# Default Camera
# =============================================================================

default_camera_router = APIRouter(prefix="/v1/rooms/{room_id}", tags=["geometries"])


@default_camera_router.get(
    "/default-camera",
    responses=problem_responses(RoomNotFound),
)
async def get_default_camera(
    session: SessionDep,
    _user: OptionalUserDep,
    room_id: str,
) -> DefaultCameraResponse:
    """Get the default camera geometry key for a room."""
    room = await verify_room(session, room_id)
    return DefaultCameraResponse(default_camera=room.default_camera)


@default_camera_router.put(
    "/default-camera",
    responses=problem_responses(
        RoomNotFound, NotAuthenticated, GeometryNotFound, InvalidPayload
    ),
)
async def set_default_camera(
    session: SessionDep,
    sio: SioDep,
    _user: CurrentUserDep,
    room_id: str,
    body: DefaultCameraRequest,
) -> DefaultCameraResponse:
    """Set or unset the default camera for a room."""
    room = await verify_room(session, room_id)

    if body.default_camera is not None:
        row = await session.get(RoomGeometry, (room_id, body.default_camera))
        if row is None:
            raise GeometryNotFound.exception(
                f"Geometry '{body.default_camera}' not found"
            )
        if row.type != "Camera":
            raise InvalidPayload.exception(
                f"Geometry '{body.default_camera}' is type '{row.type}', not Camera"
            )

    room.default_camera = body.default_camera
    session.add(room)
    await session.commit()

    await sio.emit(
        DefaultCameraInvalidate(room_id=room_id, default_camera=room.default_camera),
        room=room_channel(room_id),
    )

    return DefaultCameraResponse(default_camera=room.default_camera)
