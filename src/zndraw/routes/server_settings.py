"""Server settings REST API endpoints.

Manages server-wide configuration such as the default room
that new rooms copy from when no explicit `copyFrom` is provided.
"""

from fastapi import APIRouter, status
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.dependencies import (
    AdminUserDep,
    CurrentUserDep,
    SessionDep,
    SioDep,
    StorageDep,
)
from zndraw.exceptions import Forbidden, RoomNotFound, problem_responses
from zndraw.models import Room, ServerSettings
from zndraw.routes.rooms import build_room_update
from zndraw.schemas import StatusResponse

router = APIRouter(prefix="/v1/server-settings", tags=["server-settings"])


# =============================================================================
# Schemas
# =============================================================================


class DefaultRoomResponse(BaseModel):
    """Response for the default room setting."""

    room_id: str | None


class DefaultRoomSetRequest(BaseModel):
    """Request to set the default room."""

    room_id: str


# =============================================================================
# Helpers
# =============================================================================


async def get_server_settings(session: AsyncSession) -> ServerSettings:
    """Get or create the singleton ServerSettings row."""
    settings = await session.get(ServerSettings, 1)
    if settings is None:
        settings = ServerSettings(id=1)
        session.add(settings)
        await session.flush()
    return settings


# =============================================================================
# Endpoints
# =============================================================================


@router.get(
    "/default-room",
    response_model=DefaultRoomResponse,
)
async def get_default_room(
    session: SessionDep,
    _current_user: CurrentUserDep,
) -> DefaultRoomResponse:
    """Get the default room for new room creation.

    Returns the room ID that new rooms copy from when no explicit
    ``copyFrom`` is provided. Returns null if no default is set.
    """
    settings = await get_server_settings(session)
    return DefaultRoomResponse(room_id=settings.default_room_id)


@router.put(
    "/default-room",
    response_model=DefaultRoomResponse,
    responses=problem_responses(Forbidden, RoomNotFound),
)
async def set_default_room(
    session: SessionDep,
    storage: StorageDep,
    sio: SioDep,
    _admin: AdminUserDep,
    request: DefaultRoomSetRequest,
) -> DefaultRoomResponse:
    """Set the default room for new room creation.

    Requires admin privileges. The specified room must exist.
    Broadcasts room_update events to notify clients about the change.
    """
    room = await session.get(Room, request.room_id)
    if room is None:
        raise RoomNotFound.exception(f"Room with id {request.room_id} not found")

    settings = await get_server_settings(session)
    old_default_id = settings.default_room_id

    settings.default_room_id = request.room_id
    await session.commit()

    # Broadcast full snapshots for affected rooms
    if old_default_id and old_default_id != request.room_id:
        old_room = await session.get(Room, old_default_id)
        if old_room:
            event = await build_room_update(session, storage, old_room)
            await sio.emit(event, room="room:@overview")
            await sio.emit(event, room=f"room:{old_default_id}")

    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")
    await sio.emit(event, room=f"room:{request.room_id}")

    return DefaultRoomResponse(room_id=request.room_id)


@router.delete(
    "/default-room",
    responses=problem_responses(Forbidden),
)
async def unset_default_room(
    session: SessionDep,
    storage: StorageDep,
    sio: SioDep,
    _admin: AdminUserDep,
) -> StatusResponse:
    """Unset the default room.

    Requires admin privileges. After this, new rooms without an explicit
    ``copyFrom`` will start with a single empty frame.
    """
    settings = await get_server_settings(session)
    old_default_id = settings.default_room_id

    if old_default_id:
        settings.default_room_id = None
        await session.commit()

        old_room = await session.get(Room, old_default_id)
        if old_room:
            event = await build_room_update(session, storage, old_room)
            await sio.emit(event, room="room:@overview")
            await sio.emit(event, room=f"room:{old_default_id}")

    return StatusResponse()
