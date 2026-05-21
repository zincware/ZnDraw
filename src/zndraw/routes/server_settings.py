"""Server settings REST API endpoints.

Manages server-wide configuration such as the default room
that new rooms copy from when no explicit `copyFrom` is provided.
"""

from fastapi import APIRouter
from pydantic import BaseModel
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.dependencies import (
    AdminUserDep,
    CurrentUserDep,
    FrameStorageDep,
    SessionDep,
    SioDep,
    _load_room_by_address,
)
from zndraw.exceptions import Forbidden, RoomNotFound, problem_responses
from zndraw.models import Room, ServerSettings
from zndraw.routes.rooms import broadcast_room_update
from zndraw.schemas import StatusResponse

router = APIRouter(prefix="/v1/server-settings", tags=["server-settings"])


# =============================================================================
# Schemas
# =============================================================================


class DefaultRoomResponse(BaseModel):
    """Response for the default room setting."""

    room_id: str | None  # composed address: {owner_id}/{room_name}


class DefaultRoomSetRequest(BaseModel):
    """Request to set the default room.

    Accepts either the composed address ``{owner_id}/{room_name}`` or the
    surrogate room UUID (for backwards-compatible CLI usage).
    """

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


async def _resolve_room_by_id_or_address(
    session: AsyncSession, room_id: str
) -> Room | None:
    """Resolve a Room from either a surrogate UUID or a composed address."""
    # Try direct surrogate lookup first
    room = await session.get(Room, room_id)
    if room is not None:
        return room
    # Try composed-address lookup (owner_id/room_name)
    parts = room_id.split("/", 1)
    if len(parts) == 2:
        try:
            from uuid import UUID as _UUID
            owner_uuid = _UUID(parts[0])
        except ValueError:
            return None
        return await _load_room_by_address(session, owner_uuid, parts[1])
    return None


# =============================================================================
# Endpoints
# =============================================================================


@router.get(
    "/default-room",
)
async def get_default_room(
    session: SessionDep,
    _current_user: CurrentUserDep,
) -> DefaultRoomResponse:
    """Get the default room for new room creation."""
    settings = await get_server_settings(session)
    if settings.default_room_id is None:
        return DefaultRoomResponse(room_id=None)
    room = await session.get(Room, settings.default_room_id)
    return DefaultRoomResponse(room_id=room.public_address if room else None)


@router.put(
    "/default-room",
    responses=problem_responses(Forbidden, RoomNotFound),
)
async def set_default_room(
    session: SessionDep,
    storage: FrameStorageDep,
    sio: SioDep,
    _admin: AdminUserDep,
    request: DefaultRoomSetRequest,
) -> DefaultRoomResponse:
    """Set the default room for new room creation.

    Requires admin privileges. The specified room must exist.
    Broadcasts room_update events to notify clients about the change.
    """
    room = await _resolve_room_by_id_or_address(session, request.room_id)
    if room is None:
        raise RoomNotFound.exception(f"Room {request.room_id!r} not found")

    settings = await get_server_settings(session)
    old_default_id = settings.default_room_id

    settings.default_room_id = room.id
    await session.commit()

    # Broadcast full snapshots for affected rooms
    if old_default_id and old_default_id != room.id:
        old_room = await session.get(Room, old_default_id)
        if old_room:
            await broadcast_room_update(sio, session, storage, old_room)

    await broadcast_room_update(sio, session, storage, room)

    return DefaultRoomResponse(room_id=room.public_address)


@router.delete(
    "/default-room",
    responses=problem_responses(Forbidden),
)
async def unset_default_room(
    session: SessionDep,
    storage: FrameStorageDep,
    sio: SioDep,
    _admin: AdminUserDep,
) -> StatusResponse:
    """Unset the default room."""
    settings = await get_server_settings(session)
    old_default_id = settings.default_room_id

    if old_default_id:
        settings.default_room_id = None
        await session.commit()

        old_room = await session.get(Room, old_default_id)
        if old_room:
            await broadcast_room_update(sio, session, storage, old_room)

    return StatusResponse()
