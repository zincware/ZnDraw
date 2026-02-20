"""Utility and system routes.

Handles version, global settings, health checks, and other system endpoints.
"""

import json

from fastapi import APIRouter, Depends
from pydantic import BaseModel

import zndraw
from zndraw.config import Settings, get_zndraw_settings
from zndraw.dependencies import (
    CurrentUserDep,
    RedisDep,
    SessionDep,
    SioDep,
    VerifiedSessionDep,
    WritableRoomDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    InvalidPayload,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    SessionNotFound,
    problem_responses,
)
from zndraw.redis import RedisKey
from zndraw.schemas import (
    ActiveCameraRequest,
    ActiveCameraResponse,
    FrameSelectionResponse,
    FrameSelectionUpdateRequest,
    FrameSelectionUpdateResponse,
    SessionSettingsResponse,
    StatusResponse,
)
from zndraw.settings import RoomConfig
from zndraw.socket_events import FrameSelectionUpdate

router = APIRouter(prefix="/v1", tags=["utility"])


class VersionResponse(BaseModel):
    """Response model for version endpoint."""

    version: str


class SiMGenSettings(BaseModel):
    """SiMGen feature settings."""

    enabled: bool = False


class GlobalSettings(BaseModel):
    """Global application settings for frontend configuration."""

    simgen: SiMGenSettings


class HealthResponse(BaseModel):
    """Response model for health endpoint."""

    status: str


@router.get("/health", response_model=HealthResponse)
async def get_health() -> HealthResponse:
    """Health check endpoint.

    Returns OK status if the server is running.
    """
    return HealthResponse(status="ok")


@router.get("/version", response_model=VersionResponse)
async def get_version() -> VersionResponse:
    """Get the ZnDraw server version.

    Returns the current server version from the package metadata.
    This is used by the frontend to check version compatibility.
    """
    return VersionResponse(version=zndraw.__version__)


@router.get("/config/global-settings", response_model=GlobalSettings)
async def get_global_settings() -> GlobalSettings:
    """Get global settings for frontend feature flags.

    Returns configuration for optional features like SiMGen integration.
    """
    # TODO: Load from config when simgen support is added
    return GlobalSettings(simgen=SiMGenSettings(enabled=False))


@router.get(
    "/rooms/{room_id}/sessions/{session_id}/settings",
    responses=problem_responses(NotAuthenticated, SessionNotFound),
)
async def get_session_settings(
    redis: RedisDep, room_id: str, session_id: VerifiedSessionDep
) -> SessionSettingsResponse:
    """Get session-specific settings with JSON schema."""
    key = RedisKey.session_settings(room_id, session_id)
    raw = await redis.get(key)  # type: ignore[misc]
    config = RoomConfig.model_validate_json(raw) if raw else RoomConfig()
    return SessionSettingsResponse(
        schema=RoomConfig.model_json_schema(),
        data=config.model_dump(),
    )


@router.put(
    "/rooms/{room_id}/sessions/{session_id}/settings",
    responses=problem_responses(NotAuthenticated, SessionNotFound),
)
async def update_session_settings(
    redis: RedisDep,
    room_id: str,
    session_id: VerifiedSessionDep,
    body: RoomConfig,
    settings: Settings = Depends(get_zndraw_settings),
) -> StatusResponse:
    """Update session-specific settings."""
    key = RedisKey.session_settings(room_id, session_id)
    await redis.set(key, body.model_dump_json(), ex=settings.presence_ttl)  # type: ignore[misc]
    return StatusResponse()


@router.get(
    "/rooms/{room_id}/frame-selection",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def get_frame_selection(
    session: SessionDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> FrameSelectionResponse:
    """Get selected frame indices for a room."""
    room = await verify_room(session, room_id)
    indices = json.loads(room.frame_selection) if room.frame_selection else None
    return FrameSelectionResponse(frame_selection=indices)


@router.put(
    "/rooms/{room_id}/frame-selection",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, InvalidPayload
    ),
)
async def update_frame_selection(
    session: SessionDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    body: FrameSelectionUpdateRequest,
) -> FrameSelectionUpdateResponse:
    """Set selected frame indices for a room.

    Broadcasts frame_selection_update to the room.
    """
    if any(i < 0 for i in body.indices):
        raise InvalidPayload.exception("All indices must be non-negative")

    room.frame_selection = json.dumps(body.indices) if body.indices else None
    await session.commit()

    await sio.emit(
        FrameSelectionUpdate(indices=body.indices),
        room=room_channel(room_id),
    )

    return FrameSelectionUpdateResponse()


@router.get(
    "/rooms/{room_id}/sessions/{session_id}/active-camera",
    responses=problem_responses(NotAuthenticated, SessionNotFound),
)
async def get_active_camera(
    redis: RedisDep, room_id: str, session_id: VerifiedSessionDep
) -> ActiveCameraResponse:
    """Get active camera key for a session."""
    key = await redis.hget(RedisKey.active_cameras(room_id), session_id)  # type: ignore[misc]
    return ActiveCameraResponse(active_camera=key)


@router.put(
    "/rooms/{room_id}/sessions/{session_id}/active-camera",
    responses=problem_responses(NotAuthenticated, SessionNotFound),
)
async def set_active_camera(
    redis: RedisDep,
    room_id: str,
    session_id: VerifiedSessionDep,
    body: ActiveCameraRequest,
) -> ActiveCameraResponse:
    """Set active camera key for a session."""
    await redis.hset(RedisKey.active_cameras(room_id), session_id, body.active_camera)  # type: ignore[misc]
    return ActiveCameraResponse(active_camera=body.active_camera)
