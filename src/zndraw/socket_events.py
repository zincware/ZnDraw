"""Pydantic models for Socket.IO events.

Event names are derived from class names in snake_case:
- RoomJoin -> "room_join"
- SessionJoined -> "session_joined"
"""

from datetime import datetime
from typing import Literal
from uuid import UUID

from pydantic import BaseModel

from zndraw.schemas import ProgressResponse, RoomResponse

# =============================================================================
# Request Models (client -> server)
# =============================================================================


class RoomJoin(BaseModel):
    """Join a room for real-time updates."""

    room_id: str
    client_type: Literal["frontend", "pyclient"] = "frontend"


class RoomLeave(BaseModel):
    """Leave current room."""

    room_id: str


class Heartbeat(BaseModel):
    """Keep presence alive."""

    room_id: str


class UserGet(BaseModel):
    """Get authenticated user info."""


class TypingStart(BaseModel):
    """User started typing."""

    room_id: str


class TypingStop(BaseModel):
    """User stopped typing."""

    room_id: str


# =============================================================================
# Response Models (server -> client, for RPC calls)
# =============================================================================


class RoomJoinResponse(BaseModel):
    """Response for room join."""

    room_id: str
    session_id: str
    step: int
    frame_count: int
    locked: bool
    camera_key: str | None = None
    progress_trackers: dict[str, ProgressResponse] = {}


class RoomLeaveResponse(BaseModel):
    """Response for room leave."""

    room_id: str


class HeartbeatResponse(BaseModel):
    """Response for heartbeat."""

    status: Literal["ok"] = "ok"


class UserGetResponse(BaseModel):
    """Response for user get."""

    id: UUID
    email: str
    is_superuser: bool


class TypingResponse(BaseModel):
    """Response for typing events."""

    status: Literal["ok"] = "ok"


# =============================================================================
# Broadcast Models (server -> clients)
# =============================================================================


class SessionJoined(BaseModel):
    """Broadcast when a session joins a room."""

    room_id: str
    user_id: UUID
    sid: str
    email: str | None = None


class SessionLeft(BaseModel):
    """Broadcast when a session leaves a room."""

    room_id: str
    user_id: UUID
    sid: str


class FrameUpdate(BaseModel):
    """Broadcast when current frame changes."""

    room_id: str
    frame: int


class FramesInvalidate(BaseModel):
    """Broadcast when frames need refresh."""

    room_id: str
    action: Literal["add", "delete", "modify", "clear"]
    indices: list[int] | None = None
    count: int | None = None
    reason: str | None = None


class FrameSelectionUpdate(BaseModel):
    """Broadcast when frame selection changes."""

    indices: list[int]


class GeometryInvalidate(BaseModel):
    """Broadcast when geometries changed."""

    room_id: str
    operation: Literal["set", "delete"]
    key: str


class SelectionInvalidate(BaseModel):
    """Broadcast when selections changed."""

    room_id: str


class SelectionGroupsInvalidate(BaseModel):
    """Broadcast when selection groups changed."""

    room_id: str


class BookmarksInvalidate(BaseModel):
    """Broadcast when bookmarks changed."""

    room_id: str
    index: int
    operation: Literal["set", "delete"]


class FigureInvalidate(BaseModel):
    """Broadcast when a figure changed."""

    room_id: str
    key: str
    operation: Literal["set", "delete"]


class LockUpdate(BaseModel):
    """Broadcast when a room's edit lock changes."""

    room_id: str
    action: Literal["acquired", "refreshed", "released"]
    user_id: str | None = None
    msg: str | None = None


class RoomUpdate(RoomResponse):
    """Full room snapshot broadcast on any room state change."""


class MessageNew(BaseModel):
    """Broadcast new message in room."""

    id: int
    room_id: str
    user_id: UUID
    content: str
    created_at: datetime
    updated_at: datetime | None = None
    email: str | None = None


class MessageEdited(BaseModel):
    """Broadcast message was edited."""

    id: int
    room_id: str
    content: str
    updated_at: datetime


class ScreenshotRequest(BaseModel):
    """Sent to a specific frontend session to capture a screenshot."""

    screenshot_id: int
    upload_url: str


class Typing(BaseModel):
    """Broadcast typing indicator."""

    room_id: str
    user_id: UUID
    email: str | None = None
    is_typing: bool


class ProgressStart(BaseModel):
    """Broadcast when a new progress tracker is created."""

    progress_id: str
    description: str
    unit: str = "it"


class ProgressUpdate(BaseModel):
    """Broadcast when a progress tracker is updated."""

    progress_id: str
    description: str | None = None
    n: int | None = None
    total: int | None = None
    elapsed: float | None = None
    unit: str | None = None


class ProgressComplete(BaseModel):
    """Broadcast when a progress tracker finishes."""

    progress_id: str
