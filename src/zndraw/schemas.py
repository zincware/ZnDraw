from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar
from uuid import UUID

from pydantic import BaseModel, Field

from zndraw.models import MemberRole

T = TypeVar("T")


class CollectionResponse(BaseModel, Generic[T]):
    """Non-paginated collection envelope."""

    items: list[T]


class OffsetPage(BaseModel, Generic[T]):
    """Offset-paginated collection envelope."""

    items: list[T]
    total: int
    limit: int
    offset: int


if TYPE_CHECKING:
    import plotly.graph_objects as go

# =============================================================================
# Room Schemas
# =============================================================================


class RoomCreate(BaseModel):
    """Request body for creating a new room."""

    room_id: str  # UUID string
    description: str | None = None
    copy_from: str | None = None  # Room ID, or @-prefixed preset (@empty, @none)


class RoomResponse(BaseModel):
    """Response body for room details - matches frontend Room interface."""

    id: str
    description: str | None = None
    frame_count: int = 0
    locked: bool = False
    is_default: bool = False
    metadata: dict[str, str] | None = None

    model_config = {"from_attributes": True}


class RoomCreateResponse(BaseModel):
    """Response for room creation."""

    status: Literal["ok"] = "ok"
    room_id: str
    frame_count: int
    created: bool


class RoomPatchRequest(BaseModel):
    """Request body for PATCH /rooms/{room_id}."""

    description: str | None = None
    locked: bool | None = None
    frame_count: int | None = Field(None, ge=0)


class RoomPatchResponse(BaseModel):
    """Response body for PATCH /rooms/{room_id}."""

    status: Literal["ok"] = "ok"


class RoomMemberResponse(BaseModel):
    """Response body for room member details."""

    user_id: UUID
    email: str | None
    role: MemberRole
    joined_at: datetime

    model_config = {"from_attributes": True}


class MessageCreate(BaseModel):
    """Request body for sending a message."""

    content: str = Field(min_length=1)


class MessageEditRequest(BaseModel):
    """Request body for editing a message."""

    content: str = Field(min_length=1)


class MessageResponse(BaseModel):
    """Response body for message details."""

    id: int
    room_id: str
    user_id: UUID
    content: str
    created_at: datetime
    updated_at: datetime | None = None
    email: str | None = None

    model_config = {"from_attributes": True}


class MessagesMetadata(BaseModel):
    """Pagination metadata for message listing."""

    has_more: bool
    total_count: int
    oldest_timestamp: int | None = None  # unix ms
    newest_timestamp: int | None = None  # unix ms


class MessagesResponse(BaseModel):
    """Response body for listing messages with pagination."""

    items: list[MessageResponse]
    metadata: MessagesMetadata


class PresenceSessionResponse(BaseModel):
    """Response body for a single session (connection) in presence."""

    sid: str
    user_id: UUID
    email: str | None


class PresenceResponse(BaseModel):
    """Response body for room presence (online sessions).

    Each session (tab/connection) is listed individually.
    The same user may appear multiple times with different SIDs.
    """

    items: list[PresenceSessionResponse]


class SessionsListResponse(BaseModel):
    """Response for listing user's own frontend sessions."""

    items: list[str]


class SessionSettingsResponse(BaseModel):
    """Response for session settings with JSON schema."""

    schema_: dict[str, Any] = Field(alias="schema")
    data: dict[str, Any]

    model_config = {"populate_by_name": True}


# =============================================================================
# Frame Schemas
# =============================================================================


class FrameResponse(BaseModel):
    """Single frame response."""

    index: int
    data: dict[str, Any]


class FrameBulkResponse(BaseModel):
    """Bulk frame response for range queries."""

    frames: list[dict[str, Any]]
    total: int = Field(description="Total frames in room")
    start: int = Field(description="Start index (inclusive)")
    stop: int = Field(description="Stop index (exclusive)")


class FrameCreateRequest(BaseModel):
    """Request to append frames."""

    frames: list[dict[str, Any]] = Field(min_length=1, max_length=1000)


class FrameUpdateRequest(BaseModel):
    """Request to update a single frame."""

    data: dict[str, Any]


class FrameMergeResponse(BaseModel):
    """Response for partial frame update (PATCH)."""

    index: int = Field(description="Index of the updated frame")
    updated_keys: list[str] = Field(description="Keys that were updated")


class PropertyMeta(BaseModel):
    """Metadata for a single frame property (key)."""

    dtype: str
    shape: list[int] | None = None
    type: Literal["array", "scalar"]


class FrameMetadataResponse(BaseModel):
    """Response for frame metadata: key -> shape/dtype per key."""

    frame_id: int
    metadata: dict[str, PropertyMeta]
    source_room: str


# =============================================================================
# Frame Selection Schemas
# =============================================================================


class FrameSelectionResponse(BaseModel):
    """Response for GET frame-selection endpoint."""

    frame_selection: list[int] | None = None


class FrameSelectionUpdateRequest(BaseModel):
    """Request to update frame selection."""

    indices: list[int] = Field(default_factory=list)


class FrameSelectionUpdateResponse(BaseModel):
    """Response for PUT frame-selection endpoint."""

    success: bool = True


# =============================================================================
# Step Schemas
# =============================================================================


class StepResponse(BaseModel):
    """Response for GET step endpoint."""

    step: int
    total_frames: int


class StepUpdateResponse(BaseModel):
    """Response for PUT step endpoint."""

    success: bool = True
    step: int


class StepUpdateRequest(BaseModel):
    """Request to update step."""

    step: int = Field(ge=0)


# =============================================================================
# Geometry Schemas
# =============================================================================


class GeometryData(BaseModel):
    """Single geometry data."""

    type: str
    data: dict[str, Any]
    selection: list[int] = []


class GeometryResponse(BaseModel):
    """Response for single geometry."""

    key: str
    geometry: GeometryData


class GeometryTypesInfo(BaseModel):
    """Geometry type schemas and defaults from Pydantic models."""

    schemas: dict[str, Any]
    defaults: dict[str, Any]


class GeometriesResponse(BaseModel):
    """Response for listing all geometries."""

    items: dict[str, GeometryData]
    types: GeometryTypesInfo | None = None


class GeometryCreateRequest(BaseModel):
    """Request to create or update a geometry."""

    type: str
    data: dict[str, Any]


# =============================================================================
# Selection Schemas
# =============================================================================


class GeometrySelectionResponse(BaseModel):
    """Response for single geometry selection."""

    key: str
    selection: list[int]


class SelectionUpdateRequest(BaseModel):
    """Request to update selection."""

    indices: list[int]


class SelectionGroupResponse(BaseModel):
    """Response for selection group."""

    group: dict[str, list[int]]


class SelectionGroupUpdateRequest(BaseModel):
    """Request to update selection group."""

    selections: dict[str, list[int]]


class SelectionGroupsListResponse(BaseModel):
    """Response for listing all selection groups."""

    items: dict[str, dict[str, list[int]]]


# =============================================================================
# Bookmark Schemas
# =============================================================================


class BookmarkResponse(BaseModel):
    """Response for single bookmark."""

    index: int
    label: str


class BookmarksResponse(BaseModel):
    """Response for all bookmarks."""

    items: dict[str, str]  # frame_index (as string) -> label


class BookmarkCreateRequest(BaseModel):
    """Request to create/update bookmark."""

    label: str = Field(min_length=1)


# =============================================================================
# Figure Schemas
# =============================================================================


class FigureData(BaseModel):
    """Pydantic model for Plotly figure wire format.

    Delegates entirely to Plotly's own ``to_json`` / ``from_json``
    for serialization — no custom encoding/decoding needed.
    """

    type: str = "plotly"
    data: str  # Plotly JSON string

    @classmethod
    def from_figure(cls, fig: go.Figure) -> FigureData:
        """Serialize a Plotly figure to the wire format."""
        return cls(data=fig.to_json())

    def to_figure(self) -> go.Figure:
        """Deserialize back to a Plotly figure."""
        import plotly.io as pio

        return pio.from_json(self.data)


class FigureResponse(BaseModel):
    """Response for single figure."""

    key: str
    figure: FigureData


class FigureCreateRequest(BaseModel):
    """Request to create/update figure."""

    figure: FigureData


class FigureCreateResponse(BaseModel):
    """Response for figure creation/update."""

    key: str
    created: bool


# =============================================================================
# Active Camera Schemas
# =============================================================================


class ActiveCameraRequest(BaseModel):
    """Request body for setting active camera."""

    active_camera: str


class ActiveCameraResponse(BaseModel):
    """Response body for active camera."""

    active_camera: str | None


# =============================================================================
# Edit Lock Schemas
# =============================================================================

EDIT_LOCK_TTL = 10  # seconds
EDIT_LOCK_REFRESH = 5  # seconds


class EditLockRequest(BaseModel):
    """Request to acquire or refresh an edit lock."""

    msg: str | None = None


class StatusResponse(BaseModel):
    """Generic success response."""

    status: Literal["ok"] = "ok"


class EditLockResponse(BaseModel):
    """Response for edit lock status."""

    locked: bool
    user_id: str | None = None
    msg: str | None = None
    acquired_at: float | None = None
    ttl: int | None = None


# =============================================================================
# Screenshot Schemas
# =============================================================================


class ScreenshotResponse(BaseModel):
    """Response body for a single screenshot.

    The ``data`` field contains the raw image bytes as a base64-encoded string
    (only present for completed screenshots when the image is included).
    """

    id: int
    room_id: str
    format: str
    size: int
    width: int | None
    height: int | None
    status: Literal["pending", "completed"]
    created_by_id: UUID
    created_at: datetime
    data: str | None = None


class ScreenshotListItem(BaseModel):
    """Summary item for screenshot listing (no binary data)."""

    id: int
    room_id: str
    format: str
    size: int
    width: int | None
    height: int | None
    created_by_id: UUID
    created_at: datetime


class ScreenshotCaptureCreate(BaseModel):
    """Request body for programmatic screenshot capture."""

    session_id: str


# =============================================================================
# Progress Schemas
# =============================================================================


class ProgressCreate(BaseModel):
    """Request to start a new progress tracker."""

    progress_id: str
    description: str


class ProgressPatch(BaseModel):
    """Request to update an existing progress tracker."""

    description: str | None = None
    progress: float | None = None


class ProgressResponse(BaseModel):
    """Response for a single progress tracker."""

    progress_id: str
    description: str
    progress: float | None = None
