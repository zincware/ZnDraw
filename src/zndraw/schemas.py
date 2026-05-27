from __future__ import annotations

import copy
from datetime import datetime
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar
from uuid import UUID

from pydantic import BaseModel, ConfigDict, Field

from zndraw.access import GroupRole, ShareAccess, Visibility


def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into *base*. Override wins for leaf values."""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


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
    """Request body for POST /v1/rooms."""

    owner: str = Field(pattern=r"^[a-z][a-z0-9-]{2,63}$")
    name: str = Field(pattern=r"^[a-zA-Z0-9\-_]+$", min_length=1, max_length=128)
    description: str | None = None
    copy_from: str | None = None
    visibility: Visibility | None = None


class RoomResponse(BaseModel):
    """Response body for room details — matches frontend Room interface."""

    room_id: str  # composed: {owner}/{room_name}
    id: str  # surrogate UUID (internal, read-only)
    description: str | None = None
    frame_count: int = 0
    visibility: Visibility = Visibility.PUBLIC
    owner: str  # display_name (user) OR group name
    owner_kind: Literal["user", "group"]
    owner_label: str
    is_default: bool = False
    metadata: dict[str, str] | None = None

    model_config = ConfigDict(from_attributes=True)


class RoomCreateResponse(BaseModel):
    """Response for room creation."""

    status: Literal["ok"] = "ok"
    room_id: str  # composed: {owner}/{room_name}
    frame_count: int
    created: bool


class RoomPatchRequest(BaseModel):
    """Request body for PATCH /v1/rooms/{owner}/{room_name}."""

    description: str | None = None
    frame_count: int | None = Field(None, ge=0)
    visibility: Visibility | None = None
    new_owner: str | None = Field(default=None, pattern=r"^[a-z][a-z0-9-]{2,63}$")


class RoomPatchResponse(BaseModel):
    """Response body for PATCH /v1/rooms/{owner}/{room_name}."""

    status: Literal["ok"] = "ok"
    room_id: str  # composed: {owner}/{room_name}


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

    model_config = ConfigDict(from_attributes=True)


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


class SessionItem(BaseModel):
    """A single active frontend session."""

    sid: str
    email: str
    camera_key: str


class SessionsListResponse(BaseModel):
    """Response for listing active frontend sessions in a room."""

    items: list[SessionItem]


# =============================================================================
# Group Schemas
# =============================================================================


class GroupCreate(BaseModel):
    """Request body for POST /v1/groups."""

    name: str = Field(min_length=1, max_length=64, pattern=r"^[a-zA-Z0-9_\-]+$")
    description: str | None = None


class GroupResponse(BaseModel):
    """Response body for group details."""

    id: UUID
    name: str
    description: str | None
    created_at: datetime
    created_by_id: UUID
    my_role: GroupRole | None = None

    model_config = ConfigDict(from_attributes=True)


class GroupPatchRequest(BaseModel):
    """Request body for PATCH /v1/groups/{id}."""

    name: str | None = Field(
        default=None, min_length=1, max_length=64, pattern=r"^[a-zA-Z0-9_\-]+$"
    )
    description: str | None = None


class GroupMemberResponse(BaseModel):
    """Response body for a group member."""

    user_id: UUID
    email: str | None
    role: GroupRole
    joined_at: datetime

    model_config = ConfigDict(from_attributes=True)


class GroupMemberCreateRequest(BaseModel):
    """Request body for POST /v1/groups/{id}/members."""

    user_id: UUID
    role: GroupRole = GroupRole.VIEWER


class GroupMemberPatchRequest(BaseModel):
    """Request body for PATCH /v1/groups/{id}/members/{user_id}."""

    role: GroupRole


# =============================================================================
# Share Link Schemas
# =============================================================================


class ShareLinkCreate(BaseModel):
    """Request body for POST /v1/rooms/{id}/share-links."""

    access: ShareAccess = ShareAccess.VIEW
    expires_at: datetime | None = None


class ShareLinkResponse(BaseModel):
    """Response body for a share link."""

    id: UUID
    room_id: str
    token: str
    access: ShareAccess
    created_by_id: UUID
    created_at: datetime
    expires_at: datetime | None
    revoked_at: datetime | None

    model_config = ConfigDict(from_attributes=True)


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
    """Request to create or update a geometry (full replace)."""

    type: str
    data: dict[str, Any]


class GeometryPatchRequest(BaseModel):
    """Request to partially update a geometry (deep merge)."""

    data: dict[str, Any]


class DefaultCameraResponse(BaseModel):
    """Response for default camera endpoint."""

    default_camera: str | None = None


class DefaultCameraRequest(BaseModel):
    """Request to set/unset default camera."""

    default_camera: str | None = None


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

EDIT_LOCK_REFRESH = 5  # seconds (client-side refresh interval)


class EditLockRequest(BaseModel):
    """Request to acquire or refresh an edit lock."""

    msg: str | None = None


class StatusResponse(BaseModel):
    """Generic success response."""

    status: Literal["ok"] = "ok"


class EditLockResponse(BaseModel):
    """Response for edit lock status."""

    locked: bool
    lock_token: str | None = None
    user_id: str | None = None
    sid: str | None = None
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
    unit: str = "it"


class ProgressPatch(BaseModel):
    """Request to update an existing progress tracker."""

    description: str | None = None
    n: int | None = None
    total: int | None = None
    elapsed: float | None = None
    unit: str | None = None


class ProgressResponse(BaseModel):
    """Response for a single progress tracker."""

    progress_id: str
    description: str
    n: int = 0
    total: int | None = None
    elapsed: float = 0.0
    unit: str = "it"


# =============================================================================
# Preset Schemas
# =============================================================================


class PresetRule(BaseModel):
    """A single rule that targets geometries by pattern.

    Parameters
    ----------
    pattern : str
        fnmatch pattern to match geometry keys.
    geometry_type : str | None
        Optional geometry type filter (e.g., "Sphere", "DirectionalLight").
    config : dict[str, Any]
        Partial geometry config to deep-merge into matching geometries.
    """

    pattern: str = Field(description="fnmatch pattern for geometry keys")
    geometry_type: str | None = Field(
        default=None, description="Optional geometry type filter"
    )
    config: dict[str, Any] = Field(
        description="Partial config to deep-merge into matching geometries"
    )


class Preset(BaseModel):
    """A visual preset — a named collection of rules for styling geometries.

    Parameters
    ----------
    name : str
        Unique name within a room. URL-safe characters only.
    description : str
        Human-readable description.
    rules : list[PresetRule]
        Ordered rules. Later rules override earlier ones for the same geometry.
    created_at : datetime | None
        Server-assigned creation timestamp (read-only, ignored on create/update).
    updated_at : datetime | None
        Server-assigned update timestamp (read-only, ignored on create/update).
    """

    name: str = Field(pattern=r"^[a-zA-Z0-9_-]+$", max_length=64)
    description: str = ""
    rules: list[PresetRule] = Field(default_factory=list)
    created_at: datetime | None = None
    updated_at: datetime | None = None


class PresetsListResponse(BaseModel):
    """Response for listing presets."""

    items: list[Preset]


class PresetApplyResult(BaseModel):
    """Result of applying a preset."""

    geometries_updated: list[str]
