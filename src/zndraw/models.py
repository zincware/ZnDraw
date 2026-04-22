import uuid as uuid_mod
from datetime import UTC, datetime
from uuid import UUID

from sqlalchemy import (
    CheckConstraint,
    Column,
    ForeignKey,
    String,
    TypeDecorator,
    UniqueConstraint,
)
from sqlalchemy.types import DateTime
from sqlmodel import Field, SQLModel

from zndraw.access import GroupRole, ShareAccess, Visibility
from zndraw_joblib.models import Job, Task, Worker, WorkerJobLink  # noqa: F401


class UTCDateTime(TypeDecorator):
    """SQLAlchemy type that ensures datetimes are always UTC-aware.

    SQLite strips timezone info on storage. This type decorator
    re-attaches UTC on load so consumers never see naive datetimes.
    """

    impl = DateTime(timezone=True)
    cache_ok = True

    def process_result_value(
        self, value: datetime | None, _dialect: object
    ) -> datetime | None:
        if value is not None and value.tzinfo is None:
            return value.replace(tzinfo=UTC)
        return value


class Room(SQLModel, table=True):
    """Room with polymorphic ownership and three-value visibility."""

    __table_args__ = (
        CheckConstraint(
            "(owner_user_id IS NOT NULL) <> (owner_group_id IS NOT NULL)",
            name="room_owner_exactly_one",
        ),
        CheckConstraint(
            "(visibility = 'PRIVATE' AND owner_user_id IS NOT NULL) OR "
            "(visibility = 'GROUP'   AND owner_group_id IS NOT NULL) OR "
            "(visibility = 'PUBLIC')",
            name="room_visibility_matches_owner",
        ),
    )

    id: str = Field(default_factory=lambda: str(uuid_mod.uuid4()), primary_key=True)
    description: str | None = None
    created_by_id: UUID | None = Field(default=None, index=True)
    created_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )
    owner_user_id: UUID | None = Field(default=None, foreign_key="user.id", index=True)
    owner_group_id: UUID | None = Field(default=None, foreign_key="group.id", index=True)
    visibility: Visibility = Field(default=Visibility.PUBLIC)
    step: int = Field(default=0)
    frame_selection: str | None = Field(default=None)
    default_camera: str | None = Field(default=None)


class Group(SQLModel, table=True):
    id: UUID = Field(default_factory=uuid_mod.uuid4, primary_key=True)
    name: str = Field(unique=True, index=True)
    description: str | None = None
    created_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )
    created_by_id: UUID = Field(foreign_key="user.id", index=True)


class GroupMembership(SQLModel, table=True):
    __table_args__ = (UniqueConstraint("group_id", "user_id"),)

    id: int | None = Field(default=None, primary_key=True)
    group_id: UUID = Field(foreign_key="group.id", index=True)
    user_id: UUID = Field(foreign_key="user.id", index=True)
    role: GroupRole = Field(default=GroupRole.VIEWER)
    joined_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )


class RoomShareLink(SQLModel, table=True):
    id: UUID = Field(default_factory=uuid_mod.uuid4, primary_key=True)
    room_id: str = Field(foreign_key="room.id", index=True)
    token: str = Field(unique=True, index=True)
    access: ShareAccess = Field(default=ShareAccess.VIEW)
    created_by_id: UUID = Field(foreign_key="user.id")
    created_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )
    expires_at: datetime | None = Field(default=None, sa_type=UTCDateTime())
    revoked_at: datetime | None = Field(default=None, sa_type=UTCDateTime())


class Message(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    room_id: str = Field(foreign_key="room.id", index=True)
    user_id: UUID = Field(index=True)
    content: str
    created_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )
    updated_at: datetime | None = Field(default=None, sa_type=UTCDateTime())


class RoomGeometry(SQLModel, table=True):
    """Geometry instance in a room (particles, bonds, etc.)."""

    room_id: str = Field(foreign_key="room.id", primary_key=True)
    key: str = Field(primary_key=True)
    type: str  # Discriminator: "Sphere", "Bond", "Camera", etc.
    config: str  # Pydantic model JSON (BaseGeometry subclass), includes owner
    selection: str | None = None  # JSON list[int], NULL for cameras


class RoomBookmark(SQLModel, table=True):
    """Frame bookmark in a room."""

    room_id: str = Field(foreign_key="room.id", primary_key=True)
    frame_index: int = Field(primary_key=True)
    label: str


class SelectionGroup(SQLModel, table=True):
    """Named selection group in a room."""

    room_id: str = Field(foreign_key="room.id", primary_key=True)
    name: str = Field(primary_key=True)
    selections: str  # JSON dict[str, list[int]]


class RoomFigure(SQLModel, table=True):
    """Plotly figure in a room."""

    room_id: str = Field(foreign_key="room.id", primary_key=True)
    key: str = Field(primary_key=True)
    type: str = Field(default="plotly")
    data: str  # TEXT — Plotly JSON


class Screenshot(SQLModel, table=True):
    """Screenshot captured from a frontend session."""

    id: int | None = Field(default=None, primary_key=True)
    room_id: str = Field(
        sa_column=Column(String, ForeignKey("room.id", ondelete="CASCADE"), index=True)
    )
    format: str = Field(default="png")
    size: int = Field(default=0)
    width: int | None = None
    height: int | None = None
    status: str = Field(default="completed")
    created_by_id: UUID = Field(index=True)
    created_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )


class RoomPreset(SQLModel, table=True):
    """Visual preset stored per-room.

    Follows the same pattern as RoomGeometry: ``rules`` is a JSON-serialized
    string, consistent with how ``RoomGeometry.config`` stores geometry state.
    """

    room_id: str = Field(foreign_key="room.id", primary_key=True)
    name: str = Field(primary_key=True)
    description: str = ""
    rules: str  # JSON-serialized list[PresetRule]
    created_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )
    updated_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )


class ServerSettings(SQLModel, table=True):
    """Singleton table for server-wide configuration.

    Only one row with id=1 should exist. Use get_or_create pattern via
    `get_server_settings()` helper function.
    """

    id: int = Field(default=1, primary_key=True)
    default_room_id: str | None = Field(default=None, foreign_key="room.id")
