import uuid as uuid_mod
from datetime import UTC, datetime
from enum import Enum
from uuid import UUID

from sqlalchemy import Column, ForeignKey, String, TypeDecorator
from sqlalchemy.types import DateTime
from sqlmodel import Field, SQLModel
from zndraw_joblib.models import Job, Task, Worker, WorkerJobLink  # noqa: F401


class UTCDateTime(TypeDecorator):
    """SQLAlchemy type that ensures datetimes are always UTC-aware.

    SQLite strips timezone info on storage. This type decorator
    re-attaches UTC on load so consumers never see naive datetimes.
    """

    impl = DateTime
    cache_ok = True

    def process_result_value(
        self, value: datetime | None, dialect: object
    ) -> datetime | None:
        if value is not None and value.tzinfo is None:
            return value.replace(tzinfo=UTC)
        return value


class MemberRole(str, Enum):
    MEMBER = "member"
    MODERATOR = "moderator"
    OWNER = "owner"


class Room(SQLModel, table=True):
    """Room model with string UUID as primary key.

    This matches the frontend expectation of string room IDs (UUIDs).
    Frames are stored separately in the storage layer (memory/LMDB/etc).
    """

    id: str = Field(default_factory=lambda: str(uuid_mod.uuid4()), primary_key=True)
    description: str | None = None
    created_by_id: UUID | None = Field(default=None, index=True)
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    is_public: bool = Field(default=True)
    locked: bool = Field(default=False)  # Admin lock status
    step: int = Field(default=0)
    frame_selection: str | None = Field(default=None)  # JSON list[int]
    default_camera: str | None = Field(default=None)  # Geometry key for default camera


class Message(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    room_id: str = Field(foreign_key="room.id", index=True)
    user_id: UUID = Field(index=True)
    content: str
    created_at: datetime = Field(
        default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime()
    )
    updated_at: datetime | None = Field(default=None, sa_type=UTCDateTime())


class RoomMembership(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    room_id: str = Field(foreign_key="room.id", index=True)
    user_id: UUID = Field(index=True)
    joined_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    role: MemberRole = Field(default=MemberRole.MEMBER)


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


class ServerSettings(SQLModel, table=True):
    """Singleton table for server-wide configuration.

    Only one row with id=1 should exist. Use get_or_create pattern via
    `get_server_settings()` helper function.
    """

    id: int = Field(default=1, primary_key=True)
    default_room_id: str | None = Field(default=None, foreign_key="room.id")
