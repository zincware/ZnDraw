"""Tests for Socket.IO event Pydantic models."""

from datetime import UTC, datetime
from uuid import uuid4

import pytest
from pydantic import ValidationError

from zndraw.socket_events import (
    BookmarksInvalidate,
    FigureInvalidate,
    FramesInvalidate,
    FrameUpdate,
    GeometryInvalidate,
    Heartbeat,
    HeartbeatResponse,
    LockUpdate,
    MessageEdited,
    MessageNew,
    RoomJoin,
    RoomJoinResponse,
    RoomLeave,
    RoomLeaveResponse,
    RoomUpdate,
    SelectionGroupsInvalidate,
    SelectionInvalidate,
    SessionJoined,
    SessionLeft,
    Typing,
    TypingResponse,
    TypingStart,
    TypingStop,
    UserGet,
    UserGetResponse,
)

# =============================================================================
# socket_events: Request models
# =============================================================================


def test_room_join_serialization() -> None:
    """RoomJoin serializes with default client_type."""
    event = RoomJoin(room_id="room-1")
    data = event.model_dump(mode="json")
    assert data["room_id"] == "room-1"
    assert data["client_type"] == "frontend"


def test_room_join_custom_client_type() -> None:
    """RoomJoin accepts custom client_type."""
    event = RoomJoin(room_id="room-1", client_type="pyclient")
    assert event.client_type == "pyclient"


def test_room_join_missing_room_id() -> None:
    """RoomJoin requires room_id."""
    with pytest.raises(ValidationError):
        RoomJoin()  # type: ignore[call-arg]


@pytest.mark.parametrize(
    "model_cls",
    [RoomLeave, Heartbeat, TypingStart, TypingStop],
    ids=["RoomLeave", "Heartbeat", "TypingStart", "TypingStop"],
)
def test_room_id_only_request_roundtrip(model_cls: type) -> None:
    """Request models with only room_id serialize/deserialize correctly."""
    instance = model_cls(room_id="test-room")
    data = instance.model_dump(mode="json")
    restored = model_cls.model_validate(data)
    assert restored.room_id == "test-room"


def test_user_get_no_fields() -> None:
    """UserGet has no required fields."""
    event = UserGet()
    assert event.model_dump(mode="json") == {}


# =============================================================================
# socket_events: Response models
# =============================================================================


def test_room_join_response_serialization() -> None:
    """RoomJoinResponse serializes all fields."""
    resp = RoomJoinResponse(
        room_id="room-1",
        session_id="sess-1",
        step=3,
        frame_count=10,
        locked=True,
        camera_key="cam:user@test.com:a1b2c3d4",
    )
    data = resp.model_dump(mode="json")
    assert data == {
        "room_id": "room-1",
        "session_id": "sess-1",
        "step": 3,
        "frame_count": 10,
        "locked": True,
        "camera_key": "cam:user@test.com:a1b2c3d4",
        "progress_trackers": {},
    }


def test_room_join_response_missing_required() -> None:
    """RoomJoinResponse requires all fields (no defaults)."""
    with pytest.raises(ValidationError):
        RoomJoinResponse(room_id="room-1")  # type: ignore[call-arg]


def test_room_leave_response_roundtrip() -> None:
    """RoomLeaveResponse round-trips correctly."""
    resp = RoomLeaveResponse(room_id="room-1")
    assert RoomLeaveResponse.model_validate(resp.model_dump()).room_id == "room-1"


@pytest.mark.parametrize(
    "model_cls", [HeartbeatResponse, TypingResponse], ids=["Heartbeat", "Typing"]
)
def test_ok_response_default(model_cls: type) -> None:
    """HeartbeatResponse and TypingResponse default to status='ok'."""
    resp = model_cls()
    assert resp.model_dump(mode="json") == {"status": "ok"}


def test_user_get_response_serialization() -> None:
    """UserGetResponse serializes UUID and fields."""
    uid = uuid4()
    resp = UserGetResponse(id=uid, email="user@example.com", is_superuser=False)
    data = resp.model_dump(mode="json")
    assert data["id"] == str(uid)
    assert data["email"] == "user@example.com"
    assert data["is_superuser"] is False


def test_user_get_response_deserializes_uuid_string() -> None:
    """UserGetResponse deserializes UUID from string."""
    uid = uuid4()
    resp = UserGetResponse.model_validate(
        {"id": str(uid), "email": "a@b.com", "is_superuser": True}
    )
    assert resp.id == uid


# =============================================================================
# socket_events: Broadcast models - Frame
# =============================================================================


def test_frame_update_data() -> None:
    """FrameUpdate serializes room_id and frame."""
    broadcast = FrameUpdate(room_id="test-room-uuid", frame=5)
    data = broadcast.model_dump(mode="json")
    assert data["room_id"] == "test-room-uuid"
    assert data["frame"] == 5


def test_frames_invalidate_add_action() -> None:
    """FramesInvalidate serializes add action with count."""
    broadcast = FramesInvalidate(
        room_id="test-room-uuid",
        action="add",
        count=10,
    )
    data = broadcast.model_dump(mode="json")
    assert data["action"] == "add"
    assert data["count"] == 10
    assert data["indices"] is None


def test_frames_invalidate_delete_action() -> None:
    """FramesInvalidate serializes delete action with indices."""
    broadcast = FramesInvalidate(
        room_id="test-room-uuid",
        action="delete",
        indices=[1, 2, 3],
    )
    data = broadcast.model_dump(mode="json")
    assert data["action"] == "delete"
    assert data["indices"] == [1, 2, 3]


def test_frames_invalidate_clear_action() -> None:
    """FramesInvalidate serializes clear action."""
    broadcast = FramesInvalidate(
        room_id="test-room-uuid",
        action="clear",
        count=0,
    )
    data = broadcast.model_dump(mode="json")
    assert data["action"] == "clear"
    assert data["count"] == 0


def test_frames_invalidate_invalid_action() -> None:
    """FramesInvalidate rejects invalid action literals."""
    with pytest.raises(ValidationError):
        FramesInvalidate(room_id="r", action="invalid")  # type: ignore[arg-type]


# =============================================================================
# socket_events: Broadcast models - Session
# =============================================================================


def test_session_joined_fields() -> None:
    """SessionJoined carries user identity fields."""
    uid = uuid4()
    event = SessionJoined(
        room_id="room-1", user_id=uid, sid="socket-123", email="user@x.com"
    )
    data = event.model_dump(mode="json")
    assert data["room_id"] == "room-1"
    assert data["user_id"] == str(uid)
    assert data["sid"] == "socket-123"
    assert data["email"] == "user@x.com"


def test_session_joined_email_optional() -> None:
    """SessionJoined defaults email to None."""
    uid = uuid4()
    event = SessionJoined(room_id="room-1", user_id=uid, sid="s1")
    assert event.email is None


def test_session_left_fields() -> None:
    """SessionLeft carries minimal identity fields."""
    uid = uuid4()
    event = SessionLeft(room_id="room-1", user_id=uid, sid="s1")
    data = event.model_dump(mode="json")
    assert data["room_id"] == "room-1"
    assert data["user_id"] == str(uid)
    assert data["sid"] == "s1"


# =============================================================================
# socket_events: Broadcast models - Other invalidations
# =============================================================================


@pytest.mark.parametrize(
    "model_cls",
    [SelectionInvalidate, SelectionGroupsInvalidate],
    ids=["Selection", "SelectionGroups"],
)
def test_room_only_invalidation(model_cls: type) -> None:
    """Invalidation models with only room_id."""
    event = model_cls(room_id="room-1")
    assert event.model_dump(mode="json") == {"room_id": "room-1"}


def test_geometry_invalidate_fields() -> None:
    """GeometryInvalidate carries operation and key."""
    event = GeometryInvalidate(room_id="r", operation="set", key="mol1")
    data = event.model_dump(mode="json")
    assert data["operation"] == "set"
    assert data["key"] == "mol1"


def test_bookmarks_invalidate_fields() -> None:
    """BookmarksInvalidate carries index and operation."""
    event = BookmarksInvalidate(room_id="r", index=5, operation="delete")
    data = event.model_dump(mode="json")
    assert data["index"] == 5
    assert data["operation"] == "delete"


def test_figure_invalidate_fields() -> None:
    """FigureInvalidate carries key and operation."""
    event = FigureInvalidate(room_id="r", key="fig-1", operation="set")
    data = event.model_dump(mode="json")
    assert data["key"] == "fig-1"
    assert data["operation"] == "set"


def test_lock_update_defaults() -> None:
    """LockUpdate optional fields default to None."""
    event = LockUpdate(room_id="r", action="acquired")
    assert event.user_id is None
    assert event.msg is None


def test_lock_update_full() -> None:
    """LockUpdate serializes all fields."""
    event = LockUpdate(
        room_id="r",
        action="released",
        user_id="user-123",
        msg="editing geometries",
    )
    data = event.model_dump(mode="json")
    assert data["action"] == "released"
    assert data["user_id"] == "user-123"
    assert data["msg"] == "editing geometries"


def test_room_update_snapshot() -> None:
    """RoomUpdate is a full room snapshot with all required fields."""
    event = RoomUpdate(id="r", frame_count=5, locked=True, is_default=False)
    assert event.id == "r"
    assert event.frame_count == 5
    assert event.locked is True
    assert event.is_default is False
    assert event.description is None


# =============================================================================
# socket_events: Broadcast models - Messages and Typing
# =============================================================================


def test_message_new_fields() -> None:
    """MessageNew serializes all fields including UUID."""
    uid = uuid4()
    now = datetime.now(tz=UTC)
    event = MessageNew(
        id=1, room_id="r", user_id=uid, content="hello", created_at=now, email="a@b.c"
    )
    data = event.model_dump(mode="json")
    assert data["id"] == 1
    assert data["user_id"] == str(uid)
    assert data["content"] == "hello"
    assert data["email"] == "a@b.c"
    assert data["updated_at"] is None


def test_message_edited_fields() -> None:
    """MessageEdited carries id, room_id, content, updated_at."""
    now = datetime.now(tz=UTC)
    event = MessageEdited(id=1, room_id="r", content="edited", updated_at=now)
    data = event.model_dump(mode="json")
    assert data["id"] == 1
    assert data["content"] == "edited"


def test_typing_broadcast_fields() -> None:
    """Typing broadcast carries user identity and typing state."""
    uid = uuid4()
    event = Typing(room_id="r", user_id=uid, is_typing=True, email="a@b.c")
    data = event.model_dump(mode="json")
    assert data["is_typing"] is True
    assert data["email"] == "a@b.c"


# =============================================================================
# Deserialization roundtrips
# =============================================================================


@pytest.mark.parametrize(
    ("model_cls", "kwargs"),
    [
        (FrameUpdate, {"room_id": "r", "frame": 42}),
        (FramesInvalidate, {"room_id": "r", "action": "modify", "indices": [0, 1]}),
        (
            SessionJoined,
            {
                "room_id": "r",
                "user_id": "00000000-0000-0000-0000-000000000001",
                "sid": "s",
            },
        ),
        (
            SessionLeft,
            {
                "room_id": "r",
                "user_id": "00000000-0000-0000-0000-000000000001",
                "sid": "s",
            },
        ),
        (GeometryInvalidate, {"room_id": "r", "operation": "delete", "key": "k"}),
    ],
    ids=[
        "FrameUpdate",
        "FramesInvalidate",
        "SessionJoined",
        "SessionLeft",
        "GeometryInvalidate",
    ],
)
def test_model_roundtrip(model_cls: type, kwargs: dict) -> None:
    """Models survive a dump -> validate roundtrip."""
    instance = model_cls(**kwargs)
    data = instance.model_dump(mode="json")
    restored = model_cls.model_validate(data)
    assert restored == instance
