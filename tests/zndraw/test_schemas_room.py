"""Unit tests for room request/response schemas."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from zndraw.schemas import RoomCreate, RoomCreateResponse


def test_room_create_accepts_owner_and_name() -> None:
    body = RoomCreate(owner="alice-the-explorer", name="my-room")
    assert body.name == "my-room"
    assert body.owner == "alice-the-explorer"


def test_room_create_rejects_bad_name() -> None:
    with pytest.raises(ValidationError):
        RoomCreate(owner="alice-the-explorer", name="bad name with spaces")


def test_room_create_rejects_bad_owner() -> None:
    # UUID-shaped owner segments must hit the pattern gate.
    with pytest.raises(ValidationError):
        RoomCreate(owner="00000000-0000-0000-0000-000000000000", name="my-room")


def test_room_create_response_carries_composed_address() -> None:
    resp = RoomCreateResponse(
        room_id="alice-the-explorer/my-room", frame_count=0, created=True
    )
    assert resp.room_id == "alice-the-explorer/my-room"
