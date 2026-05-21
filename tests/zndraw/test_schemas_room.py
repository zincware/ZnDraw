"""Unit tests for room request/response schemas."""

from __future__ import annotations

from uuid import uuid4

import pytest
from pydantic import ValidationError

from zndraw.schemas import RoomCreate, RoomCreateResponse


def test_room_create_accepts_owner_id_and_name() -> None:
    body = RoomCreate(owner_id=uuid4(), name="my-room")
    assert body.name == "my-room"


def test_room_create_rejects_bad_name() -> None:
    with pytest.raises(ValidationError):
        RoomCreate(owner_id=uuid4(), name="bad name with spaces")


def test_room_create_response_carries_composed_address() -> None:
    owner = uuid4()
    resp = RoomCreateResponse(
        room_id=f"{owner}/my-room", frame_count=0, created=True
    )
    assert resp.room_id == f"{owner}/my-room"
