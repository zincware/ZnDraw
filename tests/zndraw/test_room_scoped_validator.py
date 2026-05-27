"""Regression: validator rejects NIL room_id + composed room_address (review #2)."""

from uuid import UUID

import pytest
from pydantic import ValidationError

from zndraw.socket_events import FramesInvalidate


def test_nil_room_id_with_composed_address_rejected() -> None:
    """Direct construction with NIL room_id but a composed address is a footgun."""
    composed = "11111111-1111-1111-1111-111111111111/demo"
    with pytest.raises(ValidationError, match="room_id"):
        FramesInvalidate(
            room_id=UUID(int=0),
            room_address=composed,
            action="add",
        )


def test_nil_room_id_with_sigil_address_allowed() -> None:
    """Sigil addresses (``@global``/``@internal``) may pair with NIL room_id."""
    event = FramesInvalidate(
        room_id=UUID(int=0),
        room_address="@global",
        action="clear",
    )
    assert event.room_id == UUID(int=0)
    assert event.room_address == "@global"


def test_real_room_id_with_composed_address_allowed() -> None:
    real = UUID("22222222-2222-2222-2222-222222222222")
    composed = f"{real}/demo"
    event = FramesInvalidate(
        room_id=real,
        room_address=composed,
        action="modify",
        indices=[0],
    )
    assert event.room_id == real
