"""Layer 1: static schema audit for RoomScopedEvent subclasses.

Fails fast at collect time if any subclass forgets to declare both
``room_id: UUID`` and ``room_address: str``.
"""

from __future__ import annotations

import enum
import typing
from datetime import UTC, datetime
from uuid import UUID

import pytest

# Triggers import of every module that defines RoomScopedEvent subclasses.
import zndraw.socket_events  # noqa: F401
import zndraw_joblib.events  # noqa: F401
from zndraw.socket_events import RoomScopedEvent


def _all_subclasses(cls: type) -> set[type]:
    out: set[type] = set()
    stack: list[type] = list(cls.__subclasses__())
    while stack:
        sub = stack.pop()
        if sub in out:
            continue
        out.add(sub)
        stack.extend(sub.__subclasses__())
    return out


def _sample_value(annotation: object) -> object:
    """Best-effort default for required non-id fields used only by the audit test."""
    if annotation is int:
        return 0
    if annotation is str:
        return "x"
    if annotation is bool:
        return False
    if annotation is float:
        return 0.0
    if annotation is UUID:
        return UUID("22222222-2222-2222-2222-222222222222")
    if annotation is datetime:
        return datetime.now(UTC)
    origin = typing.get_origin(annotation)
    if origin is typing.Literal:
        return typing.get_args(annotation)[0]
    if origin is list:
        return []
    if isinstance(annotation, type) and issubclass(annotation, enum.Enum):
        return next(iter(annotation))
    # Fallback — most enums/literals/optionals validate from a short string.
    return "x"


@pytest.mark.protected
def test_all_room_scoped_events_have_both_ids() -> None:
    subclasses = _all_subclasses(RoomScopedEvent)
    assert subclasses, "RoomScopedEvent has no subclasses — import surface broken"
    for cls in subclasses:
        fields = cls.model_fields
        assert "room_id" in fields, f"{cls.__name__} missing room_id"
        assert "room_address" in fields, f"{cls.__name__} missing room_address"
        assert fields["room_id"].annotation is UUID, (
            f"{cls.__name__}.room_id must be UUID, got {fields['room_id'].annotation}"
        )
        assert fields["room_address"].annotation is str, (
            f"{cls.__name__}.room_address must be str, "
            f"got {fields['room_address'].annotation}"
        )


@pytest.mark.protected
def test_for_room_round_trip_uses_room_attributes() -> None:
    """``.for_room(room, **kw)`` must populate both ids from the Room instance."""

    class _RoomStub:
        id = "00000000-0000-0000-0000-000000000001"
        public_address = "11111111-1111-1111-1111-111111111111/demo"

    for cls in _all_subclasses(RoomScopedEvent):
        extra: dict[str, object] = {}
        for field_name, field in cls.model_fields.items():
            if field_name in {"room_id", "room_address"}:
                continue
            if not field.is_required():
                continue
            extra[field_name] = _sample_value(field.annotation)
        event = cls.for_room(_RoomStub(), **extra)
        assert str(event.room_id) == _RoomStub.id
        assert event.room_address == _RoomStub.public_address
