"""Unit tests for broadcast_room_update helper routing logic.

These tests cover cases that cannot be driven via REST today:
private rooms (no REST way to create them) and the exact set of
channels a broadcast targets. Real socket integration is covered by
test_socketio_rooms.py.
"""

from uuid import uuid4

import pytest
from helpers import MockSioServer
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.models import MemberRole, Room, RoomMembership
from zndraw.routes.rooms import broadcast_room_update
from zndraw.storage import FrameStorage


@pytest.mark.asyncio
async def test_broadcast_public_room_targets_feed(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Public rooms broadcast to the shared rooms:feed channel only."""
    room = Room(id="pub", is_public=True)
    session.add(room)
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    rooms_targeted = [call["room"] for call in sio.emitted]
    assert rooms_targeted == ["rooms:feed"]


@pytest.mark.asyncio
async def test_broadcast_private_room_targets_each_member(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Private rooms broadcast to each member's user:{uid} channel and
    nothing else — never rooms:feed, never non-members."""
    member_a = uuid4()
    member_b = uuid4()
    non_member = uuid4()  # noqa: F841 — intentionally unreferenced

    room = Room(id="priv", is_public=False)
    session.add(room)
    session.add(
        RoomMembership(room_id="priv", user_id=member_a, role=MemberRole.OWNER)
    )
    session.add(
        RoomMembership(room_id="priv", user_id=member_b, role=MemberRole.MEMBER)
    )
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    rooms_targeted = sorted(call["room"] for call in sio.emitted)
    assert rooms_targeted == sorted(
        [f"user:{member_a}", f"user:{member_b}"]
    )
    assert "rooms:feed" not in rooms_targeted


@pytest.mark.asyncio
async def test_broadcast_private_room_with_no_members_is_noop(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Private room with no memberships emits nothing — no channel is
    authorized, so no broadcast occurs."""
    room = Room(id="orphan", is_public=False)
    session.add(room)
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    assert sio.emitted == []
