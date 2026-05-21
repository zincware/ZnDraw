"""Unit tests for broadcast_room_update helper routing logic.

These tests cover the scoped delivery introduced in Task 18:
- PUBLIC rooms broadcast to the shared rooms:feed channel.
- GROUP rooms fan out to each group member's user:{uid} channel.
- PRIVATE rooms deliver only to the owner's user:{uid} channel.

Real socket integration is covered by test_socketio_rooms.py.
"""

from uuid import uuid4

import pytest
from helpers import MockSioServer
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import Visibility
from zndraw.models import Group, GroupMembership, Room
from zndraw.routes.rooms import broadcast_room_update
from zndraw.storage import FrameStorage


@pytest.mark.asyncio
async def test_broadcast_public_room_targets_feed(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Public rooms broadcast to the shared rooms:feed channel only."""
    owner_id = uuid4()
    room = Room(
        room_name="pub",
        visibility=Visibility.PUBLIC,
        owner_user_id=owner_id,
        created_by_id=owner_id,
    )
    session.add(room)
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    rooms_targeted = [call["room"] for call in sio.emitted]
    assert rooms_targeted == ["rooms:feed"]


@pytest.mark.asyncio
async def test_broadcast_private_room_targets_owner(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Private rooms broadcast only to the owner's user:{uid} channel."""
    owner_id = uuid4()
    non_owner = uuid4()  # noqa: F841 — intentionally unreferenced

    room = Room(
        room_name="priv",
        visibility=Visibility.PRIVATE,
        owner_user_id=owner_id,
        created_by_id=owner_id,
    )
    session.add(room)
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    rooms_targeted = [call["room"] for call in sio.emitted]
    assert rooms_targeted == [f"user:{owner_id}"]
    assert "rooms:feed" not in rooms_targeted


@pytest.mark.asyncio
async def test_broadcast_group_room_targets_each_member(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Group rooms broadcast to each group-member's user:{uid} channel."""
    member_a = uuid4()
    member_b = uuid4()
    non_member = uuid4()  # noqa: F841 — intentionally unreferenced
    group_id = uuid4()
    creator_id = uuid4()

    group = Group(id=group_id, name="test-group", created_by_id=creator_id)
    session.add(group)
    session.add(GroupMembership(group_id=group_id, user_id=member_a))
    session.add(GroupMembership(group_id=group_id, user_id=member_b))
    room = Room(
        room_name="grp",
        visibility=Visibility.GROUP,
        owner_group_id=group_id,
        created_by_id=creator_id,
    )
    session.add(room)
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    rooms_targeted = sorted(call["room"] for call in sio.emitted)
    assert rooms_targeted == sorted([f"user:{member_a}", f"user:{member_b}"])
    assert "rooms:feed" not in rooms_targeted
