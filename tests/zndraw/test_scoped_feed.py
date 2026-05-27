"""Room-update feed delivery is scoped to authorized viewers.

These tests drive `broadcast_room_update` directly (no socket client), which
keeps the assertions isolated to the routing logic — the MockSioServer
captures emit calls and the tests check which rooms received the event.
"""

from uuid import uuid4

import pytest
from helpers import MockSioServer
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import GroupRole, Visibility
from zndraw.models import Group, GroupMembership, Room
from zndraw.routes.rooms import broadcast_room_update
from zndraw.storage import FrameStorage


@pytest.mark.asyncio
async def test_public_room_broadcasts_to_rooms_feed(
    session: AsyncSession, frame_storage: FrameStorage
) -> None:
    sio = MockSioServer()
    room = Room(
        room_name="pub-fan",
        owner_user_id=uuid4(),
        created_by_id=uuid4(),
        visibility=Visibility.PUBLIC,
    )
    session.add(room)
    await session.commit()

    await broadcast_room_update(sio, session, frame_storage, room)

    targets = {call["room"] for call in sio.emitted}
    assert targets == {"rooms:feed"}


@pytest.mark.asyncio
async def test_private_room_emits_only_to_owner(
    session: AsyncSession, frame_storage: FrameStorage
) -> None:
    sio = MockSioServer()
    owner_id = uuid4()
    room = Room(
        room_name="prv-fan",
        owner_user_id=owner_id,
        created_by_id=owner_id,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    await session.commit()

    await broadcast_room_update(sio, session, frame_storage, room)

    targets = {call["room"] for call in sio.emitted}
    assert targets == {f"user:{owner_id}"}


@pytest.mark.asyncio
async def test_group_room_emits_to_each_member(
    session: AsyncSession, frame_storage: FrameStorage
) -> None:
    sio = MockSioServer()
    admin_id = uuid4()
    member_id = uuid4()
    outsider_id = uuid4()  # NOT added to group — must NOT receive

    group = Group(name=f"fan-test-{uuid4().hex[:6]}", created_by_id=admin_id)
    session.add(group)
    await session.commit()
    session.add_all(
        [
            GroupMembership(group_id=group.id, user_id=admin_id, role=GroupRole.ADMIN),
            GroupMembership(
                group_id=group.id, user_id=member_id, role=GroupRole.VIEWER
            ),
        ]
    )
    await session.commit()

    room = Room(
        room_name="grp-fan",
        owner_group_id=group.id,
        created_by_id=admin_id,
        visibility=Visibility.GROUP,
    )
    session.add(room)
    await session.commit()

    await broadcast_room_update(sio, session, frame_storage, room)

    targets = {call["room"] for call in sio.emitted}
    assert f"user:{admin_id}" in targets
    assert f"user:{member_id}" in targets
    assert f"user:{outsider_id}" not in targets
    assert "rooms:feed" not in targets
