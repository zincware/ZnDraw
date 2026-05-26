"""Regression: RoomRenamed fans out to every previous group member (review #3)."""

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.access import GroupRole, Visibility
from zndraw.models import Group, GroupMembership, Room


@pytest.mark.asyncio
async def test_room_renamed_fans_out_to_previous_group_members(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    alice, alice_token = await create_test_user_in_db(
        session, email="alice@local.test", is_superuser=True
    )
    bob, _ = await create_test_user_in_db(session, email="bob@local.test")
    carol, _ = await create_test_user_in_db(session, email="carol@local.test")

    group = Group(name="G", created_by_id=alice.id)
    session.add(group)
    await session.flush()
    session.add(
        GroupMembership(user_id=alice.id, group_id=group.id, role=GroupRole.ADMIN)
    )
    session.add(
        GroupMembership(user_id=bob.id, group_id=group.id, role=GroupRole.MEMBER)
    )
    session.add(
        GroupMembership(user_id=carol.id, group_id=group.id, role=GroupRole.MEMBER)
    )
    await session.commit()

    room = Room(
        room_name="foo",
        owner_group_id=group.id,
        created_by_id=alice.id,
        visibility=Visibility.GROUP,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    response = await client.patch(
        f"/v1/rooms/{room.public_address}",
        json={"new_owner_id": str(alice.id), "visibility": Visibility.PRIVATE.value},
        headers=auth_header(alice_token),
    )
    assert response.status_code == 200, response.text

    renamed = [e for e in mock_sio.emitted if e["event"] == "room_renamed"]
    rooms_seen = {e["room"] for e in renamed}
    for member in (alice, bob, carol):
        assert f"user:{member.id}" in rooms_seen, (
            f"RoomRenamed missing for previous group member {member.email}; "
            f"saw {rooms_seen}"
        )
