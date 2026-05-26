"""Regression: RoomRenamed reaches the previous owner's user channel (finding #14).

Alice owns room "foo". A Group G has member Bob. A third client is "in the room"
on the surrogate channel. Alice (as superuser) transfers foo from her namespace
to G's namespace. Expectations:

  - Alice receives RoomRenamed on user:<alice>
  - Bob receives RoomUpdate on user:<bob>     (existing broadcast_room_update path)
  - The in-room channel room:<surrogate> receives RoomRenamed
"""

from __future__ import annotations

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.access import GroupRole, Visibility
from zndraw.models import Group, GroupMembership


@pytest.mark.asyncio
async def test_room_transfer_notifies_previous_owner_user_channel(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    alice, alice_token = await create_test_user_in_db(
        session, email="alice@local.test", is_superuser=True
    )
    bob, _bob_token = await create_test_user_in_db(session, email="bob@local.test")
    room = await create_test_room(session, alice, room_name="foo")

    group = Group(name="G", created_by_id=alice.id)
    session.add(group)
    await session.flush()
    session.add(
        GroupMembership(user_id=alice.id, group_id=group.id, role=GroupRole.ADMIN)
    )
    session.add(
        GroupMembership(user_id=bob.id, group_id=group.id, role=GroupRole.MEMBER)
    )
    await session.commit()

    response = await client.patch(
        f"/v1/rooms/{room.public_address}",
        json={"new_owner_id": str(group.id), "visibility": Visibility.GROUP.value},
        headers=auth_header(alice_token),
    )
    assert response.status_code == 200, response.text

    renamed = [e for e in mock_sio.emitted if e["event"] == "room_renamed"]
    assert renamed, "room_renamed not emitted"

    rooms = {e["room"] for e in renamed}
    assert f"user:{alice.id}" in rooms, (
        f"RoomRenamed did not reach previous owner channel; saw {rooms}"
    )
    assert f"room:{room.id}" in rooms, (
        f"RoomRenamed did not reach in-room channel; saw {rooms}"
    )

    room_updates = [e for e in mock_sio.emitted if e["event"] == "room_update"]
    update_targets = {e["room"] for e in room_updates}
    assert f"user:{bob.id}" in update_targets, (
        "RoomUpdate did not fan out to the new owner's group members; "
        f"saw {update_targets}"
    )

    payload = renamed[0]["data"]
    assert payload["old_address"] == f"{alice.id}/foo"
    assert payload["room_address"] == f"{group.id}/foo"
    # MockSioServer captures via model_dump(); UUID fields surface as UUID instances.
    assert str(payload["room_id"]) == room.id
