"""Regression: FrameSelectionUpdate carries room_address and routes via broadcast_to_room (review #1)."""

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession


@pytest.mark.asyncio
async def test_frame_selection_update_carries_room_address(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.public_address}/frame-selection",
        json={"indices": [0, 1, 2]},
        headers=auth_header(token),
    )
    assert response.status_code == 200, response.text

    emits = [e for e in mock_sio.emitted if e["event"] == "frame_selection_update"]
    assert len(emits) == 1, f"expected 1 frame_selection_update, got {len(emits)}"
    captured = emits[0]
    data = captured["data"]
    assert "room_id" in data, (
        f"FrameSelectionUpdate missing room_id; got keys: {list(data)}"
    )
    assert "room_address" in data, (
        f"FrameSelectionUpdate missing room_address; got keys: {list(data)}"
    )
    assert str(data["room_id"]) == room.id
    assert data["room_address"] == room.public_address
    assert captured["room"] == f"room:{room.id}"
