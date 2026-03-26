"""Tests for frame selection REST API endpoints."""

import json

import pytest
from helpers import MockSioServer, auth_header, create_test_room, create_test_user_in_db
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession



# =============================================================================
# GET /v1/rooms/{room_id}/frame-selection
# =============================================================================


@pytest.mark.asyncio
async def test_get_returns_null_when_empty(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET returns null frameSelection for a new room."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=auth_header(token)
    )
    assert response.status_code == 200
    assert response.json()["frame_selection"] is None


@pytest.mark.asyncio
async def test_get_returns_stored_indices(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET returns stored indices when frame_selection column is set."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Set column directly
    room.frame_selection = json.dumps([2, 5, 10])
    session.add(room)
    await session.commit()

    response = await client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=auth_header(token)
    )
    assert response.status_code == 200
    assert response.json()["frame_selection"] == [2, 5, 10]


@pytest.mark.asyncio
async def test_get_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/nonexistent/frame-selection", headers=auth_header(token)
    )
    assert response.status_code == 404


# =============================================================================
# PUT /v1/rooms/{room_id}/frame-selection
# =============================================================================


@pytest.mark.asyncio
async def test_put_stores_indices(
    client: AsyncClient, session: AsyncSession
) -> None:
    """PUT stores indices and GET returns them."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    put_resp = await client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [1, 3, 7]},
        headers=auth_header(token),
    )
    assert put_resp.status_code == 200
    assert put_resp.json()["success"] is True

    get_resp = await client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=auth_header(token)
    )
    assert get_resp.json()["frame_selection"] == [1, 3, 7]


@pytest.mark.asyncio
async def test_put_broadcasts_socket_event(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """PUT broadcasts FrameSelectionUpdate socket event."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    mock_sio.emitted.clear()

    await client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [0, 4]},
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) >= 1
    evt = mock_sio.emitted[-1]
    assert evt["event"] == "frame_selection_update"
    assert evt["data"]["indices"] == [0, 4]


@pytest.mark.asyncio
async def test_put_empty_list_clears_selection(
    client: AsyncClient, session: AsyncSession
) -> None:
    """PUT with empty list clears selection (GET returns null)."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Set some indices first
    await client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [1, 2]},
        headers=auth_header(token),
    )

    # Clear with empty list
    await client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": []},
        headers=auth_header(token),
    )

    get_resp = await client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=auth_header(token)
    )
    assert get_resp.json()["frame_selection"] is None


@pytest.mark.asyncio
async def test_put_rejects_negative_indices(
    client: AsyncClient, session: AsyncSession
) -> None:
    """PUT with negative indices returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": [1, -2, 3]},
        headers=auth_header(token),
    )
    assert response.status_code == 400


@pytest.mark.asyncio
async def test_put_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """PUT for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.put(
        "/v1/rooms/nonexistent/frame-selection",
        json={"indices": [0]},
        headers=auth_header(token),
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_roundtrip(client: AsyncClient, session: AsyncSession) -> None:
    """PUT then GET returns consistent data."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    indices = [0, 2, 4, 6, 8]

    await client.put(
        f"/v1/rooms/{room.id}/frame-selection",
        json={"indices": indices},
        headers=auth_header(token),
    )

    get_resp = await client.get(
        f"/v1/rooms/{room.id}/frame-selection", headers=auth_header(token)
    )
    assert get_resp.status_code == 200
    assert get_resp.json()["frame_selection"] == indices
