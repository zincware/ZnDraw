"""Tests for Selection Groups REST API endpoints."""

import json

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.exceptions import SelectionGroupNotFound
from zndraw.models import SelectionGroup
from zndraw.schemas import StatusResponse


async def _add_selection_group(
    session: AsyncSession,
    room_id: str,
    name: str,
    selections: dict[str, list[int]],
) -> None:
    """Add a selection group to the database."""
    session.add(
        SelectionGroup(
            room_id=room_id,
            name=name,
            selections=json.dumps(selections),
        )
    )
    await session.commit()


# =============================================================================
# List Selection Groups Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_selection_groups_returns_empty(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns empty groups for new room."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/selection-groups",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["items"] == {}


@pytest.mark.asyncio
async def test_list_selection_groups_returns_all(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns all stored groups."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    group_a = {"sphere": [1, 2], "cube": [3]}
    group_b = {"sphere": [4, 5]}
    await _add_selection_group(session, room.id, "group_a", group_a)
    await _add_selection_group(session, room.id, "group_b", group_b)

    response = await client.get(
        f"/v1/rooms/{room.id}/selection-groups",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    data = response.json()
    assert data["items"]["group_a"] == group_a
    assert data["items"]["group_b"] == group_b


# =============================================================================
# Get Selection Group Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_nonexistent_selection_group_returns_404(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET group returns 404 when not exists."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == SelectionGroupNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_selection_group_returns_stored(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET group returns stored selections."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    group_data = {"sphere": [1, 2], "cube": [3, 4]}
    await _add_selection_group(session, room.id, "mygroup", group_data)

    response = await client.get(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["group"] == group_data


# =============================================================================
# Update Selection Group Tests
# =============================================================================


@pytest.mark.asyncio
async def test_update_selection_group_stores_data(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT group stores selection data."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    group_data = {"sphere": [0, 1], "cube": [2, 3]}
    response = await client.put(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        json={"selections": group_data},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["status"] == "ok"

    # Verify stored in DB
    row = await session.get(SelectionGroup, (room.id, "mygroup"))
    assert row is not None
    assert json.loads(row.selections) == group_data


@pytest.mark.asyncio
async def test_update_selection_group_broadcasts(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """Test PUT group broadcasts selection_groups:invalidate event."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await client.put(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        json={"selections": {"sphere": [0]}},
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "selection_groups_invalidate"


# =============================================================================
# Delete Selection Group Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_nonexistent_selection_group_returns_404(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE on nonexistent selection group returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(
        f"/v1/rooms/{room.id}/selection-groups/nonexistent",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == SelectionGroupNotFound.type_uri()


@pytest.mark.asyncio
async def test_delete_selection_group(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE removes group."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_selection_group(session, room.id, "mygroup", {})

    response = await client.delete(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await session.get(SelectionGroup, (room.id, "mygroup"))
    assert row is None


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_selection_groups_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(f"/v1/rooms/{room.id}/selection-groups")
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_selection_groups_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/selection-groups",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
