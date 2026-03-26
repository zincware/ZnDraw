"""Tests for Figures REST API endpoints."""

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.exceptions import FigureNotFound
from zndraw.models import RoomFigure
from zndraw.schemas import FigureData, StatusResponse


async def _add_figure(
    session: AsyncSession, room_id: str, key: str, data: str, fig_type: str = "plotly"
) -> None:
    """Add a figure directly to the database."""
    session.add(RoomFigure(room_id=room_id, key=key, type=fig_type, data=data))
    await session.commit()


# =============================================================================
# List Figures Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_figures_returns_empty_initially(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns empty figures list for new room."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/figures",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["items"] == []


@pytest.mark.asyncio
async def test_list_figures_returns_all_keys(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns all figure keys."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_figure(session, room.id, "chart1", '{"data": []}')
    await _add_figure(session, room.id, "chart2", '{"data": []}')

    response = await client.get(
        f"/v1/rooms/{room.id}/figures",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert set(data["items"]) == {"chart1", "chart2"}


# =============================================================================
# Get Figure Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_figure_returns_data(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET single figure returns the figure data."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_figure(
        session, room.id, "my_chart", '{"data": [1, 2, 3], "layout": {}}'
    )

    response = await client.get(
        f"/v1/rooms/{room.id}/figures/my_chart",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["key"] == "my_chart"
    assert data["figure"]["type"] == "plotly"
    assert data["figure"]["data"] == '{"data": [1, 2, 3], "layout": {}}'


@pytest.mark.asyncio
async def test_get_figure_returns_404_for_nonexistent(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns 404 for nonexistent figure."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/figures/nonexistent",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "figure-not-found" in response.json()["type"]


# =============================================================================
# Create/Update Figure Tests
# =============================================================================


@pytest.mark.asyncio
async def test_create_figure_stores_data(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test POST creates a new figure."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    figure_data = FigureData(
        type="plotly", data='{"data": [], "layout": {"title": "Test"}}'
    )
    response = await client.post(
        f"/v1/rooms/{room.id}/figures/new_chart",
        json={"figure": figure_data.model_dump()},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    body = response.json()
    assert body["key"] == "new_chart"
    assert body["created"] is True

    # Verify persisted in DB
    row = await session.get(RoomFigure, (room.id, "new_chart"))
    assert row is not None
    assert row.type == "plotly"


@pytest.mark.asyncio
async def test_update_figure_overwrites_data(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test POST to existing key overwrites data."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_figure(session, room.id, "chart", '{"version": 1}')

    new_data = FigureData(type="plotly", data='{"version": 2}')
    response = await client.post(
        f"/v1/rooms/{room.id}/figures/chart",
        json={"figure": new_data.model_dump()},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    body = response.json()
    assert body["key"] == "chart"
    assert body["created"] is False

    # Verify updated in DB
    row = await session.get(RoomFigure, (room.id, "chart"))
    assert row is not None
    assert "2" in row.data


@pytest.mark.asyncio
async def test_create_figure_broadcasts(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """Test POST broadcasts figure:invalidate event."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    figure_data = FigureData(type="plotly", data="{}")
    await client.post(
        f"/v1/rooms/{room.id}/figures/chart",
        json={"figure": figure_data.model_dump()},
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "figure_invalidate"
    assert mock_sio.emitted[0]["room"] == f"room:{room.id}"


# =============================================================================
# Delete Figure Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_figure_removes_data(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE removes a figure."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_figure(session, room.id, "to_delete", "{}")

    response = await client.delete(
        f"/v1/rooms/{room.id}/figures/to_delete",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await session.get(RoomFigure, (room.id, "to_delete"))
    assert row is None


@pytest.mark.asyncio
async def test_delete_figure_broadcasts(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """Test DELETE broadcasts figure:invalidate event."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_figure(session, room.id, "chart", "{}")

    await client.delete(
        f"/v1/rooms/{room.id}/figures/chart",
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "figure_invalidate"


@pytest.mark.asyncio
async def test_delete_nonexistent_figure_returns_404(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE on nonexistent figure returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(
        f"/v1/rooms/{room.id}/figures/nonexistent",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == FigureNotFound.type_uri()


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_figures_public(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET without auth succeeds (public endpoint)."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(f"/v1/rooms/{room.id}/figures")
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_create_figure_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test POST without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    figure_data = FigureData(type="plotly", data="{}")
    response = await client.post(
        f"/v1/rooms/{room.id}/figures/chart",
        json={"figure": figure_data.model_dump()},
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_delete_figure_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test DELETE without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(f"/v1/rooms/{room.id}/figures/chart")
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_figures_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/figures",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_create_figure_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test POST for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    figure_data = FigureData(type="plotly", data="{}")
    response = await client.post(
        "/v1/rooms/99999/figures/chart",
        json={"figure": figure_data.model_dump()},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_delete_figure_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test DELETE for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.delete(
        "/v1/rooms/99999/figures/chart",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
