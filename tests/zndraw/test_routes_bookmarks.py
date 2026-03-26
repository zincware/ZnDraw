"""Tests for Bookmarks REST API endpoints."""

import pytest
import pytest_asyncio
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.exceptions import BookmarkNotFound
from zndraw.models import RoomBookmark
from zndraw.schemas import StatusResponse
from zndraw.socket_events import BookmarksInvalidate


async def _add_bookmark(
    session: AsyncSession, room_id: str, index: int, label: str
) -> None:
    """Add a bookmark directly to the database."""
    session.add(RoomBookmark(room_id=room_id, frame_index=index, label=label))
    await session.commit()


# =============================================================================
# List Bookmarks Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_bookmarks_returns_empty_initially(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns empty bookmarks for new room."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/bookmarks",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["items"] == {}


@pytest.mark.asyncio
async def test_list_bookmarks_returns_all_bookmarks(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns all bookmarks."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_bookmark(session, room.id, 0, "Start")
    await _add_bookmark(session, room.id, 5, "Middle")
    await _add_bookmark(session, room.id, 10, "End")

    response = await client.get(
        f"/v1/rooms/{room.id}/bookmarks",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["items"]["0"] == "Start"
    assert data["items"]["5"] == "Middle"
    assert data["items"]["10"] == "End"


# =============================================================================
# Get Bookmark Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_bookmark_returns_label(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET single bookmark returns label."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_bookmark(session, room.id, 5, "Important Frame")

    response = await client.get(
        f"/v1/rooms/{room.id}/bookmarks/5",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["index"] == 5
    assert data["label"] == "Important Frame"


@pytest.mark.asyncio
async def test_get_bookmark_returns_404_for_nonexistent(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns 404 for nonexistent bookmark."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/bookmarks/999",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "bookmark-not-found" in response.json()["type"]


# =============================================================================
# Set Bookmark Tests
# =============================================================================


@pytest.mark.asyncio
async def test_set_bookmark_creates_bookmark(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """Test PUT creates a new bookmark."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": "New Bookmark"},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["index"] == 3
    assert data["label"] == "New Bookmark"

    # Verify persisted in DB
    row = await session.get(RoomBookmark, (room.id, 3))
    assert row is not None
    assert row.label == "New Bookmark"


@pytest.mark.asyncio
async def test_set_bookmark_broadcasts(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """Test PUT broadcasts bookmarks:invalidate event."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": "Test"},
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "bookmarks_invalidate"
    assert mock_sio.emitted[0]["room"] == f"room:{room.id}"


@pytest.mark.asyncio
async def test_set_bookmark_rejects_empty_label(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT rejects empty label."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": ""},
        headers=auth_header(token),
    )
    assert response.status_code == 422  # Validation error


# =============================================================================
# Delete Bookmark Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_bookmark_removes_bookmark(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """Test DELETE removes a bookmark."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_bookmark(session, room.id, 5, "To Delete")

    response = await client.delete(
        f"/v1/rooms/{room.id}/bookmarks/5",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await session.get(RoomBookmark, (room.id, 5))
    assert row is None


@pytest.mark.asyncio
async def test_delete_nonexistent_bookmark_returns_404(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE on nonexistent bookmark returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(
        f"/v1/rooms/{room.id}/bookmarks/999",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == BookmarkNotFound.type_uri()


@pytest.mark.asyncio
async def test_delete_bookmark_broadcasts(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    """Test DELETE broadcasts bookmarks:invalidate event."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_bookmark(session, room.id, 5, "Test")

    await client.delete(
        f"/v1/rooms/{room.id}/bookmarks/5",
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "bookmarks_invalidate"


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_bookmarks_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(f"/v1/rooms/{room.id}/bookmarks")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_set_bookmark_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": "Test"},
    )
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_bookmarks_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/bookmarks",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_set_bookmark_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.put(
        "/v1/rooms/99999/bookmarks/3",
        json={"label": "Test"},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
