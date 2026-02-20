"""Tests for Bookmarks REST API endpoints."""

from collections.abc import AsyncIterator
from unittest.mock import AsyncMock, MagicMock

import pytest
import pytest_asyncio
from conftest import create_test_token, create_test_user_model
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from zndraw_auth import User
from zndraw_auth.settings import AuthSettings

from zndraw.config import Settings
from zndraw.exceptions import BookmarkNotFound
from zndraw.models import MemberRole, Room, RoomBookmark, RoomMembership
from zndraw.schemas import StatusResponse
from zndraw.socket_events import BookmarksInvalidate

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="bm_session")
async def bm_session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session for each test."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )

    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)

        async_session_factory = async_sessionmaker(
            bind=engine,
            class_=AsyncSession,
            expire_on_commit=False,
        )

        async with async_session_factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MagicMock:
    """Create a mock Socket.IO server for testing."""
    sio_mock = MagicMock()
    sio_mock.emit = AsyncMock()
    return sio_mock


@pytest_asyncio.fixture(name="bm_client")
async def bm_client_fixture(
    bm_session: AsyncSession,
    mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with dependencies overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield bm_session

    def get_sio_override() -> MagicMock:
        return mock_sio

    # Mock Redis for WritableRoomDep (returns None = no edit lock)
    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_tsio] = get_sio_override
    app.dependency_overrides[get_redis] = lambda: mock_redis
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client

    app.dependency_overrides.clear()


async def _create_user(
    session: AsyncSession, email: str = "testuser@local.test"
) -> tuple[User, str]:
    """Create a user and return the user and access token."""
    user = create_test_user_model(email=email)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    token = create_test_token(user)
    return user, token


async def _create_room(
    session: AsyncSession, user: User, description: str = "Test Room"
) -> Room:
    """Create a room with user as owner."""
    room = Room(
        description=description,
        created_by_id=user.id,  # type: ignore
        is_public=True,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore
        user_id=user.id,  # type: ignore
        role=MemberRole.OWNER,
    )
    session.add(membership)
    await session.commit()

    return room


async def _add_bookmark(
    session: AsyncSession, room_id: str, index: int, label: str
) -> None:
    """Add a bookmark directly to the database."""
    session.add(RoomBookmark(room_id=room_id, frame_index=index, label=label))
    await session.commit()


def _auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# List Bookmarks Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_bookmarks_returns_empty_initially(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
) -> None:
    """Test GET returns empty bookmarks for new room."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.get(
        f"/v1/rooms/{room.id}/bookmarks",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["items"] == {}


@pytest.mark.asyncio
async def test_list_bookmarks_returns_all_bookmarks(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
) -> None:
    """Test GET returns all bookmarks."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    await _add_bookmark(bm_session, room.id, 0, "Start")
    await _add_bookmark(bm_session, room.id, 5, "Middle")
    await _add_bookmark(bm_session, room.id, 10, "End")

    response = await bm_client.get(
        f"/v1/rooms/{room.id}/bookmarks",
        headers=_auth_header(token),
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
    bm_client: AsyncClient,
    bm_session: AsyncSession,
) -> None:
    """Test GET single bookmark returns label."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    await _add_bookmark(bm_session, room.id, 5, "Important Frame")

    response = await bm_client.get(
        f"/v1/rooms/{room.id}/bookmarks/5",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["index"] == 5
    assert data["label"] == "Important Frame"


@pytest.mark.asyncio
async def test_get_bookmark_returns_404_for_nonexistent(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
) -> None:
    """Test GET returns 404 for nonexistent bookmark."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.get(
        f"/v1/rooms/{room.id}/bookmarks/999",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "bookmark-not-found" in response.json()["type"]


# =============================================================================
# Set Bookmark Tests
# =============================================================================


@pytest.mark.asyncio
async def test_set_bookmark_creates_bookmark(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT creates a new bookmark."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": "New Bookmark"},
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["index"] == 3
    assert data["label"] == "New Bookmark"

    # Verify persisted in DB
    row = await bm_session.get(RoomBookmark, (room.id, 3))
    assert row is not None
    assert row.label == "New Bookmark"


@pytest.mark.asyncio
async def test_set_bookmark_broadcasts(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT broadcasts bookmarks:invalidate event."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    await bm_client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": "Test"},
        headers=_auth_header(token),
    )

    mock_sio.emit.assert_called()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, BookmarksInvalidate)
    assert call_args[1]["room"] == f"room:{room.id}"


@pytest.mark.asyncio
async def test_set_bookmark_rejects_empty_label(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
) -> None:
    """Test PUT rejects empty label."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": ""},
        headers=_auth_header(token),
    )
    assert response.status_code == 422  # Validation error


# =============================================================================
# Delete Bookmark Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_bookmark_removes_bookmark(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE removes a bookmark."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    await _add_bookmark(bm_session, room.id, 5, "To Delete")

    response = await bm_client.delete(
        f"/v1/rooms/{room.id}/bookmarks/5",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await bm_session.get(RoomBookmark, (room.id, 5))
    assert row is None


@pytest.mark.asyncio
async def test_delete_nonexistent_bookmark_returns_404(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
) -> None:
    """Test DELETE on nonexistent bookmark returns 404."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.delete(
        f"/v1/rooms/{room.id}/bookmarks/999",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == BookmarkNotFound.type_uri()


@pytest.mark.asyncio
async def test_delete_bookmark_broadcasts(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE broadcasts bookmarks:invalidate event."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    await _add_bookmark(bm_session, room.id, 5, "Test")

    await bm_client.delete(
        f"/v1/rooms/{room.id}/bookmarks/5",
        headers=_auth_header(token),
    )

    mock_sio.emit.assert_called()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, BookmarksInvalidate)


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_bookmarks_requires_auth(
    bm_client: AsyncClient, bm_session: AsyncSession
) -> None:
    """Test GET without auth returns 401."""
    user, _ = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.get(f"/v1/rooms/{room.id}/bookmarks")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_set_bookmark_requires_auth(
    bm_client: AsyncClient, bm_session: AsyncSession
) -> None:
    """Test PUT without auth returns 401."""
    user, _ = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.put(
        f"/v1/rooms/{room.id}/bookmarks/3",
        json={"label": "Test"},
    )
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_bookmarks_returns_404_for_nonexistent_room(
    bm_client: AsyncClient, bm_session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await _create_user(bm_session)

    response = await bm_client.get(
        "/v1/rooms/99999/bookmarks",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_set_bookmark_returns_404_for_nonexistent_room(
    bm_client: AsyncClient, bm_session: AsyncSession
) -> None:
    """Test PUT for non-existent room returns 404."""
    _, token = await _create_user(bm_session)

    response = await bm_client.put(
        "/v1/rooms/99999/bookmarks/3",
        json={"label": "Test"},
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
