"""Tests for Figures REST API endpoints."""

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
from zndraw.exceptions import FigureNotFound
from zndraw.models import MemberRole, Room, RoomFigure, RoomMembership
from zndraw.schemas import FigureData, StatusResponse
from zndraw.socket_events import FigureInvalidate

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="fig_session")
async def fig_session_fixture() -> AsyncIterator[AsyncSession]:
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


@pytest_asyncio.fixture(name="fig_mock_sio")
async def fig_mock_sio_fixture() -> MagicMock:
    """Create a mock Socket.IO server for testing."""
    sio_mock = MagicMock()
    sio_mock.emit = AsyncMock()
    return sio_mock


@pytest_asyncio.fixture(name="fig_client")
async def fig_client_fixture(
    fig_session: AsyncSession,
    fig_mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with dependencies overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield fig_session

    def get_sio_override() -> MagicMock:
        return fig_mock_sio

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


async def _add_figure(
    session: AsyncSession, room_id: str, key: str, data: str, fig_type: str = "plotly"
) -> None:
    """Add a figure directly to the database."""
    session.add(RoomFigure(room_id=room_id, key=key, type=fig_type, data=data))
    await session.commit()


def _auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# List Figures Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_figures_returns_empty_initially(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
) -> None:
    """Test GET returns empty figures list for new room."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    response = await fig_client.get(
        f"/v1/rooms/{room.id}/figures",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["items"] == []


@pytest.mark.asyncio
async def test_list_figures_returns_all_keys(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
) -> None:
    """Test GET returns all figure keys."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    await _add_figure(fig_session, room.id, "chart1", '{"data": []}')
    await _add_figure(fig_session, room.id, "chart2", '{"data": []}')

    response = await fig_client.get(
        f"/v1/rooms/{room.id}/figures",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert set(data["items"]) == {"chart1", "chart2"}


# =============================================================================
# Get Figure Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_figure_returns_data(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
) -> None:
    """Test GET single figure returns the figure data."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    await _add_figure(
        fig_session, room.id, "my_chart", '{"data": [1, 2, 3], "layout": {}}'
    )

    response = await fig_client.get(
        f"/v1/rooms/{room.id}/figures/my_chart",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["key"] == "my_chart"
    assert data["figure"]["type"] == "plotly"
    assert data["figure"]["data"] == '{"data": [1, 2, 3], "layout": {}}'


@pytest.mark.asyncio
async def test_get_figure_returns_404_for_nonexistent(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
) -> None:
    """Test GET returns 404 for nonexistent figure."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    response = await fig_client.get(
        f"/v1/rooms/{room.id}/figures/nonexistent",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "figure-not-found" in response.json()["type"]


# =============================================================================
# Create/Update Figure Tests
# =============================================================================


@pytest.mark.asyncio
async def test_create_figure_stores_data(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
    fig_mock_sio: MagicMock,
) -> None:
    """Test POST creates a new figure."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    figure_data = FigureData(
        type="plotly", data='{"data": [], "layout": {"title": "Test"}}'
    )
    response = await fig_client.post(
        f"/v1/rooms/{room.id}/figures/new_chart",
        json={"figure": figure_data.model_dump()},
        headers=_auth_header(token),
    )
    assert response.status_code == 201
    body = response.json()
    assert body["key"] == "new_chart"
    assert body["created"] is True

    # Verify persisted in DB
    row = await fig_session.get(RoomFigure, (room.id, "new_chart"))
    assert row is not None
    assert row.type == "plotly"


@pytest.mark.asyncio
async def test_update_figure_overwrites_data(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
) -> None:
    """Test POST to existing key overwrites data."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    await _add_figure(fig_session, room.id, "chart", '{"version": 1}')

    new_data = FigureData(type="plotly", data='{"version": 2}')
    response = await fig_client.post(
        f"/v1/rooms/{room.id}/figures/chart",
        json={"figure": new_data.model_dump()},
        headers=_auth_header(token),
    )
    assert response.status_code == 201
    body = response.json()
    assert body["key"] == "chart"
    assert body["created"] is False

    # Verify updated in DB
    row = await fig_session.get(RoomFigure, (room.id, "chart"))
    assert row is not None
    assert "2" in row.data


@pytest.mark.asyncio
async def test_create_figure_broadcasts(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
    fig_mock_sio: MagicMock,
) -> None:
    """Test POST broadcasts figure:invalidate event."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    figure_data = FigureData(type="plotly", data="{}")
    await fig_client.post(
        f"/v1/rooms/{room.id}/figures/chart",
        json={"figure": figure_data.model_dump()},
        headers=_auth_header(token),
    )

    fig_mock_sio.emit.assert_called()
    call_args = fig_mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, FigureInvalidate)
    assert call_args[1]["room"] == f"room:{room.id}"


# =============================================================================
# Delete Figure Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_figure_removes_data(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
    fig_mock_sio: MagicMock,
) -> None:
    """Test DELETE removes a figure."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    await _add_figure(fig_session, room.id, "to_delete", "{}")

    response = await fig_client.delete(
        f"/v1/rooms/{room.id}/figures/to_delete",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await fig_session.get(RoomFigure, (room.id, "to_delete"))
    assert row is None


@pytest.mark.asyncio
async def test_delete_figure_broadcasts(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
    fig_mock_sio: MagicMock,
) -> None:
    """Test DELETE broadcasts figure:invalidate event."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    await _add_figure(fig_session, room.id, "chart", "{}")

    await fig_client.delete(
        f"/v1/rooms/{room.id}/figures/chart",
        headers=_auth_header(token),
    )

    fig_mock_sio.emit.assert_called()
    call_args = fig_mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, FigureInvalidate)


@pytest.mark.asyncio
async def test_delete_nonexistent_figure_returns_404(
    fig_client: AsyncClient,
    fig_session: AsyncSession,
) -> None:
    """Test DELETE on nonexistent figure returns 404."""
    user, token = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    response = await fig_client.delete(
        f"/v1/rooms/{room.id}/figures/nonexistent",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == FigureNotFound.type_uri()


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_figures_public(
    fig_client: AsyncClient, fig_session: AsyncSession
) -> None:
    """Test GET without auth succeeds (public endpoint)."""
    user, _ = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    response = await fig_client.get(f"/v1/rooms/{room.id}/figures")
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_create_figure_requires_auth(
    fig_client: AsyncClient, fig_session: AsyncSession
) -> None:
    """Test POST without auth returns 401."""
    user, _ = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    figure_data = FigureData(type="plotly", data="{}")
    response = await fig_client.post(
        f"/v1/rooms/{room.id}/figures/chart",
        json={"figure": figure_data.model_dump()},
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_delete_figure_requires_auth(
    fig_client: AsyncClient, fig_session: AsyncSession
) -> None:
    """Test DELETE without auth returns 401."""
    user, _ = await _create_user(fig_session)
    room = await _create_room(fig_session, user)

    response = await fig_client.delete(f"/v1/rooms/{room.id}/figures/chart")
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_figures_returns_404_for_nonexistent_room(
    fig_client: AsyncClient, fig_session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await _create_user(fig_session)

    response = await fig_client.get(
        "/v1/rooms/99999/figures",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_create_figure_returns_404_for_nonexistent_room(
    fig_client: AsyncClient, fig_session: AsyncSession
) -> None:
    """Test POST for non-existent room returns 404."""
    _, token = await _create_user(fig_session)

    figure_data = FigureData(type="plotly", data="{}")
    response = await fig_client.post(
        "/v1/rooms/99999/figures/chart",
        json={"figure": figure_data.model_dump()},
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_delete_figure_returns_404_for_nonexistent_room(
    fig_client: AsyncClient, fig_session: AsyncSession
) -> None:
    """Test DELETE for non-existent room returns 404."""
    _, token = await _create_user(fig_session)

    response = await fig_client.delete(
        "/v1/rooms/99999/figures/chart",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
