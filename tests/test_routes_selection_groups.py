"""Tests for Selection Groups REST API endpoints."""

import json
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
from zndraw.exceptions import SelectionGroupNotFound
from zndraw.models import MemberRole, Room, RoomMembership, SelectionGroup
from zndraw.schemas import StatusResponse
from zndraw.socket_events import SelectionGroupsInvalidate

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="sg_session")
async def sg_session_fixture() -> AsyncIterator[AsyncSession]:
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


@pytest_asyncio.fixture(name="sg_client")
async def sg_client_fixture(
    sg_session: AsyncSession,
    mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with dependencies overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield sg_session

    def get_sio_override() -> MagicMock:
        return mock_sio

    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)
    mock_redis.hget = AsyncMock(return_value=None)

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


def _auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# List Selection Groups Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_selection_groups_returns_empty(
    sg_client: AsyncClient,
    sg_session: AsyncSession,
) -> None:
    """Test GET returns empty groups for new room."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    response = await sg_client.get(
        f"/v1/rooms/{room.id}/selection-groups",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["items"] == {}


@pytest.mark.asyncio
async def test_list_selection_groups_returns_all(
    sg_client: AsyncClient,
    sg_session: AsyncSession,
) -> None:
    """Test GET returns all stored groups."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    group_a = {"sphere": [1, 2], "cube": [3]}
    group_b = {"sphere": [4, 5]}
    await _add_selection_group(sg_session, room.id, "group_a", group_a)
    await _add_selection_group(sg_session, room.id, "group_b", group_b)

    response = await sg_client.get(
        f"/v1/rooms/{room.id}/selection-groups",
        headers=_auth_header(token),
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
    sg_client: AsyncClient,
    sg_session: AsyncSession,
) -> None:
    """Test GET group returns 404 when not exists."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    response = await sg_client.get(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == SelectionGroupNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_selection_group_returns_stored(
    sg_client: AsyncClient,
    sg_session: AsyncSession,
) -> None:
    """Test GET group returns stored selections."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    group_data = {"sphere": [1, 2], "cube": [3, 4]}
    await _add_selection_group(sg_session, room.id, "mygroup", group_data)

    response = await sg_client.get(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["group"] == group_data


# =============================================================================
# Update Selection Group Tests
# =============================================================================


@pytest.mark.asyncio
async def test_update_selection_group_stores_data(
    sg_client: AsyncClient,
    sg_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT group stores selection data."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    group_data = {"sphere": [0, 1], "cube": [2, 3]}
    response = await sg_client.put(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        json={"selections": group_data},
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["status"] == "ok"

    # Verify stored in DB
    row = await sg_session.get(SelectionGroup, (room.id, "mygroup"))
    assert row is not None
    assert json.loads(row.selections) == group_data


@pytest.mark.asyncio
async def test_update_selection_group_broadcasts(
    sg_client: AsyncClient,
    sg_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT group broadcasts selection_groups:invalidate event."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    await sg_client.put(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        json={"selections": {"sphere": [0]}},
        headers=_auth_header(token),
    )

    mock_sio.emit.assert_called()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, SelectionGroupsInvalidate)


# =============================================================================
# Delete Selection Group Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_nonexistent_selection_group_returns_404(
    sg_client: AsyncClient,
    sg_session: AsyncSession,
) -> None:
    """Test DELETE on nonexistent selection group returns 404."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    response = await sg_client.delete(
        f"/v1/rooms/{room.id}/selection-groups/nonexistent",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == SelectionGroupNotFound.type_uri()


@pytest.mark.asyncio
async def test_delete_selection_group(
    sg_client: AsyncClient,
    sg_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE removes group."""
    user, token = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    await _add_selection_group(sg_session, room.id, "mygroup", {})

    response = await sg_client.delete(
        f"/v1/rooms/{room.id}/selection-groups/mygroup",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await sg_session.get(SelectionGroup, (room.id, "mygroup"))
    assert row is None


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_selection_groups_requires_auth(
    sg_client: AsyncClient, sg_session: AsyncSession
) -> None:
    """Test GET without auth returns 401."""
    user, _ = await _create_user(sg_session)
    room = await _create_room(sg_session, user)

    response = await sg_client.get(f"/v1/rooms/{room.id}/selection-groups")
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_selection_groups_returns_404_for_nonexistent_room(
    sg_client: AsyncClient, sg_session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await _create_user(sg_session)

    response = await sg_client.get(
        "/v1/rooms/99999/selection-groups",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
