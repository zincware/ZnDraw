"""Tests for Step REST API endpoints."""

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

from zndraw.models import MemberRole, Room, RoomMembership
from zndraw.schemas import StepResponse, StepUpdateResponse
from zndraw.socket_events import FrameUpdate
from zndraw.storage import InMemoryStorage

# =============================================================================
# Test-specific Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="step_session")
async def step_session_fixture() -> AsyncIterator[AsyncSession]:
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


@pytest_asyncio.fixture(name="step_storage")
async def step_storage_fixture() -> AsyncIterator[InMemoryStorage]:
    """Create a fresh InMemoryStorage instance for each test."""
    storage = InMemoryStorage()
    yield storage
    await storage.close()


@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MagicMock:
    """Create a mock Socket.IO server for testing."""
    sio_mock = MagicMock()
    sio_mock.emit = AsyncMock()
    return sio_mock


@pytest_asyncio.fixture(name="step_client")
async def step_client_fixture(
    step_session: AsyncSession,
    step_storage: InMemoryStorage,
    mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with dependencies overridden."""
    from zndraw_auth import get_session
    from zndraw_auth.settings import AuthSettings

    from zndraw.app import app
    from zndraw.config import Settings
    from zndraw.dependencies import get_redis, get_storage, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield step_session

    def get_storage_override() -> InMemoryStorage:
        return step_storage

    def get_sio_override() -> MagicMock:
        return mock_sio

    # Mock Redis for WritableRoomDep (returns None = no edit lock)
    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_storage] = get_storage_override
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
    session: AsyncSession, user: User, description: str = "Test Room", step: int = 0
) -> Room:
    """Create a room with user as owner."""
    room = Room(
        description=description,
        created_by_id=user.id,  # type: ignore
        is_public=True,
        step=step,
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


def _auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# GET Step Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_step_returns_zero_initially(
    step_client: AsyncClient,
    step_session: AsyncSession,
    step_storage: InMemoryStorage,
) -> None:
    """Test GET returns step=0 for new room with no step set."""
    user, token = await _create_user(step_session)
    room = await _create_room(step_session, user)

    # Add some frames to the room
    await step_storage.extend(room.id, [{"a": 1}, {"b": 2}, {"c": 3}])  # type: ignore

    response = await step_client.get(
        f"/v1/rooms/{room.id}/step",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = StepResponse.model_validate(response.json())
    assert result.step == 0
    assert result.total_frames == 3


@pytest.mark.asyncio
async def test_get_step_returns_current_step(
    step_client: AsyncClient,
    step_session: AsyncSession,
    step_storage: InMemoryStorage,
) -> None:
    """Test GET returns previously set step."""
    user, token = await _create_user(step_session)
    room = await _create_room(step_session, user, step=2)

    # Add frames
    await step_storage.extend(room.id, [{"a": 1}, {"b": 2}, {"c": 3}])  # type: ignore

    response = await step_client.get(
        f"/v1/rooms/{room.id}/step",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = StepResponse.model_validate(response.json())
    assert result.step == 2
    assert result.total_frames == 3


# =============================================================================
# PUT Step Tests
# =============================================================================


@pytest.mark.asyncio
async def test_set_step_updates_and_returns(
    step_client: AsyncClient,
    step_session: AsyncSession,
    step_storage: InMemoryStorage,
    mock_sio: MagicMock,
) -> None:
    """Test PUT updates step and returns new value."""
    user, token = await _create_user(step_session)
    room = await _create_room(step_session, user)

    # Add frames
    await step_storage.extend(room.id, [{"a": 1}, {"b": 2}, {"c": 3}])  # type: ignore

    response = await step_client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 1},
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = StepUpdateResponse.model_validate(response.json())
    assert result.success is True
    assert result.step == 1

    # Verify step persisted in DB
    await step_session.refresh(room)
    assert room.step == 1

    # Verify Socket.IO broadcast was sent
    mock_sio.emit.assert_called_once()
    call_args = mock_sio.emit.call_args
    assert isinstance(call_args[0][0], FrameUpdate)
    assert call_args[1]["room"] == f"room:{room.id}"


@pytest.mark.asyncio
async def test_set_step_out_of_bounds_returns_422(
    step_client: AsyncClient,
    step_session: AsyncSession,
    step_storage: InMemoryStorage,
) -> None:
    """Test PUT with step > total_frames returns 422."""
    user, token = await _create_user(step_session)
    room = await _create_room(step_session, user)

    # Add 3 frames (indices 0, 1, 2)
    await step_storage.extend(room.id, [{"a": 1}, {"b": 2}, {"c": 3}])  # type: ignore

    # Request step=100 — should return 422
    response = await step_client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 100},
        headers=_auth_header(token),
    )
    assert response.status_code == 422
    body = response.json()
    assert body["type"] == "/v1/problems/step-out-of-bounds"
    assert "out of bounds" in body["detail"]


@pytest.mark.asyncio
async def test_set_step_empty_room_rejects_nonzero(
    step_client: AsyncClient,
    step_session: AsyncSession,
    step_storage: InMemoryStorage,
) -> None:
    """Test PUT to room with no frames rejects non-zero step."""
    user, token = await _create_user(step_session)
    room = await _create_room(step_session, user)

    # Room has no frames — step=5 should be rejected
    response = await step_client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 5},
        headers=_auth_header(token),
    )
    assert response.status_code == 422
    body = response.json()
    assert body["type"] == "/v1/problems/step-out-of-bounds"


@pytest.mark.asyncio
async def test_set_step_negative_returns_422(
    step_client: AsyncClient,
    step_session: AsyncSession,
    step_storage: InMemoryStorage,
) -> None:
    """Test PUT with negative step returns 422 (Pydantic ge=0 validation)."""
    user, token = await _create_user(step_session)
    room = await _create_room(step_session, user)

    await step_storage.extend(room.id, [{"a": 1}])  # type: ignore

    response = await step_client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": -1},
        headers=_auth_header(token),
    )
    assert response.status_code == 422
    body = response.json()
    assert body["type"] == "/v1/problems/unprocessable-content"
    assert "body.step" in body["detail"]


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_step_requires_auth(
    step_client: AsyncClient, step_session: AsyncSession
) -> None:
    """Test GET without auth returns 401."""
    user, _ = await _create_user(step_session)
    room = await _create_room(step_session, user)

    response = await step_client.get(f"/v1/rooms/{room.id}/step")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_set_step_requires_auth(
    step_client: AsyncClient, step_session: AsyncSession
) -> None:
    """Test PUT without auth returns 401."""
    user, _ = await _create_user(step_session)
    room = await _create_room(step_session, user)

    response = await step_client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 1},
    )
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_step_returns_404_for_nonexistent_room(
    step_client: AsyncClient, step_session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await _create_user(step_session)

    response = await step_client.get(
        "/v1/rooms/99999/step",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_set_step_returns_404_for_nonexistent_room(
    step_client: AsyncClient, step_session: AsyncSession
) -> None:
    """Test PUT for non-existent room returns 404."""
    _, token = await _create_user(step_session)

    response = await step_client.put(
        "/v1/rooms/99999/step",
        json={"step": 1},
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
