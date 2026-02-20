"""Tests for Geometry REST API endpoints."""

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
from zndraw.models import MemberRole, Room, RoomGeometry, RoomMembership
from zndraw.schemas import (
    GeometriesResponse,
    GeometryResponse,
    GeometrySelectionResponse,
    StatusResponse,
)
from zndraw.socket_events import GeometryInvalidate, SelectionInvalidate

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="geometry_session")
async def geometry_session_fixture() -> AsyncIterator[AsyncSession]:
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


@pytest_asyncio.fixture(name="geometry_client")
async def geometry_client_fixture(
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with dependencies overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield geometry_session

    def get_sio_override() -> MagicMock:
        return mock_sio

    # Mock Redis for WritableGeometryDep + session camera hash
    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)
    mock_redis.hgetall = AsyncMock(return_value={})
    mock_redis.hget = AsyncMock(return_value=None)
    mock_redis.hdel = AsyncMock(return_value=0)

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


async def _add_geometry(
    session: AsyncSession,
    room_id: str,
    key: str,
    geo_type: str,
    data: dict,
    selection: list[int] | None = None,
    owner: str | None = None,
) -> None:
    """Add a geometry directly to the database.

    Owner is stored inside the config JSON (Pydantic model field).
    """
    config = {**data}
    if owner is not None:
        config["owner"] = owner
    session.add(
        RoomGeometry(
            room_id=room_id,
            key=key,
            type=geo_type,
            config=json.dumps(config),
            selection=json.dumps(selection) if selection is not None else None,
        )
    )
    await session.commit()


def _auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# GET List Geometries Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_geometries_returns_empty_initially(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET returns empty geometries for new room."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = GeometriesResponse.model_validate(response.json())
    assert result.items == {}


@pytest.mark.asyncio
async def test_list_geometries_returns_all_geometries(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET returns all geometries in a room."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(geometry_session, room.id, "sphere1", "Sphere", {"radius": 1.0})
    await _add_geometry(geometry_session, room.id, "box1", "Box", {"width": 2.0})

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = GeometriesResponse.model_validate(response.json())
    assert len(result.items) == 2
    assert "sphere1" in result.items
    assert "box1" in result.items
    assert result.items["sphere1"].type == "Sphere"
    assert result.items["box1"].type == "Box"


@pytest.mark.asyncio
async def test_list_geometries_includes_type_schemas(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET includes geometry type schemas and defaults."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = GeometriesResponse.model_validate(response.json())
    assert result.types is not None
    assert "Sphere" in result.types.schemas
    assert "Sphere" in result.types.defaults


@pytest.mark.asyncio
async def test_list_geometries_includes_owner(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET returns owner field on geometries."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(
        geometry_session, room.id, "owned", "Sphere", {}, owner=str(user.id)
    )
    await _add_geometry(geometry_session, room.id, "shared", "Box", {})

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=_auth_header(token),
    )
    result = GeometriesResponse.model_validate(response.json())
    assert result.items["owned"].data.get("owner") == str(user.id)
    assert result.items["shared"].data.get("owner") is None


# =============================================================================
# GET Single Geometry Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_geometry_returns_geometry(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET returns a single geometry by key."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(
        geometry_session, room.id, "mysphere", "Sphere", {"radius": 1.5, "color": "red"}
    )

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries/mysphere",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = GeometryResponse.model_validate(response.json())
    assert result.key == "mysphere"
    assert result.geometry.type == "Sphere"
    assert result.geometry.data["radius"] == 1.5
    assert result.geometry.data["color"] == "red"


@pytest.mark.asyncio
async def test_get_geometry_returns_404_for_nonexistent(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET returns 404 for non-existent geometry."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries/nonexistent",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "geometry-not-found" in response.json()["type"]


# =============================================================================
# PUT Upsert Geometry Tests
# =============================================================================


@pytest.mark.asyncio
async def test_upsert_geometry_creates_new(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT creates a geometry when it doesn't exist."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/mysphere",
        json={"type": "Sphere", "data": {"radius": [2.0]}},
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["status"] == "ok"

    # Verify persisted in DB
    row = await geometry_session.get(RoomGeometry, (room.id, "mysphere"))
    assert row is not None
    assert row.type == "Sphere"


@pytest.mark.asyncio
async def test_upsert_geometry_validates_via_pydantic(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT validates geometry data through Pydantic model."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/particles",
        json={
            "type": "Sphere",
            "data": {"active": True, "radius": "arrays.radii"},
        },
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    # Verify the stored config is Pydantic-serialized (includes all defaults)
    row = await geometry_session.get(RoomGeometry, (room.id, "particles"))
    assert row is not None
    config = json.loads(row.config)
    assert "active" in config
    assert config["active"] is True


@pytest.mark.asyncio
async def test_upsert_geometry_broadcasts_set_operation(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT broadcasts geometry:invalidate with set operation."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/testkey",
        json={"type": "Sphere", "data": {}},
        headers=_auth_header(token),
    )

    mock_sio.emit.assert_called_once()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, GeometryInvalidate)
    assert model.operation == "set"
    assert model.key == "testkey"
    assert model.room_id == room.id


@pytest.mark.asyncio
async def test_upsert_geometry_updates_existing(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT to existing key updates the geometry."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(
        geometry_session, room.id, "sphere", "Sphere", {"radius": [1.0]}
    )

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/sphere",
        json={"type": "Sphere", "data": {"radius": [5.0]}},
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    # Verify updated in DB
    row = await geometry_session.get(RoomGeometry, (room.id, "sphere"))
    assert row is not None
    config = json.loads(row.config)
    assert config["radius"] == [5.0]


@pytest.mark.asyncio
async def test_upsert_geometry_rejects_invalid_data(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test PUT returns 400 when geometry data fails Pydantic validation."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/bad",
        json={"type": "Sphere", "data": {"resolution": -999}},
        headers=_auth_header(token),
    )
    assert response.status_code == 400


# =============================================================================
# DELETE Geometry Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_geometry_returns_204(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE removes a geometry and returns 204."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(geometry_session, room.id, "to_delete", "Sphere", {})

    response = await geometry_client.delete(
        f"/v1/rooms/{room.id}/geometries/to_delete",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await geometry_session.get(RoomGeometry, (room.id, "to_delete"))
    assert row is None


@pytest.mark.asyncio
async def test_delete_geometry_broadcasts_delete_operation(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE broadcasts geometry:invalidate with delete operation."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(geometry_session, room.id, "deletekey", "Sphere", {})

    await geometry_client.delete(
        f"/v1/rooms/{room.id}/geometries/deletekey",
        headers=_auth_header(token),
    )

    mock_sio.emit.assert_called_once()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, GeometryInvalidate)
    assert model.operation == "delete"
    assert model.key == "deletekey"
    assert model.room_id == room.id


@pytest.mark.asyncio
async def test_delete_nonexistent_geometry_succeeds(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE on nonexistent geometry succeeds silently (idempotent)."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.delete(
        f"/v1/rooms/{room.id}/geometries/nonexistent",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())


@pytest.mark.asyncio
async def test_delete_rejects_active_camera(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE returns 403 when camera is in use by another session."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)
    cam_key = f"cam:{user.email}:{str(user.id)[:8]}"

    # Configure mock Redis to report this camera as active
    from zndraw.app import app
    from zndraw.dependencies import get_redis

    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)
    mock_redis.hgetall = AsyncMock(return_value={"some-sid": cam_key})
    mock_redis.hget = AsyncMock(return_value=None)
    app.dependency_overrides[get_redis] = lambda: mock_redis

    response = await geometry_client.delete(
        f"/v1/rooms/{room.id}/geometries/{cam_key}",
        headers=_auth_header(token),
    )
    assert response.status_code == 403
    assert "camera that is in use" in response.json()["detail"]


# =============================================================================
# PUT Selection Tests
# =============================================================================


@pytest.mark.asyncio
async def test_update_selection_sets_indices(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT selection updates the geometry's selection column."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(geometry_session, room.id, "particles", "Sphere", {})

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/particles/selection",
        json={"indices": [0, 2, 5]},
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["status"] == "ok"

    # Verify selection in DB
    row = await geometry_session.get(RoomGeometry, (room.id, "particles"))
    assert row is not None
    assert row.selection is not None
    assert json.loads(row.selection) == [0, 2, 5]


@pytest.mark.asyncio
async def test_update_selection_broadcasts(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT selection broadcasts selection:invalidate event."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(geometry_session, room.id, "particles", "Sphere", {})

    await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/particles/selection",
        json={"indices": [1]},
        headers=_auth_header(token),
    )

    mock_sio.emit.assert_called_once()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    assert isinstance(model, SelectionInvalidate)
    assert call_args[1]["room"] == f"room:{room.id}"


@pytest.mark.asyncio
async def test_update_selection_returns_404_for_nonexistent(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test PUT selection returns 404 for non-existent geometry."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/nonexistent/selection",
        json={"indices": [0]},
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "geometry-not-found" in response.json()["type"]


# =============================================================================
# GET Selection Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_selection_returns_indices(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET selection returns stored indices."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(
        geometry_session, room.id, "sphere", "Sphere", {}, selection=[1, 2, 3]
    )

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries/sphere/selection",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = GeometrySelectionResponse.model_validate(response.json())
    assert result.key == "sphere"
    assert result.selection == [1, 2, 3]


@pytest.mark.asyncio
async def test_get_selection_returns_empty_for_no_selection(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET selection returns empty list when not set."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(geometry_session, room.id, "sphere", "Sphere", {})

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries/sphere/selection",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = GeometrySelectionResponse.model_validate(response.json())
    assert result.key == "sphere"
    assert result.selection == []


@pytest.mark.asyncio
async def test_get_selection_returns_404_for_nonexistent(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET selection returns 404 for non-existent geometry."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries/nonexistent/selection",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "geometry-not-found" in response.json()["type"]


# =============================================================================
# List Geometries Includes Selection
# =============================================================================


@pytest.mark.asyncio
async def test_list_geometries_includes_selection(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET list includes selection field on each geometry."""
    user, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(
        geometry_session, room.id, "sphere", "Sphere", {}, selection=[0, 1, 2]
    )
    await _add_geometry(geometry_session, room.id, "box", "Box", {})

    response = await geometry_client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = GeometriesResponse.model_validate(response.json())
    assert result.items["sphere"].selection == [0, 1, 2]
    assert result.items["box"].selection == []


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_geometries_public(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET without auth succeeds (public endpoint)."""
    user, _ = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.get(f"/v1/rooms/{room.id}/geometries")
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_get_geometry_public(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
) -> None:
    """Test GET single geometry without auth succeeds (public endpoint)."""
    user, _ = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    await _add_geometry(
        geometry_session, room.id, "somekey", "Sphere", {"radius": [1.0]}
    )

    response = await geometry_client.get(f"/v1/rooms/{room.id}/geometries/somekey")
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_upsert_geometry_requires_auth(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test PUT without auth returns 401."""
    user, _ = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/test",
        json={"type": "Sphere", "data": {}},
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_delete_geometry_requires_auth(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test DELETE without auth returns 401."""
    user, _ = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.delete(f"/v1/rooms/{room.id}/geometries/somekey")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_update_selection_requires_auth(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test PUT selection without auth returns 401."""
    user, _ = await _create_user(geometry_session)
    room = await _create_room(geometry_session, user)

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/particles/selection",
        json={"indices": [0]},
    )
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_geometries_returns_404_for_nonexistent_room(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await _create_user(geometry_session)

    response = await geometry_client.get(
        "/v1/rooms/99999/geometries",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_get_geometry_returns_404_for_nonexistent_room(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test GET single geometry for non-existent room returns 404."""
    _, token = await _create_user(geometry_session)

    response = await geometry_client.get(
        "/v1/rooms/99999/geometries/somekey",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_upsert_geometry_returns_404_for_nonexistent_room(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test PUT for non-existent room returns 404."""
    _, token = await _create_user(geometry_session)

    response = await geometry_client.put(
        "/v1/rooms/99999/geometries/test",
        json={"type": "Sphere", "data": {}},
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_delete_geometry_returns_404_for_nonexistent_room(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test DELETE for non-existent room returns 404."""
    _, token = await _create_user(geometry_session)

    response = await geometry_client.delete(
        "/v1/rooms/99999/geometries/somekey",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_update_selection_returns_404_for_nonexistent_room(
    geometry_client: AsyncClient, geometry_session: AsyncSession
) -> None:
    """Test PUT selection for non-existent room returns 404."""
    _, token = await _create_user(geometry_session)

    response = await geometry_client.put(
        "/v1/rooms/99999/geometries/particles/selection",
        json={"indices": [0]},
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


# =============================================================================
# Ownership Enforcement Tests
# =============================================================================


@pytest.mark.asyncio
async def test_upsert_rejects_non_owner(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT returns 403 when modifying a geometry owned by another user."""
    owner, _ = await _create_user(geometry_session, email="owner@local.test")
    _other, other_token = await _create_user(geometry_session, email="other@local.test")
    room = await _create_room(geometry_session, owner)

    await _add_geometry(
        geometry_session,
        room.id,
        "owned_sphere",
        "Sphere",
        {"radius": [1.0]},
        owner=str(owner.id),
    )

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/owned_sphere",
        json={"type": "Sphere", "data": {"radius": [99.0]}},
        headers=_auth_header(other_token),
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_upsert_allows_owner(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT succeeds when the owner modifies their own geometry."""
    owner, token = await _create_user(geometry_session)
    room = await _create_room(geometry_session, owner)

    await _add_geometry(
        geometry_session,
        room.id,
        "owned_sphere",
        "Sphere",
        {"radius": [1.0]},
        owner=str(owner.id),
    )

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/owned_sphere",
        json={"type": "Sphere", "data": {"radius": [5.0]}},
        headers=_auth_header(token),
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_upsert_allows_unowned(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT succeeds for geometry with no owner."""
    user_a, _ = await _create_user(geometry_session, email="a@local.test")
    _user_b, token_b = await _create_user(geometry_session, email="b@local.test")
    room = await _create_room(geometry_session, user_a)

    await _add_geometry(
        geometry_session, room.id, "shared_sphere", "Sphere", {"radius": [1.0]}
    )

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/shared_sphere",
        json={"type": "Sphere", "data": {"radius": [5.0]}},
        headers=_auth_header(token_b),
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_delete_rejects_non_owner(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test DELETE returns 403 when deleting geometry owned by another user."""
    owner, _ = await _create_user(geometry_session, email="owner@local.test")
    _other, other_token = await _create_user(geometry_session, email="other@local.test")
    room = await _create_room(geometry_session, owner)

    await _add_geometry(
        geometry_session,
        room.id,
        "owned_sphere",
        "Sphere",
        {},
        owner=str(owner.id),
    )

    response = await geometry_client.delete(
        f"/v1/rooms/{room.id}/geometries/owned_sphere",
        headers=_auth_header(other_token),
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_selection_update_rejects_non_owner(
    geometry_client: AsyncClient,
    geometry_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Test PUT selection returns 403 when geometry is owned by another user."""
    owner, _ = await _create_user(geometry_session, email="owner@local.test")
    _other, other_token = await _create_user(geometry_session, email="other@local.test")
    room = await _create_room(geometry_session, owner)

    await _add_geometry(
        geometry_session,
        room.id,
        "owned_sphere",
        "Sphere",
        {},
        owner=str(owner.id),
    )

    response = await geometry_client.put(
        f"/v1/rooms/{room.id}/geometries/owned_sphere/selection",
        json={"indices": [0, 1]},
        headers=_auth_header(other_token),
    )
    assert response.status_code == 403
