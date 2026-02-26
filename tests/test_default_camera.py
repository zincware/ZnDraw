"""Tests for default camera GET/PUT endpoints, delete cleanup, and Python client."""

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

from zndraw.client import ZnDraw
from zndraw.config import Settings
from zndraw.geometries import Sphere
from zndraw.geometries.camera import Camera
from zndraw.models import MemberRole, Room, RoomGeometry, RoomMembership

# =============================================================================
# Fixtures (same pattern as test_routes_geometries.py)
# =============================================================================


@pytest_asyncio.fixture(name="session")
async def session_fixture() -> AsyncIterator[AsyncSession]:
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)

        factory = async_sessionmaker(
            bind=engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as s:
            yield s
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MagicMock:
    sio_mock = MagicMock()
    sio_mock.emit = AsyncMock()
    return sio_mock


@pytest_asyncio.fixture(name="client")
async def client_fixture(
    session: AsyncSession, mock_sio: MagicMock
) -> AsyncIterator[AsyncClient]:
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield session

    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)
    mock_redis.hgetall = AsyncMock(return_value={})
    mock_redis.hget = AsyncMock(return_value=None)
    mock_redis.hdel = AsyncMock(return_value=0)

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_tsio] = lambda: mock_sio
    app.dependency_overrides[get_redis] = lambda: mock_redis
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as c:
        yield c

    app.dependency_overrides.clear()


async def _create_user(
    session: AsyncSession, email: str = "testuser@local.test"
) -> tuple[User, str]:
    user = create_test_user_model(email=email)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user, create_test_token(user)


async def _create_room(session: AsyncSession, user: User) -> Room:
    room = Room(created_by_id=user.id, is_public=True)  # type: ignore[arg-type]
    session.add(room)
    await session.commit()
    await session.refresh(room)
    membership = RoomMembership(
        room_id=room.id,
        user_id=user.id,
        role=MemberRole.OWNER,  # type: ignore[arg-type]
    )
    session.add(membership)
    await session.commit()
    return room


async def _create_geometry(
    session: AsyncSession, room_id: str, key: str, geo_type: str = "Camera"
) -> RoomGeometry:
    config = json.dumps({"owner": None})
    row = RoomGeometry(room_id=room_id, key=key, type=geo_type, config=config)
    session.add(row)
    await session.commit()
    return row


# =============================================================================
# Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_default_camera_none(
    session: AsyncSession, client: AsyncClient
) -> None:
    """New room returns null default_camera."""
    user, token = await _create_user(session)
    room = await _create_room(session, user)

    resp = await client.get(
        f"/v1/rooms/{room.id}/default-camera",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    assert resp.json()["default_camera"] is None


@pytest.mark.asyncio
async def test_set_default_camera(session: AsyncSession, client: AsyncClient) -> None:
    """PUT with valid Camera key, then GET returns it."""
    user, token = await _create_user(session)
    room = await _create_room(session, user)
    await _create_geometry(session, room.id, "template-cam", "Camera")

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "template-cam"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    assert resp.json()["default_camera"] == "template-cam"

    resp = await client.get(
        f"/v1/rooms/{room.id}/default-camera",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.json()["default_camera"] == "template-cam"


@pytest.mark.asyncio
async def test_set_default_camera_not_found(
    session: AsyncSession, client: AsyncClient
) -> None:
    """PUT with nonexistent key returns 404."""
    user, token = await _create_user(session)
    room = await _create_room(session, user)

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "nonexistent"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 404


@pytest.mark.asyncio
async def test_set_default_camera_wrong_type(
    session: AsyncSession, client: AsyncClient
) -> None:
    """PUT with non-Camera geometry returns 400."""
    user, token = await _create_user(session)
    room = await _create_room(session, user)
    await _create_geometry(session, room.id, "my-sphere", "Sphere")

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "my-sphere"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 400


@pytest.mark.asyncio
async def test_unset_default_camera(session: AsyncSession, client: AsyncClient) -> None:
    """PUT null after setting unsets the default."""
    user, token = await _create_user(session)
    room = await _create_room(session, user)
    await _create_geometry(session, room.id, "template-cam", "Camera")

    await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "template-cam"},
        headers={"Authorization": f"Bearer {token}"},
    )

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": None},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    assert resp.json()["default_camera"] is None

    resp = await client.get(
        f"/v1/rooms/{room.id}/default-camera",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.json()["default_camera"] is None


@pytest.mark.asyncio
async def test_delete_geometry_clears_default(
    session: AsyncSession, client: AsyncClient
) -> None:
    """Deleting the default camera geometry clears the default."""
    user, token = await _create_user(session)
    room = await _create_room(session, user)
    await _create_geometry(session, room.id, "template-cam", "Camera")
    headers = {"Authorization": f"Bearer {token}"}

    # Set as default
    await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "template-cam"},
        headers=headers,
    )

    # Delete the geometry
    resp = await client.delete(
        f"/v1/rooms/{room.id}/geometries/template-cam",
        headers=headers,
    )
    assert resp.status_code == 200

    # Default should be cleared
    resp = await client.get(
        f"/v1/rooms/{room.id}/default-camera",
        headers=headers,
    )
    assert resp.json()["default_camera"] is None


# =============================================================================
# Integration tests (real server + Python client)
# =============================================================================


def test_default_camera_property(server_auth: str) -> None:
    """Test vis.default_camera get/set/unset."""
    vis = ZnDraw(url=server_auth)

    # Initially None
    assert vis.default_camera is None

    # Create a camera geometry
    vis.geometries["template-cam"] = Camera(position=(100, 10, 10))

    # Set default
    vis.default_camera = "template-cam"
    assert vis.default_camera == "template-cam"

    # Unset
    vis.default_camera = None
    assert vis.default_camera is None

    vis.disconnect()


def test_default_camera_validation(server_auth: str) -> None:
    """Test vis.default_camera validation."""
    vis = ZnDraw(url=server_auth)

    vis.geometries["my-sphere"] = Sphere(position=[(0, 0, 0)], radius=[1.0])

    # Non-existent key
    with pytest.raises(KeyError):
        vis.default_camera = "nonexistent"

    # Wrong type
    with pytest.raises(TypeError):
        vis.default_camera = "my-sphere"

    vis.disconnect()
