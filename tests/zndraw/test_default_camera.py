"""Tests for default camera GET/PUT endpoints, delete cleanup, and Python client."""

import json

import pytest
from helpers import auth_header, create_test_room, create_test_user_in_db
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.client import ZnDraw
from zndraw.geometries import Sphere
from zndraw.geometries.camera import Camera
from zndraw.models import RoomGeometry

# =============================================================================
# Helpers unique to this test file
# =============================================================================


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
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    resp = await client.get(
        f"/v1/rooms/{room.id}/default-camera",
        headers=auth_header(token),
    )
    assert resp.status_code == 200
    assert resp.json()["default_camera"] is None


@pytest.mark.asyncio
async def test_set_default_camera(session: AsyncSession, client: AsyncClient) -> None:
    """PUT with valid Camera key, then GET returns it."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _create_geometry(session, room.id, "template-cam", "Camera")

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "template-cam"},
        headers=auth_header(token),
    )
    assert resp.status_code == 200
    assert resp.json()["default_camera"] == "template-cam"

    resp = await client.get(
        f"/v1/rooms/{room.id}/default-camera",
        headers=auth_header(token),
    )
    assert resp.json()["default_camera"] == "template-cam"


@pytest.mark.asyncio
async def test_set_default_camera_not_found(
    session: AsyncSession, client: AsyncClient
) -> None:
    """PUT with nonexistent key returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "nonexistent"},
        headers=auth_header(token),
    )
    assert resp.status_code == 404


@pytest.mark.asyncio
async def test_set_default_camera_wrong_type(
    session: AsyncSession, client: AsyncClient
) -> None:
    """PUT with non-Camera geometry returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _create_geometry(session, room.id, "my-sphere", "Sphere")

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "my-sphere"},
        headers=auth_header(token),
    )
    assert resp.status_code == 400


@pytest.mark.asyncio
async def test_unset_default_camera(session: AsyncSession, client: AsyncClient) -> None:
    """PUT null after setting unsets the default."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _create_geometry(session, room.id, "template-cam", "Camera")

    await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": "template-cam"},
        headers=auth_header(token),
    )

    resp = await client.put(
        f"/v1/rooms/{room.id}/default-camera",
        json={"default_camera": None},
        headers=auth_header(token),
    )
    assert resp.status_code == 200
    assert resp.json()["default_camera"] is None

    resp = await client.get(
        f"/v1/rooms/{room.id}/default-camera",
        headers=auth_header(token),
    )
    assert resp.json()["default_camera"] is None


@pytest.mark.asyncio
async def test_delete_geometry_clears_default(
    session: AsyncSession, client: AsyncClient
) -> None:
    """Deleting the default camera geometry clears the default."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _create_geometry(session, room.id, "template-cam", "Camera")
    headers = auth_header(token)

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
