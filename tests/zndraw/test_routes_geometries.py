"""Tests for Geometry REST API endpoints."""

import json

import pytest
from helpers import auth_header, create_test_room, create_test_user_in_db
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.models import RoomGeometry
from zndraw.schemas import (
    GeometriesResponse,
    GeometryResponse,
    GeometrySelectionResponse,
    StatusResponse,
)
from zndraw.socket_events import GeometryInvalidate, SelectionInvalidate


# =============================================================================
# Geometry-specific helpers (kept from original)
# =============================================================================


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


# =============================================================================
# GET List Geometries Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_geometries_returns_empty_initially(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns empty geometries for new room."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = GeometriesResponse.model_validate(response.json())
    assert result.items == {}


@pytest.mark.asyncio
async def test_list_geometries_returns_all_geometries(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns all geometries in a room."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(session, room.id, "sphere1", "Sphere", {"radius": 1.0})
    await _add_geometry(session, room.id, "box1", "Box", {"width": 2.0})

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=auth_header(token),
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
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET includes geometry type schemas and defaults."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = GeometriesResponse.model_validate(response.json())
    assert result.types is not None
    assert "Sphere" in result.types.schemas
    assert "Sphere" in result.types.defaults


@pytest.mark.asyncio
async def test_list_geometries_includes_owner(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns owner field on geometries."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(
        session, room.id, "owned", "Sphere", {}, owner=str(user.id)
    )
    await _add_geometry(session, room.id, "shared", "Box", {})

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=auth_header(token),
    )
    result = GeometriesResponse.model_validate(response.json())
    assert result.items["owned"].data.get("owner") == str(user.id)
    assert result.items["shared"].data.get("owner") is None


# =============================================================================
# GET Single Geometry Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_geometry_returns_geometry(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns a single geometry by key."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(
        session, room.id, "mysphere", "Sphere", {"radius": 1.5, "color": "red"}
    )

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries/mysphere",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = GeometryResponse.model_validate(response.json())
    assert result.key == "mysphere"
    assert result.geometry.type == "Sphere"
    assert result.geometry.data["radius"] == 1.5
    assert result.geometry.data["color"] == "red"


@pytest.mark.asyncio
async def test_get_geometry_returns_404_for_nonexistent(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET returns 404 for non-existent geometry."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries/nonexistent",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "geometry-not-found" in response.json()["type"]


# =============================================================================
# PUT Upsert Geometry Tests
# =============================================================================


@pytest.mark.asyncio
async def test_upsert_geometry_creates_new(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT creates a geometry when it doesn't exist."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/mysphere",
        json={"type": "Sphere", "data": {"radius": [2.0]}},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["status"] == "ok"

    # Verify persisted in DB
    row = await session.get(RoomGeometry, (room.id, "mysphere"))
    assert row is not None
    assert row.type == "Sphere"


@pytest.mark.asyncio
async def test_upsert_geometry_validates_via_pydantic(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT validates geometry data through Pydantic model."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/particles",
        json={
            "type": "Sphere",
            "data": {"active": True, "radius": "arrays.radii"},
        },
        headers=auth_header(token),
    )
    assert response.status_code == 200

    # Verify the stored config is Pydantic-serialized (includes all defaults)
    row = await session.get(RoomGeometry, (room.id, "particles"))
    assert row is not None
    config = json.loads(row.config)
    assert "active" in config
    assert config["active"] is True


@pytest.mark.asyncio
async def test_upsert_geometry_broadcasts_set_operation(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio,
) -> None:
    """Test PUT broadcasts geometry:invalidate with set operation."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await client.put(
        f"/v1/rooms/{room.id}/geometries/testkey",
        json={"type": "Sphere", "data": {}},
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    emitted = mock_sio.emitted[0]
    assert emitted["event"] == "geometry_invalidate"
    model = GeometryInvalidate.model_validate(emitted["data"])
    assert model.operation == "set"
    assert model.key == "testkey"
    assert model.room_id == room.id


@pytest.mark.asyncio
async def test_upsert_geometry_updates_existing(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT to existing key updates the geometry."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(
        session, room.id, "sphere", "Sphere", {"radius": [1.0]}
    )

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/sphere",
        json={"type": "Sphere", "data": {"radius": [5.0]}},
        headers=auth_header(token),
    )
    assert response.status_code == 200

    # Verify updated in DB
    row = await session.get(RoomGeometry, (room.id, "sphere"))
    assert row is not None
    config = json.loads(row.config)
    assert config["radius"] == [5.0]


@pytest.mark.asyncio
async def test_upsert_geometry_rejects_invalid_data(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT returns 400 when geometry data fails Pydantic validation."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/bad",
        json={"type": "Sphere", "data": {"resolution": -999}},
        headers=auth_header(token),
    )
    assert response.status_code == 400


# =============================================================================
# DELETE Geometry Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_geometry_returns_204(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE removes a geometry and returns 204."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(session, room.id, "to_delete", "Sphere", {})

    response = await client.delete(
        f"/v1/rooms/{room.id}/geometries/to_delete",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify deleted from DB
    row = await session.get(RoomGeometry, (room.id, "to_delete"))
    assert row is None


@pytest.mark.asyncio
async def test_delete_geometry_broadcasts_delete_operation(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio,
) -> None:
    """Test DELETE broadcasts geometry:invalidate with delete operation."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(session, room.id, "deletekey", "Sphere", {})

    await client.delete(
        f"/v1/rooms/{room.id}/geometries/deletekey",
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    emitted = mock_sio.emitted[0]
    assert emitted["event"] == "geometry_invalidate"
    model = GeometryInvalidate.model_validate(emitted["data"])
    assert model.operation == "delete"
    assert model.key == "deletekey"
    assert model.room_id == room.id


@pytest.mark.asyncio
async def test_delete_nonexistent_geometry_succeeds(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE on nonexistent geometry succeeds silently (idempotent)."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(
        f"/v1/rooms/{room.id}/geometries/nonexistent",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())


@pytest.mark.asyncio
async def test_delete_rejects_active_camera(
    client: AsyncClient,
    session: AsyncSession,
    redis_client,
) -> None:
    """Test DELETE returns 403 when camera is in use by another session."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    cam_key = f"cam:{user.email}:{str(user.id)[:8]}"

    # Seed real Redis to report this camera as active.
    # The route checks redis.hgetall(RedisKey.active_cameras(room_id))
    # which resolves to "room:{room_id}:active-cameras".
    await redis_client.hset(f"room:{room.id}:active-cameras", "some-sid", cam_key)

    response = await client.delete(
        f"/v1/rooms/{room.id}/geometries/{cam_key}",
        headers=auth_header(token),
    )
    assert response.status_code == 403
    assert "camera that is in use" in response.json()["detail"]


# =============================================================================
# PUT Selection Tests
# =============================================================================


@pytest.mark.asyncio
async def test_update_selection_sets_indices(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT selection updates the geometry's selection column."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(session, room.id, "particles", "Sphere", {})

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/particles/selection",
        json={"indices": [0, 2, 5]},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["status"] == "ok"

    # Verify selection in DB
    row = await session.get(RoomGeometry, (room.id, "particles"))
    assert row is not None
    assert row.selection is not None
    assert json.loads(row.selection) == [0, 2, 5]


@pytest.mark.asyncio
async def test_update_selection_broadcasts(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio,
) -> None:
    """Test PUT selection broadcasts selection:invalidate event."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(session, room.id, "particles", "Sphere", {})

    await client.put(
        f"/v1/rooms/{room.id}/geometries/particles/selection",
        json={"indices": [1]},
        headers=auth_header(token),
    )

    assert len(mock_sio.emitted) == 1
    emitted = mock_sio.emitted[0]
    assert emitted["event"] == "selection_invalidate"
    SelectionInvalidate.model_validate(emitted["data"])
    assert emitted["room"] == f"room:{room.id}"


@pytest.mark.asyncio
async def test_update_selection_returns_404_for_nonexistent(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT selection returns 404 for non-existent geometry."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/nonexistent/selection",
        json={"indices": [0]},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "geometry-not-found" in response.json()["type"]


# =============================================================================
# GET Selection Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_selection_returns_indices(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET selection returns stored indices."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(
        session, room.id, "sphere", "Sphere", {}, selection=[1, 2, 3]
    )

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries/sphere/selection",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = GeometrySelectionResponse.model_validate(response.json())
    assert result.key == "sphere"
    assert result.selection == [1, 2, 3]


@pytest.mark.asyncio
async def test_get_selection_returns_empty_for_no_selection(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET selection returns empty list when not set."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(session, room.id, "sphere", "Sphere", {})

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries/sphere/selection",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = GeometrySelectionResponse.model_validate(response.json())
    assert result.key == "sphere"
    assert result.selection == []


@pytest.mark.asyncio
async def test_get_selection_returns_404_for_nonexistent(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET selection returns 404 for non-existent geometry."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries/nonexistent/selection",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "geometry-not-found" in response.json()["type"]


# =============================================================================
# List Geometries Includes Selection
# =============================================================================


@pytest.mark.asyncio
async def test_list_geometries_includes_selection(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET list includes selection field on each geometry."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(
        session, room.id, "sphere", "Sphere", {}, selection=[0, 1, 2]
    )
    await _add_geometry(session, room.id, "box", "Box", {})

    response = await client.get(
        f"/v1/rooms/{room.id}/geometries",
        headers=auth_header(token),
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
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET without auth succeeds (public endpoint)."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(f"/v1/rooms/{room.id}/geometries")
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_get_geometry_public(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test GET single geometry without auth succeeds (public endpoint)."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_geometry(
        session, room.id, "somekey", "Sphere", {"radius": [1.0]}
    )

    response = await client.get(f"/v1/rooms/{room.id}/geometries/somekey")
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_upsert_geometry_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/test",
        json={"type": "Sphere", "data": {}},
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_delete_geometry_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test DELETE without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(f"/v1/rooms/{room.id}/geometries/somekey")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_update_selection_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT selection without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/particles/selection",
        json={"indices": [0]},
    )
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_geometries_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/geometries",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_get_geometry_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET single geometry for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/geometries/somekey",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_upsert_geometry_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.put(
        "/v1/rooms/99999/geometries/test",
        json={"type": "Sphere", "data": {}},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_delete_geometry_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test DELETE for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.delete(
        "/v1/rooms/99999/geometries/somekey",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_update_selection_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT selection for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.put(
        "/v1/rooms/99999/geometries/particles/selection",
        json={"indices": [0]},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


# =============================================================================
# Ownership Enforcement Tests
# =============================================================================


@pytest.mark.asyncio
async def test_upsert_rejects_non_owner(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT returns 403 when modifying a geometry owned by another user."""
    owner, _ = await create_test_user_in_db(session, email="owner@local.test")
    _other, other_token = await create_test_user_in_db(session, email="other@local.test")
    room = await create_test_room(session, owner)

    await _add_geometry(
        session,
        room.id,
        "owned_sphere",
        "Sphere",
        {"radius": [1.0]},
        owner=str(owner.id),
    )

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/owned_sphere",
        json={"type": "Sphere", "data": {"radius": [99.0]}},
        headers=auth_header(other_token),
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_upsert_allows_owner(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT succeeds when the owner modifies their own geometry."""
    owner, token = await create_test_user_in_db(session)
    room = await create_test_room(session, owner)

    await _add_geometry(
        session,
        room.id,
        "owned_sphere",
        "Sphere",
        {"radius": [1.0]},
        owner=str(owner.id),
    )

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/owned_sphere",
        json={"type": "Sphere", "data": {"radius": [5.0]}},
        headers=auth_header(token),
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_upsert_allows_unowned(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT succeeds for geometry with no owner."""
    user_a, _ = await create_test_user_in_db(session, email="a@local.test")
    _user_b, token_b = await create_test_user_in_db(session, email="b@local.test")
    room = await create_test_room(session, user_a)

    await _add_geometry(
        session, room.id, "shared_sphere", "Sphere", {"radius": [1.0]}
    )

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/shared_sphere",
        json={"type": "Sphere", "data": {"radius": [5.0]}},
        headers=auth_header(token_b),
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_delete_rejects_non_owner(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test DELETE returns 403 when deleting geometry owned by another user."""
    owner, _ = await create_test_user_in_db(session, email="owner@local.test")
    _other, other_token = await create_test_user_in_db(session, email="other@local.test")
    room = await create_test_room(session, owner)

    await _add_geometry(
        session,
        room.id,
        "owned_sphere",
        "Sphere",
        {},
        owner=str(owner.id),
    )

    response = await client.delete(
        f"/v1/rooms/{room.id}/geometries/owned_sphere",
        headers=auth_header(other_token),
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_selection_update_rejects_non_owner(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT selection returns 403 when geometry is owned by another user."""
    owner, _ = await create_test_user_in_db(session, email="owner@local.test")
    _other, other_token = await create_test_user_in_db(session, email="other@local.test")
    room = await create_test_room(session, owner)

    await _add_geometry(
        session,
        room.id,
        "owned_sphere",
        "Sphere",
        {},
        owner=str(owner.id),
    )

    response = await client.put(
        f"/v1/rooms/{room.id}/geometries/owned_sphere/selection",
        json={"indices": [0, 1]},
        headers=auth_header(other_token),
    )
    assert response.status_code == 403
