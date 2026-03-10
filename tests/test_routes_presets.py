"""Tests for Presets REST API endpoints."""

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
from zndraw.exceptions import (
    InvalidPresetRule,
    PresetAlreadyExists,
    PresetNotFound,
)
from zndraw.models import (
    MemberRole,
    Room,
    RoomGeometry,
    RoomMembership,
    RoomPreset,
)

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="pr_session")
async def pr_session_fixture() -> AsyncIterator[AsyncSession]:
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


@pytest_asyncio.fixture(name="pr_client")
async def pr_client_fixture(
    pr_session: AsyncSession,
    mock_sio: MagicMock,
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with dependencies overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield pr_session

    def get_sio_override() -> MagicMock:
        return mock_sio

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


async def _add_preset(
    session: AsyncSession,
    room_id: str,
    name: str,
    description: str = "",
    rules: list[dict] | None = None,
) -> None:
    """Add a preset directly to the database."""
    session.add(
        RoomPreset(
            room_id=room_id,
            name=name,
            description=description,
            rules=json.dumps(rules or []),
        )
    )
    await session.commit()


async def _add_geometry(
    session: AsyncSession,
    room_id: str,
    key: str,
    geometry_type: str,
    config: dict,
) -> None:
    """Add a geometry directly to the database."""
    session.add(
        RoomGeometry(
            room_id=room_id,
            key=key,
            type=geometry_type,
            config=json.dumps(config),
        )
    )
    await session.commit()


def _auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# List Presets Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_presets_includes_bundled(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """GET /presets includes bundled presets even with no DB rows."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    response = await pr_client.get(
        f"/v1/rooms/{room.id}/presets",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    items = response.json()["items"]
    names = {p["name"] for p in items}
    assert "matt" in names
    assert "flat" in names


@pytest.mark.asyncio
async def test_list_presets_db_overrides_bundled(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """DB preset with same name as bundled takes precedence in list."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    # Override the bundled "matt" preset at room level
    await _add_preset(pr_session, room.id, "matt", "My custom matt")

    response = await pr_client.get(
        f"/v1/rooms/{room.id}/presets",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    items = response.json()["items"]
    names = [p["name"] for p in items]
    # No duplicates
    assert names.count("matt") == 1
    # The DB version wins
    pub = next(p for p in items if p["name"] == "matt")
    assert pub["description"] == "My custom matt"
    # Bundled "flat" still present
    assert any(p["name"] == "flat" for p in items)


@pytest.mark.asyncio
async def test_list_presets_returns_all(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """GET /presets returns all presets in the room."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    await _add_preset(pr_session, room.id, "custom-a", "Custom A")
    await _add_preset(pr_session, room.id, "custom-b", "Custom B")

    response = await pr_client.get(
        f"/v1/rooms/{room.id}/presets",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    items = response.json()["items"]
    names = {p["name"] for p in items}
    # DB presets + bundled presets (matt, flat, glossy, pathtracing)
    assert {"custom-a", "custom-b", "matt", "flat"}.issubset(names)


@pytest.mark.asyncio
async def test_list_presets_room_not_found(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """GET /presets returns 404 for nonexistent room."""
    user, token = await _create_user(pr_session)

    response = await pr_client.get(
        "/v1/rooms/nonexistent/presets",
        headers=_auth_header(token),
    )
    assert response.status_code == 404


# =============================================================================
# Get Preset Tests (including bundled fallback)
# =============================================================================


@pytest.mark.asyncio
async def test_get_bundled_preset(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """GET /presets/{name} returns a bundled preset when no DB row exists."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    response = await pr_client.get(
        f"/v1/rooms/{room.id}/presets/matt",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "matt"
    assert len(data["rules"]) > 0


@pytest.mark.asyncio
async def test_get_preset_returns_data(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """GET /presets/{name} returns the preset."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    rules = [{"pattern": "fog", "config": {"active": True}}]
    await _add_preset(pr_session, room.id, "test-preset", "Test", rules)

    response = await pr_client.get(
        f"/v1/rooms/{room.id}/presets/test-preset",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "test-preset"
    assert data["description"] == "Test"
    assert len(data["rules"]) == 1
    assert data["rules"][0]["pattern"] == "fog"


@pytest.mark.asyncio
async def test_get_preset_returns_404_for_nonexistent(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """GET /presets/{name} returns 404 for nonexistent preset."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    response = await pr_client.get(
        f"/v1/rooms/{room.id}/presets/nonexistent",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert PresetNotFound.type_uri() in response.json()["type"]


# =============================================================================
# Create Preset Tests (POST)
# =============================================================================


@pytest.mark.asyncio
async def test_create_preset(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """POST /presets creates a new preset."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    body = {
        "name": "my-preset",
        "description": "A test preset",
        "rules": [
            {"pattern": "fog", "config": {"active": True}},
        ],
    }

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets",
        json=body,
        headers=_auth_header(token),
    )
    assert response.status_code == 201
    data = response.json()
    assert data["name"] == "my-preset"
    assert data["description"] == "A test preset"

    # Verify persisted in DB
    row = await pr_session.get(RoomPreset, (room.id, "my-preset"))
    assert row is not None
    assert row.description == "A test preset"


@pytest.mark.asyncio
async def test_create_preset_returns_409_if_exists(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """POST /presets returns 409 if preset already exists."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    await _add_preset(pr_session, room.id, "existing")

    body = {"name": "existing", "rules": []}

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets",
        json=body,
        headers=_auth_header(token),
    )
    assert response.status_code == 409
    assert PresetAlreadyExists.type_uri() in response.json()["type"]


@pytest.mark.asyncio
async def test_create_preset_validates_geometry_type(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """POST /presets returns 422 for unknown geometry_type in rules."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    body = {
        "name": "invalid",
        "rules": [
            {
                "pattern": "*",
                "geometry_type": "NonExistentType",
                "config": {"active": True},
            },
        ],
    }

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets",
        json=body,
        headers=_auth_header(token),
    )
    assert response.status_code == 422
    assert InvalidPresetRule.type_uri() in response.json()["type"]


@pytest.mark.asyncio
async def test_create_preset_validates_config_keys(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """POST /presets returns 422 for invalid config keys on known type."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    body = {
        "name": "bad-keys",
        "rules": [
            {
                "pattern": "fog",
                "geometry_type": "Fog",
                "config": {"nonexistent_field": 42},
            },
        ],
    }

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets",
        json=body,
        headers=_auth_header(token),
    )
    assert response.status_code == 422
    assert InvalidPresetRule.type_uri() in response.json()["type"]


# =============================================================================
# Upsert Preset Tests (PUT)
# =============================================================================


@pytest.mark.asyncio
async def test_put_preset_creates_new(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """PUT /presets/{name} creates a preset if it doesn't exist."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    body = {
        "name": "new-preset",
        "description": "Created via PUT",
        "rules": [],
    }

    response = await pr_client.put(
        f"/v1/rooms/{room.id}/presets/new-preset",
        json=body,
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["name"] == "new-preset"

    row = await pr_session.get(RoomPreset, (room.id, "new-preset"))
    assert row is not None


@pytest.mark.asyncio
async def test_put_preset_updates_existing(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """PUT /presets/{name} updates an existing preset."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    await _add_preset(pr_session, room.id, "existing", "Old description")

    body = {
        "name": "existing",
        "description": "New description",
        "rules": [{"pattern": "*", "config": {"active": False}}],
    }

    response = await pr_client.put(
        f"/v1/rooms/{room.id}/presets/existing",
        json=body,
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["description"] == "New description"

    await pr_session.refresh(
        await pr_session.get(RoomPreset, (room.id, "existing"))  # type: ignore
    )
    row = await pr_session.get(RoomPreset, (room.id, "existing"))
    assert row is not None
    assert row.description == "New description"


# =============================================================================
# Delete Preset Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_preset(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """DELETE /presets/{name} removes the preset."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    await _add_preset(pr_session, room.id, "to-delete")

    response = await pr_client.delete(
        f"/v1/rooms/{room.id}/presets/to-delete",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["status"] == "ok"

    row = await pr_session.get(RoomPreset, (room.id, "to-delete"))
    assert row is None


@pytest.mark.asyncio
async def test_delete_preset_returns_404_for_nonexistent(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """DELETE /presets/{name} returns 404 if preset doesn't exist."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    response = await pr_client.delete(
        f"/v1/rooms/{room.id}/presets/nonexistent",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert PresetNotFound.type_uri() in response.json()["type"]


@pytest.mark.asyncio
async def test_delete_bundled_preset_returns_404(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """DELETE on a bundled-only preset returns 404 (no DB row to delete)."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    response = await pr_client.delete(
        f"/v1/rooms/{room.id}/presets/matt",
        headers=_auth_header(token),
    )
    assert response.status_code == 404


# =============================================================================
# Apply Preset Tests (including bundled)
# =============================================================================


@pytest.mark.asyncio
async def test_apply_bundled_preset(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """POST /presets/{name}/apply works for bundled presets without DB row."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    # Add a fog geometry so the bundled "matt" preset has something to match
    await _add_geometry(
        pr_session,
        room.id,
        "fog",
        "Fog",
        {"active": False, "near": 10.0, "far": 100.0, "color": "#000000"},
    )

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets/matt/apply",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    result = response.json()
    assert "fog" in result["geometries_updated"]

    # Verify the fog config was updated per the bundled matt preset
    geom = await pr_session.get(RoomGeometry, (room.id, "fog"))
    assert geom is not None
    config = json.loads(geom.config)
    assert config["active"] is False


# =============================================================================
# Apply Preset Tests
# =============================================================================


@pytest.mark.asyncio
async def test_apply_preset_updates_matching_geometries(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """POST /presets/{name}/apply deep-merges config into matching geometries."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    # Add a Fog geometry
    await _add_geometry(
        pr_session,
        room.id,
        "fog",
        "Fog",
        {"active": False, "near": 10.0, "far": 100.0, "color": "#000000"},
    )

    # Create preset that modifies fog
    rules = [{"pattern": "fog", "config": {"active": True, "color": "#ffffff"}}]
    await _add_preset(pr_session, room.id, "test-apply", rules=rules)

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets/test-apply/apply",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    result = response.json()
    assert "fog" in result["geometries_updated"]

    # Verify geometry was updated in DB
    geom = await pr_session.get(RoomGeometry, (room.id, "fog"))
    assert geom is not None
    config = json.loads(geom.config)
    assert config["active"] is True
    assert config["color"] == "#ffffff"
    # Verify existing fields preserved
    assert config["near"] == 10.0
    assert config["far"] == 100.0


@pytest.mark.asyncio
async def test_apply_preset_deep_merges_config(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """Apply performs deep merge — overrides specified fields, preserves others."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    # Fog geometry with multiple fields
    await _add_geometry(
        pr_session,
        room.id,
        "fog",
        "Fog",
        {"active": True, "near": 50.0, "far": 200.0, "color": "#000000"},
    )

    # Preset changes near and color, preserves active and far
    rules = [{"pattern": "fog", "config": {"near": 100.0, "color": "#ffffff"}}]
    await _add_preset(pr_session, room.id, "adjust-fog", rules=rules)

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets/adjust-fog/apply",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    geom = await pr_session.get(RoomGeometry, (room.id, "fog"))
    assert geom is not None
    config = json.loads(geom.config)
    assert config["near"] == 100.0
    assert config["color"] == "#ffffff"
    # Preserved fields
    assert config["active"] is True
    assert config["far"] == 200.0


@pytest.mark.asyncio
async def test_apply_preset_filters_by_geometry_type(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """Apply only targets geometries matching both pattern and geometry_type."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    await _add_geometry(
        pr_session, room.id, "main-light", "DirectionalLight", {"intensity": 1.0}
    )
    await _add_geometry(
        pr_session, room.id, "ambient-light", "AmbientLight", {"intensity": 0.5}
    )

    # Only target AmbientLight, even though both match *light*
    rules = [
        {
            "pattern": "*light*",
            "geometry_type": "AmbientLight",
            "config": {"intensity": 0.2},
        }
    ]
    await _add_preset(pr_session, room.id, "dim-ambient", rules=rules)

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets/dim-ambient/apply",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    result = response.json()
    assert "ambient-light" in result["geometries_updated"]
    assert "main-light" not in result["geometries_updated"]

    # Verify ambient changed
    ambient = await pr_session.get(RoomGeometry, (room.id, "ambient-light"))
    assert ambient is not None
    assert json.loads(ambient.config)["intensity"] == 0.2

    # Verify directional unchanged
    directional = await pr_session.get(RoomGeometry, (room.id, "main-light"))
    assert directional is not None
    assert json.loads(directional.config)["intensity"] == 1.0


@pytest.mark.asyncio
async def test_apply_preset_emits_geometry_invalidate(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
    mock_sio: MagicMock,
) -> None:
    """Apply emits GeometryInvalidate for each updated geometry."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    await _add_geometry(pr_session, room.id, "fog", "Fog", {"active": False})

    rules = [{"pattern": "fog", "config": {"active": True}}]
    await _add_preset(pr_session, room.id, "enable-fog", rules=rules)

    await pr_client.post(
        f"/v1/rooms/{room.id}/presets/enable-fog/apply",
        headers=_auth_header(token),
    )

    mock_sio.emit.assert_called()
    call_args = mock_sio.emit.call_args
    model = call_args[0][0]
    from zndraw.socket_events import GeometryInvalidate

    assert isinstance(model, GeometryInvalidate)
    assert model.key == "fog"
    assert model.operation == "set"


@pytest.mark.asyncio
async def test_apply_preset_returns_404_for_nonexistent(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """POST /presets/{name}/apply returns 404 for nonexistent preset."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets/nonexistent/apply",
        headers=_auth_header(token),
    )
    assert response.status_code == 404
    assert PresetNotFound.type_uri() in response.json()["type"]


@pytest.mark.asyncio
async def test_apply_preset_skips_non_matching_geometries(
    pr_client: AsyncClient,
    pr_session: AsyncSession,
) -> None:
    """Apply returns empty list when no geometries match the pattern."""
    user, token = await _create_user(pr_session)
    room = await _create_room(pr_session, user)

    await _add_geometry(pr_session, room.id, "particles", "Sphere", {"size": 1.0})

    rules = [{"pattern": "fog", "config": {"active": True}}]
    await _add_preset(pr_session, room.id, "fog-only", rules=rules)

    response = await pr_client.post(
        f"/v1/rooms/{room.id}/presets/fog-only/apply",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["geometries_updated"] == []


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_presets_requires_auth(
    pr_client: AsyncClient,
) -> None:
    """GET /presets returns 401 without auth."""
    response = await pr_client.get("/v1/rooms/any-room/presets")
    assert response.status_code == 401
