import uuid

import requests

from zndraw import ZnDraw
from zndraw.settings import RoomConfig


def test_settings_endpoint_returns_schema_and_defaults(
    server, join_room_and_get_headers
):
    """GET /api/rooms/{room_id}/sessions/{session_id}/settings returns schema and data.

    Note: Camera is no longer part of room settings - it's now a geometry
    with per-session state accessed via session.camera.
    """
    room_id = "test-settings-defaults"
    headers = join_room_and_get_headers(server, room_id, "test-user")
    session_id = headers["X-Session-ID"]

    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=headers,
        timeout=10,
    )
    assert response.status_code == 200

    data = response.json()

    # Check schema is present (camera is no longer in settings)
    assert "schema" in data
    assert "properties" in data["schema"]
    assert "studio_lighting" in data["schema"]["properties"]
    assert "pathtracing" in data["schema"]["properties"]
    assert "property_inspector" in data["schema"]["properties"]
    # Camera is NOT in settings anymore
    assert "camera" not in data["schema"]["properties"]

    # Check data is present with defaults
    assert "data" in data
    assert "studio_lighting" in data["data"]
    assert "pathtracing" in data["data"]
    assert "property_inspector" in data["data"]
    # Camera is NOT in settings data anymore
    assert "camera" not in data["data"]

    # Verify defaults match RoomConfig defaults
    default_config = RoomConfig()
    assert (
        data["data"]["studio_lighting"]["key_light"]
        == default_config.studio_lighting.key_light
    )
    assert data["data"]["pathtracing"]["enabled"] == default_config.pathtracing.enabled


def test_settings_endpoint_partial_update(server, connect_room):
    """PUT /api/rooms/{room_id}/sessions/{session_id}/settings supports partial updates."""
    room_id = "test-settings-partial"
    conn = connect_room(room_id, user="test-user")
    session_id = conn.session_id

    # Get initial defaults
    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    initial = response.json()["data"]
    initial_pathtracing_enabled = initial["pathtracing"]["enabled"]

    # Update only studio_lighting.key_light
    response = requests.put(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=conn.headers,
        json={"studio_lighting": {"key_light": 1.5}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify the updated field changed
    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    updated = response.json()["data"]

    # Updated field should have new value
    assert updated["studio_lighting"]["key_light"] == 1.5
    # Other categories should be unchanged (still have defaults)
    assert updated["pathtracing"]["enabled"] == initial_pathtracing_enabled


def test_settings_endpoint_multiple_categories_update(server, connect_room):
    """PUT /api/rooms/{room_id}/sessions/{session_id}/settings updates multiple categories."""
    # Use connect_room to keep socket connected during test (avoids race condition)
    conn = connect_room("test-settings-multi")

    # Update multiple categories (camera is no longer part of settings)
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/sessions/{conn.session_id}/settings",
        headers=conn.headers,
        json={
            "studio_lighting": {"key_light": 2.0, "fill_light": 1.0},
            "pathtracing": {"bounces": 5},
        },
        timeout=10,
    )
    assert response.status_code == 200

    # Verify all updates applied
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/sessions/{conn.session_id}/settings",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    updated = response.json()["data"]

    assert updated["studio_lighting"]["key_light"] == 2.0
    assert updated["studio_lighting"]["fill_light"] == 1.0
    assert updated["pathtracing"]["bounces"] == 5


def test_settings_endpoint_invalid_category_rejected(server, join_room_and_get_headers):
    """PUT /api/rooms/{room_id}/sessions/{session_id}/settings rejects unknown categories."""
    room_id = "test-settings-invalid"
    headers = join_room_and_get_headers(server, room_id, "test-user")
    session_id = headers["X-Session-ID"]

    response = requests.put(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=headers,
        json={"nonexistent_category": {"foo": "bar"}},
        timeout=10,
    )
    assert response.status_code == 400
    assert "Unknown settings categories" in response.json()["error"]


def test_settings_endpoint_camera_rejected(server, join_room_and_get_headers):
    """PUT /api/rooms/{room_id}/sessions/{session_id}/settings rejects camera category.

    Camera settings are now part of the Camera geometry model and accessed
    via session.camera, not via room settings.
    """
    room_id = "test-settings-camera-rejected"
    headers = join_room_and_get_headers(server, room_id, "test-user")
    session_id = headers["X-Session-ID"]

    response = requests.put(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=headers,
        json={"camera": {"near": 0.5}},
        timeout=10,
    )
    assert response.status_code == 400
    assert "Unknown settings categories" in response.json()["error"]


def test_settings_isolation_same_room_same_session(server, connect_room):
    """Settings changes by same session in same room should be visible."""
    room_id = "test-settings-isolation-same"
    conn = connect_room(room_id, user="test-user")
    session_id = conn.session_id

    # Set a value
    response = requests.put(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=conn.headers,
        json={"studio_lighting": {"key_light": 0.8}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify it's visible on same room same session
    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions/{session_id}/settings",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.8


def test_settings_isolation_different_room(server, join_room_and_get_headers):
    """Settings in different rooms should be isolated."""
    room_id_1 = "test-settings-isolation-room1"
    room_id_2 = "test-settings-isolation-room2"
    headers1 = join_room_and_get_headers(server, room_id_1, "test-user")
    headers2 = join_room_and_get_headers(server, room_id_2, "test-user")
    session_id_1 = headers1["X-Session-ID"]
    session_id_2 = headers2["X-Session-ID"]

    # Set a value in room 1
    response = requests.put(
        f"{server}/api/rooms/{room_id_1}/sessions/{session_id_1}/settings",
        headers=headers1,
        json={"studio_lighting": {"key_light": 0.8}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify room 1 has the value
    response = requests.get(
        f"{server}/api/rooms/{room_id_1}/sessions/{session_id_1}/settings",
        headers=headers1,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.8

    # Verify room 2 still has default
    response = requests.get(
        f"{server}/api/rooms/{room_id_2}/sessions/{session_id_2}/settings",
        headers=headers2,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.7  # Default


def test_settings_isolation_different_sessions(server, join_room_and_get_headers):
    """Settings by different sessions in same room should be isolated (per-session)."""
    room_id = "test-settings-isolation-sessions"
    headers1 = join_room_and_get_headers(server, room_id, "user1")
    headers2 = join_room_and_get_headers(server, room_id, "user2")
    session_id_1 = headers1["X-Session-ID"]
    session_id_2 = headers2["X-Session-ID"]

    # Set a value as session 1
    response = requests.put(
        f"{server}/api/rooms/{room_id}/sessions/{session_id_1}/settings",
        headers=headers1,
        json={"studio_lighting": {"key_light": 0.8}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify session 1 sees the value
    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions/{session_id_1}/settings",
        headers=headers1,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.8

    # Verify session 2 still sees default (settings are per-session)
    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions/{session_id_2}/settings",
        headers=headers2,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.7  # Default


# ==================== Python API Tests ====================


def test_session_settings_returns_roomconfig(server, connect_room):
    """vis.sessions[x].settings returns a RoomConfig instance."""
    room_id = str(uuid.uuid4())

    # Join as frontend (connect_room uses clientType: "frontend")
    _ = connect_room(room_id)

    # Create Python client to access sessions
    vis = ZnDraw(url=server, room=room_id, user=str(uuid.uuid4()))

    # Get the frontend session
    sessions = list(vis.sessions.values())
    assert len(sessions) == 1

    session = sessions[0]
    settings = session.settings

    # Verify it's a RoomConfig instance
    assert isinstance(settings, RoomConfig)

    # Verify it has expected attributes with defaults
    assert hasattr(settings, "studio_lighting")
    assert hasattr(settings, "pathtracing")
    assert hasattr(settings, "property_inspector")
    assert settings.studio_lighting.key_light == 0.7  # Default


def test_session_settings_setter_accepts_roomconfig(server, connect_room):
    """vis.sessions[x].settings setter accepts RoomConfig and persists."""
    room_id = str(uuid.uuid4())

    # Join as frontend
    _ = connect_room(room_id)

    # Create Python client
    vis = ZnDraw(url=server, room=room_id, user=str(uuid.uuid4()))

    session = list(vis.sessions.values())[0]

    # Create a new RoomConfig with modified values
    new_settings = RoomConfig()
    new_settings.studio_lighting.key_light = 1.5
    new_settings.pathtracing.enabled = True

    # Set it via the property
    session.settings = new_settings

    # Fetch again and verify persistence
    fetched = session.settings
    assert isinstance(fetched, RoomConfig)
    assert fetched.studio_lighting.key_light == 1.5
    assert fetched.pathtracing.enabled is True


def test_session_settings_attribute_modification_autosaves(server, connect_room):
    """Modifying session.settings attributes auto-saves via callback."""
    room_id = str(uuid.uuid4())

    # Join as frontend
    _ = connect_room(room_id)  # keep connection alive

    # Create Python client
    vis = ZnDraw(url=server, room=room_id, user=str(uuid.uuid4()))

    session = list(vis.sessions.values())[0]

    # Get settings and modify an attribute (should auto-save via callback)
    settings = session.settings
    settings.studio_lighting.key_light = 2.0

    # Fetch fresh and verify the change persisted
    fresh_settings = session.settings
    assert fresh_settings.studio_lighting.key_light == 2.0


def test_session_settings_multiple_attributes_autosave(server, connect_room):
    """Multiple attribute modifications each auto-save correctly."""
    room_id = str(uuid.uuid4())

    # Join as frontend
    _ = connect_room(room_id)  # keep connection alive

    # Create Python client
    vis = ZnDraw(url=server, room=room_id, user=str(uuid.uuid4()))

    session = list(vis.sessions.values())[0]

    # Modify multiple attributes
    settings = session.settings
    settings.studio_lighting.fill_light = 1.2
    settings.pathtracing.bounces = 7

    # Verify both persisted
    fresh = session.settings
    assert fresh.studio_lighting.fill_light == 1.2
    assert fresh.pathtracing.bounces == 7


def test_session_settings_validation(server, connect_room):
    """RoomConfig validates attribute values via Pydantic."""
    import pytest

    room_id = str(uuid.uuid4())

    # Join as frontend
    _ = connect_room(room_id)  # keep connection alive

    # Create Python client
    vis = ZnDraw(url=server, room=room_id, user=str(uuid.uuid4()))

    session = list(vis.sessions.values())[0]
    settings = session.settings

    # Pydantic should validate - key_light has ge=0.0, le=3.0
    with pytest.raises(Exception):
        settings.studio_lighting.key_light = 10.0  # Out of range


def test_session_settings_isolation_between_sessions(server, connect_room):
    """Different frontend sessions have isolated settings."""
    room_id = str(uuid.uuid4())

    # Join as two different frontend sessions
    _ = connect_room(room_id, user="user1")  # keep connection alive
    _ = connect_room(room_id, user="user2")  # keep connection alive

    # Create Python client
    vis = ZnDraw(url=server, room=room_id, user=str(uuid.uuid4()))

    sessions = list(vis.sessions.values())
    assert len(sessions) == 2

    # Modify settings on first session
    sessions[0].settings.studio_lighting.key_light = 1.8

    # Verify first session has the change
    assert sessions[0].settings.studio_lighting.key_light == 1.8

    # Verify second session still has default
    assert sessions[1].settings.studio_lighting.key_light == 0.7
