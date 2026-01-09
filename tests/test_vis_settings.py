import time

import requests

from zndraw.settings import RoomConfig


def test_settings_endpoint_returns_schema_and_defaults(
    server, join_room_and_get_headers
):
    """GET /api/rooms/{room_id}/settings returns schema and data with defaults.

    Note: Camera is no longer part of room settings - it's now a geometry
    with per-session state accessed via session.camera.
    """
    room_id = "test-settings-defaults"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    response = requests.get(
        f"{server}/api/rooms/{room_id}/settings", headers=headers, timeout=10
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


def test_settings_endpoint_partial_update(server, join_room_and_get_headers):
    """PUT /api/rooms/{room_id}/settings supports partial updates."""
    room_id = "test-settings-partial"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    # Get initial defaults
    response = requests.get(
        f"{server}/api/rooms/{room_id}/settings", headers=headers, timeout=10
    )
    assert response.status_code == 200
    initial = response.json()["data"]
    initial_pathtracing_enabled = initial["pathtracing"]["enabled"]

    # Update only studio_lighting.key_light
    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={"studio_lighting": {"key_light": 1.5}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify the updated field changed
    response = requests.get(
        f"{server}/api/rooms/{room_id}/settings", headers=headers, timeout=10
    )
    assert response.status_code == 200
    updated = response.json()["data"]

    # Updated field should have new value
    assert updated["studio_lighting"]["key_light"] == 1.5
    # Other categories should be unchanged (still have defaults)
    assert updated["pathtracing"]["enabled"] == initial_pathtracing_enabled


def test_settings_endpoint_multiple_categories_update(
    server, join_room_and_get_headers
):
    """PUT /api/rooms/{room_id}/settings can update multiple categories at once."""
    room_id = "test-settings-multi"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    # Update multiple categories (camera is no longer part of settings)
    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={
            "studio_lighting": {"key_light": 2.0, "fill_light": 1.0},
            "pathtracing": {"bounces": 5},
        },
        timeout=10,
    )
    assert response.status_code == 200

    # Verify all updates applied
    response = requests.get(
        f"{server}/api/rooms/{room_id}/settings", headers=headers, timeout=10
    )
    assert response.status_code == 200
    updated = response.json()["data"]

    assert updated["studio_lighting"]["key_light"] == 2.0
    assert updated["studio_lighting"]["fill_light"] == 1.0
    assert updated["pathtracing"]["bounces"] == 5


def test_settings_endpoint_invalid_category_rejected(server, join_room_and_get_headers):
    """PUT /api/rooms/{room_id}/settings rejects unknown categories."""
    room_id = "test-settings-invalid"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={"nonexistent_category": {"foo": "bar"}},
        timeout=10,
    )
    assert response.status_code == 400
    assert "Unknown settings categories" in response.json()["error"]


def test_settings_endpoint_camera_rejected(server, join_room_and_get_headers):
    """PUT /api/rooms/{room_id}/settings rejects camera category (moved to session).

    Camera settings are now part of the Camera geometry model and accessed
    via session.camera, not via room settings.
    """
    room_id = "test-settings-camera-rejected"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={"camera": {"near": 0.5}},
        timeout=10,
    )
    assert response.status_code == 400
    assert "Unknown settings categories" in response.json()["error"]


def test_settings_isolation_same_room_same_user(server, join_room_and_get_headers):
    """Settings changes by same user in same room should be visible."""
    room_id = "test-settings-isolation-same"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    # Set a value
    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={"studio_lighting": {"key_light": 0.8}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify it's visible on same room same user
    response = requests.get(
        f"{server}/api/rooms/{room_id}/settings", headers=headers, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.8


def test_settings_isolation_different_room_same_user(server, join_room_and_get_headers):
    """Settings in different rooms should be isolated."""
    room_id_1 = "test-settings-isolation-room1"
    room_id_2 = "test-settings-isolation-room2"
    headers = join_room_and_get_headers(server, room_id_1, "test-user")
    headers2 = join_room_and_get_headers(server, room_id_2, "test-user")

    # Set a value in room 1
    response = requests.put(
        f"{server}/api/rooms/{room_id_1}/settings",
        headers=headers,
        json={"studio_lighting": {"key_light": 0.8}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify room 1 has the value
    response = requests.get(
        f"{server}/api/rooms/{room_id_1}/settings", headers=headers, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.8

    # Verify room 2 still has default
    response = requests.get(
        f"{server}/api/rooms/{room_id_2}/settings", headers=headers2, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.7  # Default


def test_settings_isolation_same_room_different_user(server, join_room_and_get_headers):
    """Settings by different users in same room should be isolated (per-user)."""
    room_id = "test-settings-isolation-users"
    headers1 = join_room_and_get_headers(server, room_id, "user1")
    headers2 = join_room_and_get_headers(server, room_id, "user2")

    # Set a value as user1
    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers1,
        json={"studio_lighting": {"key_light": 0.8}},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify user1 sees the value
    response = requests.get(
        f"{server}/api/rooms/{room_id}/settings", headers=headers1, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.8

    # Verify user2 still sees default (settings are per-user)
    response = requests.get(
        f"{server}/api/rooms/{room_id}/settings", headers=headers2, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["data"]["studio_lighting"]["key_light"] == 0.7  # Default
