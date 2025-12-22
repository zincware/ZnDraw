import time

import requests

from zndraw import ZnDraw
from zndraw.settings import RoomConfig


def test_settings_endpoint_returns_schema_and_defaults(
    server, join_room_and_get_headers
):
    """GET /api/rooms/{room_id}/settings returns schema and data with defaults."""
    room_id = "test-settings-defaults"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    response = requests.get(f"{server}/api/rooms/{room_id}/settings", headers=headers)
    assert response.status_code == 200

    data = response.json()

    # Check schema is present
    assert "schema" in data
    assert "properties" in data["schema"]
    assert "camera" in data["schema"]["properties"]
    assert "studio_lighting" in data["schema"]["properties"]
    assert "pathtracing" in data["schema"]["properties"]
    assert "property_inspector" in data["schema"]["properties"]

    # Check data is present with defaults
    assert "data" in data
    assert "camera" in data["data"]
    assert "studio_lighting" in data["data"]
    assert "pathtracing" in data["data"]
    assert "property_inspector" in data["data"]

    # Verify defaults match RoomConfig defaults
    default_config = RoomConfig()
    assert (
        data["data"]["studio_lighting"]["key_light"]
        == default_config.studio_lighting.key_light
    )
    assert data["data"]["camera"]["near_plane"] == default_config.camera.near_plane
    assert data["data"]["pathtracing"]["enabled"] == default_config.pathtracing.enabled


def test_settings_endpoint_partial_update(server, join_room_and_get_headers):
    """PUT /api/rooms/{room_id}/settings supports partial updates."""
    room_id = "test-settings-partial"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    # Get initial defaults
    response = requests.get(f"{server}/api/rooms/{room_id}/settings", headers=headers)
    assert response.status_code == 200
    initial = response.json()["data"]
    initial_camera = initial["camera"]["near_plane"]
    initial_pathtracing = initial["pathtracing"]["enabled"]

    # Update only studio_lighting
    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={"studio_lighting": {"key_light": 1.5}},
    )
    assert response.status_code == 200

    # Verify only studio_lighting changed
    response = requests.get(f"{server}/api/rooms/{room_id}/settings", headers=headers)
    assert response.status_code == 200
    updated = response.json()["data"]

    assert updated["studio_lighting"]["key_light"] == 1.5
    assert updated["camera"]["near_plane"] == initial_camera  # Unchanged
    assert updated["pathtracing"]["enabled"] == initial_pathtracing  # Unchanged


def test_settings_endpoint_multiple_categories_update(
    server, join_room_and_get_headers
):
    """PUT /api/rooms/{room_id}/settings can update multiple categories at once."""
    room_id = "test-settings-multi"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    # Update multiple categories
    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={
            "studio_lighting": {"key_light": 2.0, "fill_light": 1.0},
            "camera": {"near_plane": 0.5, "far_plane": 500},
        },
    )
    assert response.status_code == 200

    # Verify all updates applied
    response = requests.get(f"{server}/api/rooms/{room_id}/settings", headers=headers)
    assert response.status_code == 200
    updated = response.json()["data"]

    assert updated["studio_lighting"]["key_light"] == 2.0
    assert updated["studio_lighting"]["fill_light"] == 1.0
    assert updated["camera"]["near_plane"] == 0.5
    assert updated["camera"]["far_plane"] == 500


def test_settings_endpoint_invalid_category_rejected(server, join_room_and_get_headers):
    """PUT /api/rooms/{room_id}/settings rejects unknown categories."""
    room_id = "test-settings-invalid"
    headers = join_room_and_get_headers(server, room_id, "test-user")

    response = requests.put(
        f"{server}/api/rooms/{room_id}/settings",
        headers=headers,
        json={"nonexistent_category": {"foo": "bar"}},
    )
    assert response.status_code == 400
    assert "Unknown settings categories" in response.json()["error"]


def test_vis_settings_same_room_same_user(server):
    client1 = ZnDraw(url=server, room="testroom", user="user1")
    client2 = ZnDraw(url=server, room="testroom", user="user1")

    client1.settings.studio_lighting.key_light = 0.8
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.8
    assert client2.settings.studio_lighting.key_light == 0.8

    client1.settings.studio_lighting.key_light = 0.5
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.5

    client2.settings.studio_lighting.key_light = 0.9
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.9
    assert client2.settings.studio_lighting.key_light == 0.9


def test_vis_settings_different_room_same_user(server):
    client1 = ZnDraw(url=server, room="room1", user="user1")
    client2 = ZnDraw(url=server, room="room2", user="user1")

    client1.settings.studio_lighting.key_light = 0.8
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.8
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client1.settings.studio_lighting.key_light = 0.5
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client2.settings.studio_lighting.key_light = 0.9
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.9


def test_vis_settings_same_room_different_user(server):
    client1 = ZnDraw(url=server, room="testroom", user="user1")
    client2 = ZnDraw(url=server, room="testroom", user="user2")

    client1.settings.studio_lighting.key_light = 0.8
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.8
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client1.settings.studio_lighting.key_light = 0.5
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client2.settings.studio_lighting.key_light = 0.9
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.9
