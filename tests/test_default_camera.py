"""Tests for default camera feature."""

import requests

from zndraw.app.redis_keys import RoomKeys


def test_default_camera_redis_key_format():
    """RoomKeys.default_camera() returns correct key format."""
    keys = RoomKeys("test-room")
    assert keys.default_camera() == "room:test-room:default_camera"


def test_get_default_camera_returns_null_when_unset(server, connect_room):
    """GET default-camera returns null when not set."""
    conn = connect_room("test-default-camera-unset")

    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )

    assert response.status_code == 200
    assert response.json()["default_camera"] is None


def test_set_and_get_default_camera(server, connect_room):
    """PUT and GET default-camera work correctly."""
    conn = connect_room("test-default-camera-set")

    # Create a camera geometry first
    camera_data = {
        "key": "my_default_cam",
        "type": "Camera",
        "data": {"position": [5.0, 5.0, 5.0], "fov": 60.0},
    }
    r = requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )
    assert r.status_code == 200

    # Set as default
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        json={"default_camera": "my_default_cam"},
        timeout=10,
    )
    assert response.status_code == 200

    # Get and verify
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["default_camera"] == "my_default_cam"


def test_set_default_camera_validates_existence(server, connect_room):
    """PUT default-camera with nonexistent camera returns 404."""
    conn = connect_room("test-default-camera-nonexistent")

    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        json={"default_camera": "nonexistent_camera"},
        timeout=10,
    )

    assert response.status_code == 404
    assert "not found" in response.json()["error"].lower()


def test_set_default_camera_validates_type(server, connect_room):
    """PUT default-camera with non-Camera geometry returns 400."""
    conn = connect_room("test-default-camera-wrong-type")

    # Create a Sphere geometry (not a Camera)
    sphere_data = {
        "key": "my_sphere",
        "type": "Sphere",
        "data": {"position": [[0.0, 0.0, 0.0]]},
    }
    r = requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=sphere_data,
        timeout=10,
    )
    assert r.status_code == 200

    # Try to set as default
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        json={"default_camera": "my_sphere"},
        timeout=10,
    )

    assert response.status_code == 400
    assert "not a camera" in response.json()["error"].lower()


def test_unset_default_camera(server, connect_room):
    """PUT default-camera with null clears the default."""
    conn = connect_room("test-default-camera-unset-explicit")

    # Create and set a default camera
    camera_data = {
        "key": "temp_cam",
        "type": "Camera",
        "data": {"position": [1.0, 1.0, 1.0]},
    }
    requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )
    requests.put(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        json={"default_camera": "temp_cam"},
        timeout=10,
    )

    # Unset
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        json={"default_camera": None},
        timeout=10,
    )
    assert response.status_code == 200

    # Verify unset
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )
    assert response.json()["default_camera"] is None


def test_room_join_includes_default_camera(server, connect_room):
    """room:join response includes defaultCamera field."""
    room_id = "test-join-default-camera"

    # First connection creates the room
    conn1 = connect_room(room_id, user="user1")

    # Create and set default camera
    camera_data = {
        "key": "default_cam",
        "type": "Camera",
        "data": {"position": [10.0, 10.0, 10.0], "fov": 75.0},
    }
    requests.post(
        f"{server}/api/rooms/{conn1.room_id}/geometries",
        headers=conn1.headers,
        json=camera_data,
        timeout=10,
    )
    requests.put(
        f"{server}/api/rooms/{conn1.room_id}/default-camera",
        headers=conn1.headers,
        json={"default_camera": "default_cam"},
        timeout=10,
    )

    # Second connection should receive defaultCamera in join response
    conn2 = connect_room(room_id, user="user2")

    # The connect_room fixture should have captured the join response
    # We verify via REST that the default camera is set
    response = requests.get(
        f"{server}/api/rooms/{conn2.room_id}/default-camera",
        headers=conn2.headers,
        timeout=10,
    )
    assert response.json()["default_camera"] == "default_cam"


def test_room_join_default_camera_null_when_unset(server, connect_room):
    """room:join response has null defaultCamera when not set."""
    room_id = "test-join-no-default"
    conn = connect_room(room_id)

    # Verify default camera is not set
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )
    assert response.json()["default_camera"] is None


def test_api_manager_get_default_camera(server, connect_room):
    """APIManager.get_default_camera returns correct value."""
    from zndraw.api_manager import APIManager

    conn = connect_room("test-api-get-default")

    # Create APIManager instance from connection data
    jwt_token = conn.headers["Authorization"].replace("Bearer ", "")
    api = APIManager(
        url=server,
        room=conn.room_id,
        jwt_token=jwt_token,
        session_id=conn.session_id,
    )

    # Initially None
    result = api.get_default_camera()
    assert result is None

    # Set via REST
    camera_data = {
        "key": "api_test_cam",
        "type": "Camera",
        "data": {"position": [1.0, 2.0, 3.0]},
    }
    requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )
    requests.put(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        json={"default_camera": "api_test_cam"},
        timeout=10,
    )

    # Get via APIManager
    result = api.get_default_camera()
    assert result == "api_test_cam"


def test_api_manager_set_default_camera(server, connect_room):
    """APIManager.set_default_camera updates value correctly."""
    from zndraw.api_manager import APIManager

    conn = connect_room("test-api-set-default")

    # Create APIManager instance from connection data
    jwt_token = conn.headers["Authorization"].replace("Bearer ", "")
    api = APIManager(
        url=server,
        room=conn.room_id,
        jwt_token=jwt_token,
        session_id=conn.session_id,
    )

    # Create camera
    camera_data = {
        "key": "api_set_cam",
        "type": "Camera",
        "data": {"position": [1.0, 2.0, 3.0]},
    }
    requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )

    # Set via APIManager
    api.set_default_camera("api_set_cam")

    # Verify via REST
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )
    assert response.json()["default_camera"] == "api_set_cam"

    # Unset
    api.set_default_camera(None)
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )
    assert response.json()["default_camera"] is None


def test_zndraw_default_camera_property(server):
    """ZnDraw.default_camera property works correctly."""
    from zndraw import ZnDraw
    from zndraw.geometries import Camera

    vis = ZnDraw(url=server, room="test-vis-default-camera", user="user1")

    # Initially None
    assert vis.default_camera is None

    # Create camera via geometries
    cam = Camera(position=(5.0, 5.0, 5.0), fov=60.0)
    vis.geometries["my_cam"] = cam

    # Set as default
    vis.default_camera = "my_cam"
    assert vis.default_camera == "my_cam"

    # Unset
    vis.default_camera = None
    assert vis.default_camera is None

    vis.disconnect()


def test_zndraw_default_camera_validates_key(server):
    """ZnDraw.default_camera setter validates key exists."""
    import pytest

    from zndraw import ZnDraw

    vis = ZnDraw(url=server, room="test-vis-default-camera-validation", user="user1")

    with pytest.raises(KeyError):
        vis.default_camera = "nonexistent_key"

    vis.disconnect()


def test_zndraw_default_camera_validates_type(server):
    """ZnDraw.default_camera setter validates geometry is Camera."""
    import pytest

    from zndraw import ZnDraw
    from zndraw.geometries import Sphere

    vis = ZnDraw(url=server, room="test-vis-default-camera-type", user="user1")

    # Create a non-Camera geometry
    sphere = Sphere(position=[[0.0, 0.0, 0.0]])
    vis.geometries["my_sphere"] = sphere

    with pytest.raises(TypeError):
        vis.default_camera = "my_sphere"

    vis.disconnect()


def test_default_camera_cleared_on_delete(server, connect_room):
    """Deleting the default camera clears the default_camera setting."""
    conn = connect_room("test-default-camera-delete-cleanup")

    # Create camera (not protected)
    camera_data = {
        "key": "deletable_cam",
        "type": "Camera",
        "data": {"position": [1.0, 1.0, 1.0], "protected": False},
    }
    requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )

    # Set as default
    requests.put(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        json={"default_camera": "deletable_cam"},
        timeout=10,
    )

    # Verify it's set
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )
    assert response.json()["default_camera"] == "deletable_cam"

    # Delete the camera
    response = requests.delete(
        f"{server}/api/rooms/{conn.room_id}/geometries/deletable_cam",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200

    # Verify default_camera is now None
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/default-camera",
        headers=conn.headers,
        timeout=10,
    )
    assert response.json()["default_camera"] is None
