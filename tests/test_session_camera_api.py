"""Integration tests for session camera as geometry.

Session cameras are now regular geometries with key 'cam:session:{session_id}'.
They are accessed via the standard geometry endpoints.
"""

import time
import uuid

import pytest
import requests


def _camera_key(session_id: str) -> str:
    """Get geometry key for a session camera."""
    return f"cam:session:{session_id}"


def test_session_camera_created_on_join(server, connect_room):
    """Session camera geometry is created when frontend joins a room."""
    conn = connect_room("test-camera-created-on-join")
    camera_key = _camera_key(conn.session_id)

    # Get camera via geometry endpoint
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/geometries/{camera_key}",
        headers=conn.headers,
        timeout=10,
    )

    assert response.status_code == 200
    data = response.json()
    assert "geometry" in data
    assert data["geometry"]["type"] == "Camera"

    # Verify default values
    cam_data = data["geometry"]["data"]
    assert cam_data["position"] == [-10.0, 10.0, 30.0]
    assert cam_data["target"] == [0.0, 0.0, 0.0]
    assert cam_data["fov"] == 50.0
    assert cam_data["near"] == 0.1
    assert cam_data["far"] == 1000.0
    # Session cameras are protected by default
    assert cam_data["protected"] is True
    # Session cameras have helper_visible=False (don't clutter the scene)
    assert cam_data["helper_visible"] is False


def test_update_session_camera_via_geometry(server, connect_room):
    """Session camera can be updated via geometry POST endpoint."""
    conn = connect_room("test-camera-update")
    camera_key = _camera_key(conn.session_id)

    # Update camera via geometry POST endpoint
    camera_data = {
        "key": camera_key,
        "type": "Camera",
        "data": {
            "position": [10.0, 20.0, 30.0],
            "target": [1.0, 2.0, 3.0],
            "fov": 60.0,
        },
    }
    response = requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )

    assert response.status_code == 200
    assert response.json()["status"] == "success"

    # Verify the update
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/geometries/{camera_key}",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    cam = response.json()["geometry"]["data"]
    assert cam["position"] == [10.0, 20.0, 30.0]
    assert cam["target"] == [1.0, 2.0, 3.0]
    assert cam["fov"] == 60.0


def test_session_camera_invalid_data_rejected(server, connect_room):
    """POST with invalid camera data returns 400."""
    conn = connect_room("test-camera-invalid")
    camera_key = _camera_key(conn.session_id)

    # Try to set invalid camera (far < near)
    invalid_data = {
        "key": camera_key,
        "type": "Camera",
        "data": {
            "near": 100.0,
            "far": 10.0,  # Invalid: far < near
        },
    }
    response = requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=invalid_data,
        timeout=10,
    )

    assert response.status_code == 400
    assert "error" in response.json()


def test_session_camera_not_found(server, connect_room):
    """GET for nonexistent session camera returns 404."""
    conn = connect_room("test-camera-not-found")

    # Try to access nonexistent session camera
    fake_session = f"fake-{uuid.uuid4().hex}"
    fake_camera_key = _camera_key(fake_session)

    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/geometries/{fake_camera_key}",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 404


def test_cannot_create_new_session_camera(server, connect_room):
    """Cannot create new session camera - reserved prefix is blocked for new keys."""
    conn = connect_room("test-camera-create-blocked")

    # Try to create a new session camera with a different (fake) session ID
    fake_camera_key = f"cam:session:{uuid.uuid4().hex}"
    camera_data = {
        "key": fake_camera_key,
        "type": "Camera",
        "data": {"position": [1.0, 1.0, 1.0]},
    }
    response = requests.post(
        f"{server}/api/rooms/{conn.room_id}/geometries",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )

    assert response.status_code == 400
    assert "reserved prefix" in response.json()["error"]


def test_multiple_sessions_have_unique_camera_keys(server, connect_room):
    """Two sessions in the same room get different camera keys."""
    room_id = "test-camera-unique-keys"
    conn1 = connect_room(room_id, user="user1")
    conn2 = connect_room(room_id, user="user2")

    camera_key_1 = _camera_key(conn1.session_id)
    camera_key_2 = _camera_key(conn2.session_id)

    # Session IDs should be different
    assert conn1.session_id != conn2.session_id

    # Camera keys should be different
    assert camera_key_1 != camera_key_2

    # Both cameras should exist
    response1 = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{camera_key_1}",
        headers=conn1.headers,
        timeout=10,
    )
    response2 = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{camera_key_2}",
        headers=conn2.headers,
        timeout=10,
    )

    assert response1.status_code == 200
    assert response2.status_code == 200


def test_session_camera_updates_are_independent(server, connect_room):
    """Updates to one session's camera don't affect another session's camera."""
    room_id = "test-camera-independence"
    conn1 = connect_room(room_id, user="user1")
    conn2 = connect_room(room_id, user="user2")

    camera_key_1 = _camera_key(conn1.session_id)
    camera_key_2 = _camera_key(conn2.session_id)

    # Update first camera
    r1 = requests.post(
        f"{server}/api/rooms/{room_id}/geometries",
        headers=conn1.headers,
        json={
            "key": camera_key_1,
            "type": "Camera",
            "data": {"fov": 60.0, "position": [1.0, 1.0, 1.0]},
        },
        timeout=10,
    )
    assert r1.status_code == 200

    # Update second camera with different values
    r2 = requests.post(
        f"{server}/api/rooms/{room_id}/geometries",
        headers=conn2.headers,
        json={
            "key": camera_key_2,
            "type": "Camera",
            "data": {"fov": 90.0, "position": [2.0, 2.0, 2.0]},
        },
        timeout=10,
    )
    assert r2.status_code == 200

    # Verify first camera has its values (not overwritten by second)
    response1 = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{camera_key_1}",
        headers=conn1.headers,
        timeout=10,
    )
    cam1 = response1.json()["geometry"]["data"]
    assert cam1["fov"] == 60.0
    assert cam1["position"] == [1.0, 1.0, 1.0]

    # Verify second camera has its values (not affected by first)
    response2 = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{camera_key_2}",
        headers=conn2.headers,
        timeout=10,
    )
    cam2 = response2.json()["geometry"]["data"]
    assert cam2["fov"] == 90.0
    assert cam2["position"] == [2.0, 2.0, 2.0]


def test_session_camera_deleted_on_disconnect(server, connect_room):
    """Session camera geometry is deleted when frontend disconnects."""
    room_id = "test-camera-deleted-on-disconnect"
    conn = connect_room(room_id)
    camera_key = _camera_key(conn.session_id)

    # Verify camera exists while connected
    response = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{camera_key}",
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200

    # Disconnect
    conn.sio.disconnect()

    # Poll until camera is deleted (max 5 seconds)
    max_wait = 5.0
    poll_interval = 0.1
    elapsed = 0.0
    while elapsed < max_wait:
        response = requests.get(
            f"{server}/api/rooms/{room_id}/geometries/{camera_key}",
            headers=conn.headers,
            timeout=10,
        )
        if response.status_code == 404:
            break
        time.sleep(poll_interval)
        elapsed += poll_interval
    else:
        pytest.fail(f"Session camera not deleted within {max_wait}s after disconnect")

    assert response.status_code == 404


def test_list_sessions(server, connect_room):
    """GET /sessions returns list of frontend sessions."""
    room_id = "test-list-sessions"
    conn = connect_room(room_id)

    # List sessions
    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions", headers=conn.headers, timeout=10
    )

    assert response.status_code == 200
    data = response.json()
    assert "sessions" in data

    session_ids = [s["session_id"] for s in data["sessions"]]
    assert conn.session_id in session_ids
