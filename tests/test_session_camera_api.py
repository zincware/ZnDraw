"""Integration tests for session camera as geometry.

Session cameras are now regular geometries with key 'cam:session:{session_id}'.
They are accessed via the standard geometry endpoints.
"""

import time
import uuid
from contextlib import contextmanager

import pytest
import requests
import socketio as sio_lib


@contextmanager
def _join_room_session(server: str, room_id: str, headers: dict):
    """Context manager that joins a room and keeps socket connected.

    The session only exists while the socket is connected.
    Yields (session_id, socket_client, headers_with_session) tuple.
    """
    jwt_token = headers["Authorization"].replace("Bearer ", "")

    # Connect socket and join room to get session ID
    sio = sio_lib.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    try:
        response = sio.call("room:join", {"roomId": room_id, "clientType": "frontend"})
        assert response["status"] == "ok"
        session_id = response["sessionId"]
        # Create headers with session ID for POST requests
        headers_with_session = {**headers, "X-Session-ID": session_id}
        yield session_id, sio, headers_with_session
    finally:
        sio.disconnect()


def _get_session_camera_key(session_id: str) -> str:
    """Get the geometry key for a session camera."""
    return f"cam:session:{session_id}"


def test_session_camera_created_on_join(server, get_jwt_auth_headers):
    """Session camera geometry is created when frontend joins a room."""
    room_id = "test-camera-created-on-join"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Join as frontend - this creates the session camera geometry
    with _join_room_session(server, room_id, headers) as (session_id, _, _):
        camera_key = _get_session_camera_key(session_id)

        # Get camera via geometry endpoint
        response = requests.get(
            f"{server}/api/rooms/{room_id}/geometries/{camera_key}",
            headers=headers,
            timeout=10,
        )

        assert response.status_code == 200
        data = response.json()
        assert "geometry" in data
        assert data["geometry"]["type"] == "Camera"

        # Verify default values
        cam_data = data["geometry"]["data"]
        assert cam_data["position"] == [0.0, 5.0, 10.0]
        assert cam_data["target"] == [0.0, 0.0, 0.0]
        assert cam_data["fov"] == 75.0
        assert cam_data["near"] == 0.1
        assert cam_data["far"] == 1000.0
        # Session cameras are protected by default
        assert cam_data["protected"] is True
        # Session cameras have helper_visible=False (don't clutter the scene)
        assert cam_data["helper_visible"] is False


def test_update_session_camera_via_geometry(server, get_jwt_auth_headers):
    """Session camera can be updated via geometry POST endpoint."""
    room_id = "test-camera-update"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    with _join_room_session(server, room_id, headers) as (
        session_id,
        _,
        session_headers,
    ):
        camera_key = _get_session_camera_key(session_id)

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
            f"{server}/api/rooms/{room_id}/geometries",
            headers=session_headers,
            json=camera_data,
            timeout=10,
        )

        assert response.status_code == 200
        assert response.json()["status"] == "success"

        # Verify the update
        response = requests.get(
            f"{server}/api/rooms/{room_id}/geometries/{camera_key}",
            headers=headers,
            timeout=10,
        )
        assert response.status_code == 200
        cam = response.json()["geometry"]["data"]
        assert cam["position"] == [10.0, 20.0, 30.0]
        assert cam["target"] == [1.0, 2.0, 3.0]
        assert cam["fov"] == 60.0


def test_session_camera_invalid_data_rejected(server, get_jwt_auth_headers):
    """POST with invalid camera data returns 400."""
    room_id = "test-camera-invalid"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    with _join_room_session(server, room_id, headers) as (
        session_id,
        _,
        session_headers,
    ):
        camera_key = _get_session_camera_key(session_id)

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
            f"{server}/api/rooms/{room_id}/geometries",
            headers=session_headers,
            json=invalid_data,
            timeout=10,
        )

        assert response.status_code == 400
        assert "error" in response.json()


def test_session_camera_not_found(server, get_jwt_auth_headers):
    """GET for nonexistent session camera returns 404."""
    room_id = "test-camera-not-found"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Try to access nonexistent session camera
    fake_session = f"fake-{uuid.uuid4().hex}"
    fake_camera_key = _get_session_camera_key(fake_session)

    response = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{fake_camera_key}",
        headers=headers,
        timeout=10,
    )
    assert response.status_code == 404


def test_cannot_create_new_session_camera(server, get_jwt_auth_headers):
    """Cannot create new session camera - reserved prefix is blocked for new keys."""
    room_id = "test-camera-create-blocked"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Join room to get valid session headers
    with _join_room_session(server, room_id, headers) as (_, _, session_headers):
        # Try to create a new session camera with a different (fake) session ID
        fake_camera_key = f"cam:session:{uuid.uuid4().hex}"
        camera_data = {
            "key": fake_camera_key,
            "type": "Camera",
            "data": {"position": [1.0, 1.0, 1.0]},
        }
        response = requests.post(
            f"{server}/api/rooms/{room_id}/geometries",
            headers=session_headers,
            json=camera_data,
            timeout=10,
        )

        assert response.status_code == 400
        assert "reserved prefix" in response.json()["error"]


def test_session_camera_isolation(server, get_jwt_auth_headers):
    """Different sessions have independent cameras."""
    room_id = "test-camera-isolation"
    headers = get_jwt_auth_headers(server)
    headers2 = get_jwt_auth_headers(server, "user2")

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Keep both sessions active
    with _join_room_session(server, room_id, headers) as (
        session_id_1,
        _,
        session_headers_1,
    ):
        with _join_room_session(server, room_id, headers2) as (
            session_id_2,
            _,
            session_headers_2,
        ):
            camera_key_1 = _get_session_camera_key(session_id_1)
            camera_key_2 = _get_session_camera_key(session_id_2)

            # Update cameras with different values
            requests.post(
                f"{server}/api/rooms/{room_id}/geometries",
                headers=session_headers_1,
                json={
                    "key": camera_key_1,
                    "type": "Camera",
                    "data": {"fov": 60.0, "position": [1.0, 1.0, 1.0]},
                },
                timeout=10,
            )

            requests.post(
                f"{server}/api/rooms/{room_id}/geometries",
                headers=session_headers_2,
                json={
                    "key": camera_key_2,
                    "type": "Camera",
                    "data": {"fov": 90.0, "position": [2.0, 2.0, 2.0]},
                },
                timeout=10,
            )

            # Verify each session has its own camera
            response1 = requests.get(
                f"{server}/api/rooms/{room_id}/geometries/{camera_key_1}",
                headers=headers,
                timeout=10,
            )
            response2 = requests.get(
                f"{server}/api/rooms/{room_id}/geometries/{camera_key_2}",
                headers=headers2,
                timeout=10,
            )

            assert response1.status_code == 200
            assert response2.status_code == 200

            cam1 = response1.json()["geometry"]["data"]
            cam2 = response2.json()["geometry"]["data"]

            assert cam1["fov"] == 60.0
            assert cam1["position"] == [1.0, 1.0, 1.0]

            assert cam2["fov"] == 90.0
            assert cam2["position"] == [2.0, 2.0, 2.0]


def test_session_camera_deleted_on_disconnect(server, get_jwt_auth_headers):
    """Session camera geometry is deleted when frontend disconnects."""
    room_id = "test-camera-deleted-on-disconnect"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Join and then disconnect
    jwt_token = headers["Authorization"].replace("Bearer ", "")
    sio = sio_lib.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)
    response = sio.call("room:join", {"roomId": room_id, "clientType": "frontend"})
    session_id = response["sessionId"]
    camera_key = _get_session_camera_key(session_id)

    # Verify camera exists while connected
    response = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{camera_key}",
        headers=headers,
        timeout=10,
    )
    assert response.status_code == 200

    # Disconnect
    sio.disconnect()

    # Give server time to process disconnect
    time.sleep(0.2)

    # Verify camera is deleted
    response = requests.get(
        f"{server}/api/rooms/{room_id}/geometries/{camera_key}",
        headers=headers,
        timeout=10,
    )
    assert response.status_code == 404


def test_list_sessions(server, get_jwt_auth_headers):
    """GET /sessions returns list of frontend sessions."""
    room_id = "test-list-sessions"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Keep session active while testing
    with _join_room_session(server, room_id, headers) as (session_id, _, _):
        # List sessions
        response = requests.get(
            f"{server}/api/rooms/{room_id}/sessions", headers=headers, timeout=10
        )

        assert response.status_code == 200
        data = response.json()
        assert "sessions" in data

        session_ids = [s["session_id"] for s in data["sessions"]]
        assert session_id in session_ids


def test_geometry_defaults_endpoint(server, get_jwt_auth_headers):
    """GET /api/schema/geometries/defaults returns Camera defaults."""
    headers = get_jwt_auth_headers(server)

    response = requests.get(
        f"{server}/api/schema/geometries/defaults", headers=headers, timeout=10
    )

    assert response.status_code == 200
    data = response.json()
    assert "defaults" in data
    assert "Camera" in data["defaults"]

    cam_defaults = data["defaults"]["Camera"]
    assert cam_defaults["fov"] == 75.0
    assert cam_defaults["near"] == 0.1
    assert cam_defaults["far"] == 1000.0
    assert cam_defaults["position"] == [0.0, 5.0, 10.0]
