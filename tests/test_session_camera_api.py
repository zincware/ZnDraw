"""Integration tests for session camera REST API.

Tests the /api/rooms/{room_id}/sessions/{session_id}/camera endpoints.
"""

import uuid
from contextlib import contextmanager

import pytest
import requests
import socketio as sio_lib


@contextmanager
def _join_room_session(server: str, room_id: str, headers: dict):
    """Context manager that joins a room and keeps socket connected.

    The session only exists while the socket is connected.
    Yields (session_id, socket_client) tuple.
    """
    jwt_token = headers["Authorization"].replace("Bearer ", "")

    # Connect socket and join room to get session ID
    sio = sio_lib.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    try:
        response = sio.call("room:join", {"roomId": room_id, "clientType": "frontend"})
        assert response["status"] == "ok"
        session_id = response["sessionId"]
        yield session_id, sio
    finally:
        sio.disconnect()


def test_get_session_camera_default(server, get_jwt_auth_headers):
    """GET /sessions/{session_id}/camera returns default Camera when not set."""
    room_id = "test-camera-default"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Keep session active while testing
    with _join_room_session(server, room_id, headers) as (session_id, _):
        # Get camera
        response = requests.get(
            f"{server}/api/rooms/{room_id}/sessions/{session_id}/camera",
            headers=headers,
            timeout=10,
        )

        assert response.status_code == 200
        data = response.json()
        assert "camera" in data

        # Verify default values
        cam = data["camera"]
        assert cam["position"] == [0.0, 5.0, 10.0]
        assert cam["target"] == [0.0, 0.0, 0.0]
        assert cam["fov"] == 75.0
        assert cam["near"] == 0.1
        assert cam["far"] == 1000.0


def test_set_session_camera(server, get_jwt_auth_headers):
    """PUT /sessions/{session_id}/camera stores camera data."""
    room_id = "test-camera-set"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    with _join_room_session(server, room_id, headers) as (session_id, _):
        # Set camera
        camera_data = {
            "position": [10.0, 20.0, 30.0],
            "target": [1.0, 2.0, 3.0],
            "fov": 60.0,
        }
        response = requests.put(
            f"{server}/api/rooms/{room_id}/sessions/{session_id}/camera",
            headers=headers,
            json=camera_data,
            timeout=10,
        )

        assert response.status_code == 200
        assert response.json()["status"] == "success"


def test_get_session_camera_after_set(server, get_jwt_auth_headers):
    """GET returns previously set camera data."""
    room_id = "test-camera-get-after-set"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    with _join_room_session(server, room_id, headers) as (session_id, _):
        # Set camera
        camera_data = {
            "position": [15.0, 25.0, 35.0],
            "target": [5.0, 6.0, 7.0],
            "fov": 90.0,
            "near": 0.5,
            "far": 500.0,
        }
        requests.put(
            f"{server}/api/rooms/{room_id}/sessions/{session_id}/camera",
            headers=headers,
            json=camera_data,
            timeout=10,
        )

        # Get camera
        response = requests.get(
            f"{server}/api/rooms/{room_id}/sessions/{session_id}/camera",
            headers=headers,
            timeout=10,
        )

        assert response.status_code == 200
        cam = response.json()["camera"]
        assert cam["position"] == [15.0, 25.0, 35.0]
        assert cam["target"] == [5.0, 6.0, 7.0]
        assert cam["fov"] == 90.0
        assert cam["near"] == 0.5
        assert cam["far"] == 500.0


def test_set_session_camera_invalid_data(server, get_jwt_auth_headers):
    """PUT with invalid camera data returns 400."""
    room_id = "test-camera-invalid"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    with _join_room_session(server, room_id, headers) as (session_id, _):
        # Try to set invalid camera (far < near)
        invalid_data = {
            "near": 100.0,
            "far": 10.0,  # Invalid: far < near
        }
        response = requests.put(
            f"{server}/api/rooms/{room_id}/sessions/{session_id}/camera",
            headers=headers,
            json=invalid_data,
            timeout=10,
        )

        assert response.status_code == 400
        assert "error" in response.json()


def test_session_camera_not_found(server, get_jwt_auth_headers):
    """GET/PUT for nonexistent session returns 404."""
    room_id = "test-camera-not-found"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Try to access nonexistent session (doesn't need active session)
    fake_session = f"fake-{uuid.uuid4().hex}"

    # GET
    response = requests.get(
        f"{server}/api/rooms/{room_id}/sessions/{fake_session}/camera",
        headers=headers,
        timeout=10,
    )
    assert response.status_code == 404

    # PUT
    response = requests.put(
        f"{server}/api/rooms/{room_id}/sessions/{fake_session}/camera",
        headers=headers,
        json={"fov": 60.0},
        timeout=10,
    )
    assert response.status_code == 404


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
    with _join_room_session(server, room_id, headers) as (session_id_1, _):
        with _join_room_session(server, room_id, headers2) as (session_id_2, _):
            # Set different cameras for each session
            requests.put(
                f"{server}/api/rooms/{room_id}/sessions/{session_id_1}/camera",
                headers=headers,
                json={"fov": 60.0, "position": [1.0, 1.0, 1.0]},
                timeout=10,
            )

            requests.put(
                f"{server}/api/rooms/{room_id}/sessions/{session_id_2}/camera",
                headers=headers2,
                json={"fov": 90.0, "position": [2.0, 2.0, 2.0]},
                timeout=10,
            )

            # Verify each session has its own camera
            response1 = requests.get(
                f"{server}/api/rooms/{room_id}/sessions/{session_id_1}/camera",
                headers=headers,
                timeout=10,
            )
            response2 = requests.get(
                f"{server}/api/rooms/{room_id}/sessions/{session_id_2}/camera",
                headers=headers2,
                timeout=10,
            )

            assert response1.status_code == 200
            assert response2.status_code == 200

            cam1 = response1.json()["camera"]
            cam2 = response2.json()["camera"]

            assert cam1["fov"] == 60.0
            assert cam1["position"] == [1.0, 1.0, 1.0]

            assert cam2["fov"] == 90.0
            assert cam2["position"] == [2.0, 2.0, 2.0]


def test_list_sessions(server, get_jwt_auth_headers):
    """GET /sessions returns list of frontend sessions."""
    room_id = "test-list-sessions"
    headers = get_jwt_auth_headers(server)

    # Create room
    requests.post(
        f"{server}/api/rooms", json={"roomId": room_id}, headers=headers, timeout=10
    )

    # Keep session active while testing
    with _join_room_session(server, room_id, headers) as (session_id, _):
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
