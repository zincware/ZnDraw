"""Integration tests for camera sync between frontend and backend.

Tests socket events for camera state updates and programmatic control.
"""

import json

import redis
import requests

from zndraw import ZnDraw
from zndraw.app.redis_keys import RoomKeys
from zndraw.geometries import Camera


def _create_frontend_session(room_id: str, session_id: str) -> None:
    """Simulate a frontend session by adding to Redis."""
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.sadd(f"room:{room_id}:frontend_sessions", session_id)


def test_camera_state_stored_in_redis(server, connect_room):
    """Camera state updates from frontend are stored in Redis."""
    room = "test-camera-redis"

    # Create room and join via socket (keep connection alive)
    conn = connect_room(room, user="test-camera-user")
    session_id = conn.session_id
    _create_frontend_session(room, session_id)

    # Set camera via REST API (simulates what socket event does)
    camera_data = {
        "position": [5.0, 10.0, 15.0],
        "target": [1.0, 2.0, 3.0],
        "fov": 60.0,
    }
    response = requests.put(
        f"{server}/api/rooms/{room}/sessions/{session_id}/camera",
        headers=conn.headers,
        json=camera_data,
        timeout=10,
    )
    assert response.status_code == 200

    # Verify directly in Redis
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    stored = r.hget(RoomKeys(room).session_cameras(), session_id)

    assert stored is not None
    stored_data = json.loads(stored)
    assert stored_data["fov"] == 60.0


def test_programmatic_camera_control(server):
    """Python client can control frontend camera via session.camera."""
    room = "test-camera-programmatic"
    vis = ZnDraw(url=server, room=room, user="tester")

    # Simulate frontend session
    session_id = "controlled-session"
    _create_frontend_session(room, session_id)

    # Get session and set camera from Python
    session = vis.sessions[session_id]
    new_camera = Camera(
        position=(100.0, 100.0, 100.0),
        target=(50.0, 50.0, 50.0),
        fov=45.0,
    )
    session.camera = new_camera

    # Verify camera was stored (frontend would receive socket event)
    retrieved = session.camera
    assert retrieved.position == (100.0, 100.0, 100.0)
    assert retrieved.target == (50.0, 50.0, 50.0)
    assert retrieved.fov == 45.0


def test_camera_defaults_endpoint_returns_pydantic_defaults(
    server, get_jwt_auth_headers
):
    """GET /api/schema/geometries/defaults returns Camera defaults from Pydantic."""
    headers = get_jwt_auth_headers(server)

    response = requests.get(
        f"{server}/api/schema/geometries/defaults", headers=headers, timeout=10
    )

    assert response.status_code == 200
    data = response.json()
    assert "defaults" in data
    assert "Camera" in data["defaults"]

    # Verify defaults match Pydantic model
    cam_defaults = data["defaults"]["Camera"]
    expected = Camera()

    assert cam_defaults["fov"] == expected.fov
    assert cam_defaults["near"] == expected.near
    assert cam_defaults["far"] == expected.far
    assert cam_defaults["zoom"] == expected.zoom
    assert cam_defaults["show_crosshair"] == expected.show_crosshair
    assert cam_defaults["preserve_drawing_buffer"] == expected.preserve_drawing_buffer


def test_camera_defaults_endpoint_includes_all_geometry_types(
    server, get_jwt_auth_headers
):
    """GET /api/schema/geometries/defaults includes all geometry types."""
    headers = get_jwt_auth_headers(server)

    response = requests.get(
        f"{server}/api/schema/geometries/defaults", headers=headers, timeout=10
    )

    assert response.status_code == 200
    defaults = response.json()["defaults"]

    # Should include common geometry types
    assert "Camera" in defaults
    assert "Box" in defaults
    assert "Curve" in defaults


def test_camera_state_isolation_between_rooms(server):
    """Camera states in different rooms are isolated."""
    vis1 = ZnDraw(url=server, room="room-a", user="tester")
    vis2 = ZnDraw(url=server, room="room-b", user="tester")

    # Create sessions in each room
    session_id = "shared-session-id"  # Same ID, different rooms
    _create_frontend_session("room-a", session_id)
    _create_frontend_session("room-b", session_id)

    # Set different cameras in each room
    session_a = vis1.sessions[session_id]
    session_b = vis2.sessions[session_id]

    session_a.camera = Camera(fov=30.0)
    session_b.camera = Camera(fov=120.0)

    # Verify isolation
    assert session_a.camera.fov == 30.0
    assert session_b.camera.fov == 120.0


def test_camera_with_curve_attachment_serialization(server):
    """Camera with CurveAttachment position/target serializes correctly."""
    from zndraw.geometries import Curve
    from zndraw.transformations import CurveAttachment

    room = "test-camera-curve"
    vis = ZnDraw(url=server, room=room, user="tester")

    # Create a curve geometry
    vis.geometries["flight_path"] = Curve(
        position=[(0, 0, 0), (10, 10, 10), (20, 0, 0)]
    )

    # Create camera with CurveAttachment
    cam = Camera(
        position=CurveAttachment(geometry_key="flight_path", progress=0.5),
        target=(0.0, 0.0, 0.0),
    )

    session_id = "curve-camera-session"
    _create_frontend_session(room, session_id)

    session = vis.sessions[session_id]
    session.camera = cam

    # Retrieve and verify
    retrieved = session.camera
    assert isinstance(retrieved.position, dict) or isinstance(
        retrieved.position, CurveAttachment
    )


def test_session_camera_update_does_not_affect_geometries(server):
    """Updating session camera doesn't affect vis.geometries cameras."""
    room = "test-camera-vs-geometry"
    vis = ZnDraw(url=server, room=room, user="tester")

    # Create a camera geometry
    vis.geometries["my_camera"] = Camera(fov=90.0, position=(5.0, 5.0, 5.0))

    # Create session and set session camera
    session_id = "separate-session"
    _create_frontend_session(room, session_id)

    session = vis.sessions[session_id]
    session.camera = Camera(fov=30.0, position=(1.0, 1.0, 1.0))

    # Verify geometry camera unchanged
    geo_camera = vis.geometries["my_camera"]
    assert geo_camera.fov == 90.0
    assert geo_camera.position == (5.0, 5.0, 5.0)

    # Verify session camera is different
    assert session.camera.fov == 30.0
    assert session.camera.position == (1.0, 1.0, 1.0)
