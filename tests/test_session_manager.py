"""Integration tests for Python client session management API.

Tests vis.sessions access to frontend browser sessions.
Note: Frontend sessions must be simulated by adding to Redis.
"""

import pytest
import redis

from zndraw import ZnDraw
from zndraw.geometries import Camera


def _create_frontend_session(room_id: str, session_id: str) -> None:
    """Simulate a frontend session by adding to Redis."""
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.sadd(f"room:{room_id}:frontend_sessions", session_id)


def test_sessions_empty_without_frontend(server):
    """vis.sessions is empty when no frontend connected."""
    room = "test-sessions-empty"
    vis = ZnDraw(url=server, room=room, user="tester")

    assert len(vis.sessions) == 0
    assert list(vis.sessions) == []


def test_sessions_lists_frontend_session(server):
    """vis.sessions lists registered frontend sessions."""
    room = "test-sessions-list"
    vis = ZnDraw(url=server, room=room, user="tester")

    # Simulate frontend session
    session_id = "frontend-abc-123"
    _create_frontend_session(room, session_id)

    assert len(vis.sessions) == 1
    assert session_id in vis.sessions


def test_sessions_getitem(server):
    """vis.sessions[session_id] returns FrontendSession."""
    room = "test-sessions-getitem"
    vis = ZnDraw(url=server, room=room, user="tester")

    session_id = "frontend-getitem-123"
    _create_frontend_session(room, session_id)

    session = vis.sessions[session_id]
    assert session.session_id == session_id


def test_sessions_getitem_not_found(server):
    """vis.sessions[invalid_id] raises KeyError."""
    room = "test-sessions-keyerror"
    vis = ZnDraw(url=server, room=room, user="tester")

    with pytest.raises(KeyError):
        _ = vis.sessions["nonexistent-session"]


def test_sessions_iteration(server):
    """Iteration over vis.sessions yields session IDs."""
    room = "test-sessions-iteration"
    vis = ZnDraw(url=server, room=room, user="tester")

    # Create multiple sessions
    session_ids = ["session-a", "session-b", "session-c"]
    for sid in session_ids:
        _create_frontend_session(room, sid)

    # Iterate and collect
    found_ids = list(vis.sessions)

    for sid in session_ids:
        assert sid in found_ids


def test_sessions_values(server):
    """vis.sessions.values() yields FrontendSession objects."""
    room = "test-sessions-values"
    vis = ZnDraw(url=server, room=room, user="tester")

    session_id = "session-values"
    _create_frontend_session(room, session_id)

    sessions = list(vis.sessions.values())
    assert len(sessions) == 1
    assert sessions[0].session_id == session_id


def test_sessions_items(server):
    """vis.sessions.items() yields (id, FrontendSession) pairs."""
    room = "test-sessions-items"
    vis = ZnDraw(url=server, room=room, user="tester")

    session_id = "session-items"
    _create_frontend_session(room, session_id)

    items = list(vis.sessions.items())
    assert len(items) == 1
    assert items[0][0] == session_id
    assert items[0][1].session_id == session_id


def test_session_camera_get_default(server):
    """session.camera returns default Camera when not set."""
    room = "test-session-camera-default"
    vis = ZnDraw(url=server, room=room, user="tester")

    session_id = "camera-default-session"
    _create_frontend_session(room, session_id)

    session = vis.sessions[session_id]
    cam = session.camera

    assert isinstance(cam, Camera)
    assert cam.fov == 75.0
    assert cam.position == (0.0, 5.0, 10.0)


def test_session_camera_set(server):
    """Setting session.camera updates the camera."""
    room = "test-session-camera-set"
    vis = ZnDraw(url=server, room=room, user="tester")

    session_id = "camera-set-session"
    _create_frontend_session(room, session_id)

    session = vis.sessions[session_id]

    # Set new camera
    new_cam = Camera(
        position=(10.0, 20.0, 30.0),
        target=(1.0, 2.0, 3.0),
        fov=60.0,
    )
    session.camera = new_cam

    # Verify it was stored
    retrieved_cam = session.camera
    assert retrieved_cam.position == (10.0, 20.0, 30.0)
    assert retrieved_cam.target == (1.0, 2.0, 3.0)
    assert retrieved_cam.fov == 60.0


def test_sessions_get_by_id(server):
    """vis.sessions.get(session_id) returns session or None."""
    room = "test-sessions-get"
    vis = ZnDraw(url=server, room=room, user="tester")

    session_id = "get-session"
    _create_frontend_session(room, session_id)

    # Found
    session = vis.sessions.get(session_id)
    assert session is not None
    assert session.session_id == session_id

    # Not found
    missing = vis.sessions.get("nonexistent")
    assert missing is None


def test_session_repr(server):
    """FrontendSession has meaningful repr."""
    room = "test-session-repr"
    vis = ZnDraw(url=server, room=room, user="tester")

    session_id = "repr-session"
    _create_frontend_session(room, session_id)

    session = vis.sessions[session_id]

    repr_str = repr(session)
    assert "FrontendSession" in repr_str
    assert session_id in repr_str


def test_sessions_repr(server):
    """FrontendSessions collection has meaningful repr."""
    room = "test-sessions-repr"
    vis = ZnDraw(url=server, room=room, user="tester")

    repr_str = repr(vis.sessions)
    assert "FrontendSessions" in repr_str


def test_multiple_sessions_independent_cameras(server):
    """Different sessions have independent camera states."""
    room = "test-multi-camera"
    vis = ZnDraw(url=server, room=room, user="tester")

    # Create two sessions
    _create_frontend_session(room, "session-a")
    _create_frontend_session(room, "session-b")

    session_a = vis.sessions["session-a"]
    session_b = vis.sessions["session-b"]

    # Set different cameras
    session_a.camera = Camera(fov=60.0, position=(1.0, 1.0, 1.0))
    session_b.camera = Camera(fov=90.0, position=(2.0, 2.0, 2.0))

    # Verify independence
    assert session_a.camera.fov == 60.0
    assert session_a.camera.position == (1.0, 1.0, 1.0)

    assert session_b.camera.fov == 90.0
    assert session_b.camera.position == (2.0, 2.0, 2.0)
