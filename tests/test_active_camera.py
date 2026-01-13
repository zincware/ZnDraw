"""Integration tests for active camera feature.

Tests the active_camera property on sessions that controls which camera
the frontend session views through.
"""

import pytest

from zndraw import ZnDraw
from zndraw.geometries import Camera, Sphere


def test_active_camera_default_is_session_camera(server, connect_room):
    """Default active_camera should be the session's own camera key."""
    room = "test-active-camera-default"
    vis = ZnDraw(url=server, room=room, user="tester")

    conn = connect_room(room, user="frontend-user")
    session = vis.sessions[conn.session_id]

    # Default should be session's own camera
    expected_key = f"cam:session:{conn.session_id}"
    assert session.active_camera == expected_key


def test_active_camera_setter(server, connect_room):
    """Setting active_camera should persist."""
    room = "test-active-camera-setter"
    vis = ZnDraw(url=server, room=room, user="tester")

    conn = connect_room(room, user="frontend-user")
    session = vis.sessions[conn.session_id]

    # Create a camera geometry to switch to
    vis.geometries["my_camera"] = Camera(position=(10.0, 20.0, 30.0))

    # Set active camera
    session.active_camera = "my_camera"

    # Verify it was stored
    assert session.active_camera == "my_camera"


def test_active_camera_switch_back_to_session_camera(server, connect_room):
    """Setting active_camera back to session camera should work."""
    room = "test-active-camera-switch-back"
    vis = ZnDraw(url=server, room=room, user="tester")

    conn = connect_room(room, user="frontend-user")
    session = vis.sessions[conn.session_id]

    # Create a camera geometry
    vis.geometries["other_camera"] = Camera(position=(5.0, 5.0, 5.0))

    # Switch to other camera
    session.active_camera = "other_camera"
    assert session.active_camera == "other_camera"

    # Switch back to session camera
    session_camera_key = session.camera_key
    session.active_camera = session_camera_key
    assert session.active_camera == session_camera_key


def test_active_camera_isolation_between_sessions(server, connect_room):
    """Different sessions have independent active_camera values."""
    room = "test-active-camera-isolation"
    vis = ZnDraw(url=server, room=room, user="tester")

    # Create two frontend sessions
    conn_a = connect_room(room, user="user-a")
    conn_b = connect_room(room, user="user-b")

    session_a = vis.sessions[conn_a.session_id]
    session_b = vis.sessions[conn_b.session_id]

    # Create a camera geometry
    vis.geometries["shared_camera"] = Camera(position=(1.0, 2.0, 3.0))

    # Set different active cameras
    session_a.active_camera = "shared_camera"
    # session_b should still have its default (own camera)

    assert session_a.active_camera == "shared_camera"
    assert session_b.active_camera == f"cam:session:{conn_b.session_id}"


def test_active_camera_persists_across_api_calls(server, connect_room):
    """Active camera should persist across multiple API calls."""
    room = "test-active-camera-persist"
    vis = ZnDraw(url=server, room=room, user="tester")

    conn = connect_room(room, user="frontend-user")

    # Set via first session object
    session1 = vis.sessions[conn.session_id]
    vis.geometries["camera1"] = Camera(position=(1.0, 1.0, 1.0))
    session1.active_camera = "camera1"

    # Get via fresh session object
    session2 = vis.sessions[conn.session_id]
    assert session2.active_camera == "camera1"


def test_active_camera_nonexistent_raises_keyerror(server, connect_room):
    """Setting active_camera to non-existent camera should raise KeyError."""
    room = "test-active-camera-nonexistent"
    vis = ZnDraw(url=server, room=room, user="tester")

    conn = connect_room(room, user="frontend-user")
    session = vis.sessions[conn.session_id]

    with pytest.raises(KeyError, match="not found in geometries"):
        session.active_camera = "nonexistent_camera"


def test_active_camera_non_camera_raises_typeerror(server, connect_room):
    """Setting active_camera to non-Camera geometry should raise TypeError."""
    room = "test-active-camera-non-camera"
    vis = ZnDraw(url=server, room=room, user="tester")

    conn = connect_room(room, user="frontend-user")
    session = vis.sessions[conn.session_id]

    # Create a non-Camera geometry (Sphere uses lists for position/radius)
    vis.geometries["my_sphere"] = Sphere()

    with pytest.raises(TypeError, match="is not a Camera"):
        session.active_camera = "my_sphere"


def test_camera_setter_updates_active_camera_geometry(server, connect_room):
    """Setting session.camera should update the active camera geometry."""
    room = "test-camera-setter-updates-active"
    vis = ZnDraw(url=server, room=room, user="tester")

    conn = connect_room(room, user="frontend-user")
    session = vis.sessions[conn.session_id]

    # Create another camera and switch to it
    vis.geometries["other_camera"] = Camera(position=(10.0, 10.0, 10.0))
    session.active_camera = "other_camera"
    assert session.active_camera == "other_camera"

    # Now set the camera - should update "other_camera", not session camera
    session.camera = Camera(position=(5.0, 5.0, 5.0))

    # The "other_camera" geometry should be updated
    updated_cam = vis.geometries["other_camera"]
    assert updated_cam.position == (5.0, 5.0, 5.0)

    # active_camera should still be "other_camera"
    assert session.active_camera == "other_camera"
