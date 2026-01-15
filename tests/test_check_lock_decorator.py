"""Tests for the @check_lock decorator.

The @check_lock decorator checks for blocking locks but does NOT require the caller
to hold a lock. Operations proceed if no lock exists (FIFO handles ordering).
Only blocks if ANOTHER session holds a lock on the target.
"""

import pytest
import redis
import requests


@pytest.fixture
def room_with_lock(server, connect_room):
    """Create a room, join it, acquire a lock, and return connection info."""
    room = "test-lock-room"

    # Use connect_room to keep socket alive
    conn = connect_room(room)

    # Acquire lock for trajectory:meta
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "testing lock"},
        headers=conn.headers,
    )
    assert response.status_code == 200
    lock_data = response.json()
    assert lock_data["success"] is True
    lock_token = lock_data["lockToken"]

    return server, room, conn.session_id, conn.headers, lock_token


def test_check_lock_missing_session_id(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that decorator rejects requests without X-Session-ID header."""
    room = "test-missing-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create room first
    create_and_join_room(server, room, auth_headers)

    # Try to create geometry without session ID header
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers=auth_headers,  # Missing X-Session-ID
    )

    assert response.status_code == 400
    data = response.json()
    assert "X-Session-ID" in data["error"]


def test_check_lock_invalid_session_id(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that decorator rejects requests with invalid session ID."""
    room = "test-invalid-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create room first
    create_and_join_room(server, room, auth_headers)

    # Try to create geometry with invalid session ID
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={**auth_headers, "X-Session-ID": "invalid-session-id-12345"},
    )

    assert response.status_code == 401
    data = response.json()
    assert "Invalid or expired session" in data["error"]


def test_check_lock_no_lock_proceeds(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that operations proceed when no lock exists (FIFO ordering).

    This is the key behavioral difference from @requires_lock.
    """
    room = "test-no-lock-proceeds"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create and join room
    session_id = create_and_join_room(server, room, auth_headers)

    # Create geometry WITHOUT acquiring lock first - should succeed
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test_sphere",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    # Verify geometry was created
    response = requests.get(f"{server}/api/rooms/{room}/geometries/test_sphere")
    assert response.status_code == 200


def test_check_lock_session_user_mismatch(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that decorator rejects requests where session doesn't match JWT user."""
    room = "test-session-mismatch"

    # User 1 creates room and joins
    user1_headers = get_jwt_auth_headers(server, "user1")
    user1_session_id = create_and_join_room(server, room, user1_headers)

    # User 2 gets auth token
    user2_headers = get_jwt_auth_headers(server, "user2")

    # Try to use user1's session with user2's JWT token
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={**user2_headers, "X-Session-ID": user1_session_id},
    )

    assert response.status_code == 403
    data = response.json()
    assert "Session/user mismatch" in data["error"]


def test_check_lock_valid_session_and_lock(room_with_lock):
    """Test that decorator allows requests with valid session when caller holds lock."""
    server, room, session_id, auth_headers, _ = room_with_lock

    # Create geometry with valid session and lock
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test_sphere",
            "type": "Sphere",
            "data": {
                "position": [[1.0, 2.0, 3.0]],
                "color": ["#FF0000"],
                "radius": [1.0],
            },
        },
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"


def test_check_lock_blocks_when_other_session_holds_lock(server, s22, connect_room):
    """Test that decorator blocks when ANOTHER session holds the lock.

    Uses /api/rooms/<room_id>/frames POST which checks target="trajectory:meta".
    """
    import ase.io

    room = "test-other-session-holds-lock"

    # Session 1: Create room, join and acquire lock (socket stays connected)
    conn1 = connect_room(room, user="test-user")

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "session 1 lock"},
        headers=conn1.headers,
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Session 2: Same user joins in a different session (simulating another browser tab)
    conn2 = connect_room(room, user="test-user")

    # Prepare frame data
    frame_bytes = ase.io.write("-", s22[0], format="extxyz")

    # Session 2 tries to append frames - should be blocked since session 1 holds lock
    response = requests.post(
        f"{server}/api/rooms/{room}/frames",
        files={"frames": ("frames.xyz", frame_bytes)},
        headers=conn2.headers,
    )

    assert response.status_code == 423  # Locked by another session
    data = response.json()
    assert "locked" in data["error"].lower() or "test-user" in data["error"]


def test_check_lock_allows_same_session_with_lock(server, connect_room):
    """Test that lock holder can perform operations on locked resource."""
    room = "test-same-session-with-lock"

    # Create room, join and acquire geometry lock
    conn = connect_room(room, user="test-user")

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/geometry:test_sphere/acquire",
        json={"msg": "testing"},
        headers=conn.headers,
    )
    assert response.status_code == 200

    # Same session can create geometry for the locked key
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test_sphere",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers=conn.headers,
    )

    assert response.status_code == 200


def test_check_lock_after_lock_release_proceeds(server, connect_room):
    """Test that operations proceed after lock is released."""
    room = "test-lock-release"

    # Session 1: Create room, join and acquire geometry lock
    conn1 = connect_room(room, user="test-user")

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/geometry:test_sphere/acquire",
        json={"msg": "session 1 lock"},
        headers=conn1.headers,
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # Session 2: Same user joins in a different session
    conn2 = connect_room(room, user="test-user")

    # Session 1 releases the lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/geometry:test_sphere/release",
        json={"lockToken": lock_token},
        headers=conn1.headers,
    )
    assert response.status_code == 200

    # Now session 2 can proceed (no lock held)
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test_sphere",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers=conn2.headers,
    )
    assert response.status_code == 200


def test_check_lock_after_lock_expiry_proceeds(server, connect_room):
    """Test that operations proceed after lock expires.

    With @check_lock, if no lock exists, operations proceed.
    """
    room = "test-lock-expiry"

    # Create room and join
    conn = connect_room(room, user="test-user")

    # Acquire geometry lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/geometry:test_sphere/acquire",
        json={"msg": "testing expiry"},
        headers=conn.headers,
    )
    assert response.status_code == 200

    # Manually expire the lock in Redis (simulate TTL expiry)
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    lock_key = f"room:{room}:lock:geometry:test_sphere"
    r.delete(lock_key)

    # Operation should now proceed (no lock held)
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test_sphere",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers=conn.headers,
    )

    # With @check_lock, this should succeed since no lock exists
    assert response.status_code == 200


def test_check_lock_multiple_operations_same_lock(room_with_lock):
    """Test that multiple operations can be performed by lock holder."""
    server, room, session_id, auth_headers, _ = room_with_lock

    # Create first geometry
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "sphere1",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Create second geometry with same lock
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "sphere2",
            "type": "Sphere",
            "data": {"position": [[5, 5, 5]], "radius": [2.0]},
        },
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Verify both geometries exist
    response = requests.get(
        f"{server}/api/rooms/{room}/geometries", headers=auth_headers
    )
    assert response.status_code == 200
    geom_list = response.json()["geometries"]
    assert "sphere1" in geom_list
    assert "sphere2" in geom_list


def test_check_lock_with_missing_jwt(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that decorator rejects requests without JWT token."""
    room = "test-missing-jwt"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create room first
    create_and_join_room(server, room, auth_headers)

    # Try to create geometry without authentication
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={"X-Session-ID": "some-session-id"},
    )

    # Should fail at JWT authentication step (before session validation)
    assert response.status_code in [401, 403]


def test_check_lock_with_invalid_jwt(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that decorator rejects requests with invalid JWT token."""
    room = "test-invalid-jwt"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create room first
    create_and_join_room(server, room, auth_headers)

    # Try to create geometry with invalid JWT
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={
            "Authorization": "Bearer invalid.jwt.token",
            "X-Session-ID": "some-session-id",
        },
    )

    # Should fail at JWT authentication step
    assert response.status_code == 401


def test_check_lock_global_lock_blocks_all(server, connect_room):
    """Test that a global lock blocks all operations by other sessions."""
    room = "test-global-lock"

    # Session 1: Create room, join and acquire GLOBAL lock
    conn1 = connect_room(room, user="test-user")

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/global/acquire",
        json={"msg": "global lock"},
        headers=conn1.headers,
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Session 2: Another user joins
    conn2 = connect_room(room, user="other-user")

    # Session 2 tries to create geometry - should be blocked by global lock
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers=conn2.headers,
    )

    assert response.status_code == 423
    data = response.json()
    assert (
        "globally locked" in data["error"].lower() or "global" in data["error"].lower()
    )


def test_check_lock_global_lock_allows_same_session(server, connect_room):
    """Test that global lock holder can perform operations."""
    room = "test-global-lock-same-session"

    # Create room, join and acquire GLOBAL lock
    conn = connect_room(room, user="test-user")

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/global/acquire",
        json={"msg": "global lock"},
        headers=conn.headers,
    )
    assert response.status_code == 200

    # Same session can still perform operations
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers=conn.headers,
    )

    assert response.status_code == 200


def test_check_lock_per_geometry_lock(server, connect_room):
    """Test that per-geometry locks only block that specific geometry."""
    room = "test-per-geometry-lock"

    # Session 1: Create room and acquire lock on geometry:sphere
    conn1 = connect_room(room, user="user1")

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/geometry:sphere/acquire",
        json={"msg": "locking sphere"},
        headers=conn1.headers,
    )
    assert response.status_code == 200

    # Session 2: Another user joins
    conn2 = connect_room(room, user="user2")

    # Session 2 can create a DIFFERENT geometry (cube) - not blocked by sphere lock
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "cube",  # Different key - not locked
            "type": "Sphere",  # Using Sphere since Box requires different data format
            "data": {"position": [[0, 0, 0]], "radius": [1.0]},
        },
        headers=conn2.headers,
    )
    assert response.status_code == 200

    # But session 2 cannot create geometry with key "sphere" (locked by session 1)
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "sphere",  # Locked key
            "type": "Sphere",
            "data": {"position": [[1, 1, 1]], "radius": [0.5]},
        },
        headers=conn2.headers,
    )
    assert response.status_code == 423
