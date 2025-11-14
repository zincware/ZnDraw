"""Tests for the @requires_lock decorator."""
import pytest
import redis
import requests

from conftest import get_jwt_auth_headers


@pytest.fixture
def room_with_lock(server):
    """Join a room, acquire a lock, and return server, room, session_id, auth_headers, lock_token."""
    room = "test-lock-room"

    # Get auth headers
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join the room to get session ID
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    join_data = response.json()
    session_id = join_data["sessionId"]

    # Acquire lock for trajectory:meta
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "testing lock"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_data = response.json()
    assert lock_data["success"] is True
    lock_token = lock_data["lockToken"]

    return server, room, session_id, auth_headers, lock_token


def test_requires_lock_missing_session_id(server):
    """Test that decorator rejects requests without X-Session-ID header."""
    room = "test-missing-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Try to create geometry without session ID header
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers=auth_headers,  # Missing X-Session-ID
    )

    assert response.status_code == 400
    data = response.json()
    assert "X-Session-ID" in data["error"]


def test_requires_lock_invalid_session_id(server):
    """Test that decorator rejects requests with invalid session ID."""
    room = "test-invalid-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Try to create geometry with invalid session ID
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={**auth_headers, "X-Session-ID": "invalid-session-id-12345"},
    )

    assert response.status_code == 401
    data = response.json()
    assert "Invalid or expired session" in data["error"]


def test_requires_lock_missing_lock_token(server):
    """Test that decorator rejects requests when lock is not held."""
    room = "test-missing-lock"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join room to get valid session ID
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    session_id = response.json()["sessionId"]

    # Try to create geometry without acquiring lock first
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 423  # Locked - lock not held
    data = response.json()
    assert "Lock not held" in data["error"]


def test_requires_lock_session_user_mismatch(server):
    """Test that decorator rejects requests where session doesn't match JWT user."""

    room = "test-session-mismatch"

    # User 1 joins and gets session
    user1_headers = get_jwt_auth_headers(server, "user1")
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=user1_headers
    )
    assert response.status_code == 200
    user1_session_id = response.json()["sessionId"]

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


def test_requires_lock_valid_session_and_token(room_with_lock):
    """Test that decorator allows requests with valid session and lock token."""
    server, room, session_id, auth_headers, lock_token = room_with_lock

    # Create geometry with valid session and lock
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "test_sphere",
            "type": "Sphere",
            "data": {"position": [[1.0, 2.0, 3.0]], "color": "#FF0000", "radius": 1.0},
        },
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    # Verify geometry was created
    response = requests.get(f"{server}/api/rooms/{room}/geometries/test_sphere")
    assert response.status_code == 200
    geom_data = response.json()
    assert geom_data["key"] == "test_sphere"
    assert geom_data["geometry"]["type"] == "Sphere"
    assert geom_data["geometry"]["data"]["position"] == [[1.0, 2.0, 3.0]]


def test_requires_lock_different_session_same_user(server):
    """Test that decorator rejects when a different session (same user) tries to use the lock."""
    room = "test-different-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Session 1: Join and acquire lock
    response1 = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response1.status_code == 200
    session1_id = response1.json()["sessionId"]

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "session 1 lock"},
        headers={**auth_headers, "X-Session-ID": session1_id},
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Session 2: Same user joins in a different session (simulating another browser tab)
    response2 = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response2.status_code == 200
    session2_id = response2.json()["sessionId"]

    # Session 2 tries to create geometry (should fail - session doesn't hold lock)
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={**auth_headers, "X-Session-ID": session2_id},
    )

    assert response.status_code == 403
    data = response.json()
    assert "Session does not hold the lock" in data["error"]


def test_requires_lock_after_lock_expiry(server):
    """Test that decorator rejects after lock expires."""

    room = "test-lock-expiry"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join room
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    session_id = response.json()["sessionId"]

    # Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "testing expiry"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # Manually expire the lock in Redis (simulate TTL expiry)
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    lock_key = f"room:{room}:lock:trajectory:meta"
    r.delete(lock_key)

    # Try to create geometry after lock expired
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 423  # Locked - lock not held
    data = response.json()
    assert "Lock not held" in data["error"]


def test_requires_lock_multiple_operations_same_lock(room_with_lock):
    """Test that multiple operations can be performed with the same lock."""
    server, room, session_id, auth_headers, lock_token = room_with_lock

    # Create first geometry
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "sphere1",
            "type": "Sphere",
            "data": {"position": [[0, 0, 0]], "radius": 1.0},
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
            "data": {"position": [[5, 5, 5]], "radius": 2.0},
        },
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Verify both geometries exist
    response = requests.get(f"{server}/api/rooms/{room}/geometries", headers=auth_headers)
    assert response.status_code == 200
    geom_list = response.json()["geometries"]
    assert "sphere1" in geom_list
    assert "sphere2" in geom_list


def test_requires_lock_with_missing_jwt(server):
    """Test that decorator rejects requests without JWT token."""
    room = "test-missing-jwt"

    # Try to create geometry without authentication
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={"key": "test", "type": "Sphere", "data": {"position": [[0, 0, 0]]}},
        headers={"X-Session-ID": "some-session-id"},
    )

    # Should fail at JWT authentication step (before session validation)
    assert response.status_code in [401, 403]


def test_requires_lock_with_invalid_jwt(server):
    """Test that decorator rejects requests with invalid JWT token."""
    room = "test-invalid-jwt"

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
