"""Tests for step/frame REST endpoints with lock-based access control."""
import pytest
import redis
import requests

from conftest import get_jwt_auth_headers


@pytest.fixture
def room_with_step_lock(server):
    """Join a room, acquire step lock, and return server, room, session_id, auth_headers, lock_token."""
    room = "test-step-lock-room"

    # Get auth headers
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join the room to get session ID
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    join_data = response.json()
    session_id = join_data["sessionId"]

    # Acquire lock for step
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "testing step lock"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_data = response.json()
    assert lock_data["success"] is True
    lock_token = lock_data["lockToken"]

    return server, room, session_id, auth_headers, lock_token


def test_get_step_without_auth(server):
    """Test that GET /step requires authentication."""
    room = "test-get-step"
    auth_headers = get_jwt_auth_headers(server)

    # Get step with authentication
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=auth_headers)
    assert response.status_code == 200
    data = response.json()
    assert "step" in data
    assert "totalFrames" in data
    # Step might be None for empty room
    assert data["step"] is None or isinstance(data["step"], int)
    assert isinstance(data["totalFrames"], int)


def test_put_step_without_lock(server):
    """Test that PUT /step fails when lock is not held."""
    room = "test-put-step-no-lock"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join room to get valid session ID
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    session_id = response.json()["sessionId"]

    # Try to set step without acquiring lock first
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 42},
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 423  # Locked - lock not held
    data = response.json()
    assert "Lock not held" in data["error"]


def test_put_step_with_lock(room_with_step_lock):
    """Test that PUT /step succeeds with valid lock."""
    server, room, session_id, auth_headers, lock_token = room_with_step_lock

    # Set step with valid lock
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 42},
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert data["step"] == 42

    # Verify step was updated
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=auth_headers)
    assert response.status_code == 200
    data = response.json()
    assert data["step"] == 42


def test_put_step_negative_value(room_with_step_lock):
    """Test that PUT /step rejects negative step values."""
    server, room, session_id, auth_headers, lock_token = room_with_step_lock

    # Try to set negative step
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": -1},
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 400
    data = response.json()
    assert "non-negative" in data["error"]


def test_put_step_invalid_type(room_with_step_lock):
    """Test that PUT /step rejects non-integer step values."""
    server, room, session_id, auth_headers, lock_token = room_with_step_lock

    # Try to set step with invalid type (string)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": "not-a-number"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 400
    data = response.json()
    assert "Invalid step value" in data["error"]


def test_put_step_missing_session_id(server):
    """Test that PUT /step rejects requests without X-Session-ID header."""
    room = "test-missing-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Try to set step without session ID header
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers=auth_headers,  # Missing X-Session-ID
    )

    assert response.status_code == 400
    data = response.json()
    assert "X-Session-ID" in data["error"]


def test_put_step_invalid_session_id(server):
    """Test that PUT /step rejects requests with invalid session ID."""
    room = "test-invalid-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Try to set step with invalid session ID
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers={**auth_headers, "X-Session-ID": "invalid-session-id-12345"},
    )

    assert response.status_code == 401
    data = response.json()
    assert "Invalid or expired session" in data["error"]


def test_put_step_continuous_updates(room_with_step_lock):
    """Test continuous step updates while holding lock (presentation mode)."""
    server, room, session_id, auth_headers, lock_token = room_with_step_lock

    # Simulate continuous updates (like dragging a slider)
    for step in [0, 10, 20, 30, 40]:
        response = requests.put(
            f"{server}/api/rooms/{room}/step",
            json={"step": step},
            headers={**auth_headers, "X-Session-ID": session_id},
        )
        assert response.status_code == 200
        assert response.json()["step"] == step

    # Verify final step
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=auth_headers)
    assert response.status_code == 200
    assert response.json()["step"] == 40


def test_put_step_after_lock_expiry(server):
    """Test that PUT /step fails after lock expires."""
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
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "testing expiry"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # Manually expire the lock in Redis (simulate TTL expiry)
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    lock_key = f"room:{room}:lock:step"
    r.delete(lock_key)

    # Try to set step after lock expired
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers={**auth_headers, "X-Session-ID": session_id},
    )

    assert response.status_code == 423  # Locked - lock not held
    data = response.json()
    assert "Lock not held" in data["error"]


def test_put_step_session_user_mismatch(server):
    """Test that PUT /step rejects requests where session doesn't match JWT user."""
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
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers={**user2_headers, "X-Session-ID": user1_session_id},
    )

    assert response.status_code == 403
    data = response.json()
    assert "Session/user mismatch" in data["error"]


def test_put_step_different_session_same_user(server):
    """Test that different session from same user can't use another session's lock."""
    room = "test-different-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Session 1: Join and acquire lock
    response1 = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response1.status_code == 200
    session1_id = response1.json()["sessionId"]

    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
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

    # Session 2 tries to set step (should fail - session doesn't hold lock)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers={**auth_headers, "X-Session-ID": session2_id},
    )

    assert response.status_code == 403
    data = response.json()
    assert "Session does not hold the lock" in data["error"]


def test_atomic_pattern_acquire_set_release(server):
    """Test atomic update pattern: acquire lock → set step → release lock."""
    room = "test-atomic-pattern"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join room
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    session_id = response.json()["sessionId"]

    # 1. Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "atomic update"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # 2. Set step
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 99},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["step"] == 99

    # 3. Release lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/release",
        json={"lockToken": lock_token},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Verify step was set
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=auth_headers)
    assert response.status_code == 200
    assert response.json()["step"] == 99


def test_step_lock_independence_from_trajectory_lock(server):
    """Test that step lock is independent from trajectory:meta lock."""
    room = "test-lock-independence"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join room
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    session_id = response.json()["sessionId"]

    # Acquire trajectory:meta lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "trajectory lock"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Try to set step without step lock (should fail)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 50},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 423
    assert "Lock not held" in response.json()["error"]

    # Now acquire step lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "step lock"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Now setting step should work
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 50},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["step"] == 50


def test_step_update_emits_frame_update_event(server):
    """Test that updating step emits frame_update socket event to other clients."""
    import socketio
    import time

    room = "test-frame-update-event"

    # Setup user 1
    user1_headers = get_jwt_auth_headers(server, "user1")
    response = requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=user1_headers)
    assert response.status_code == 200
    user1_session = response.json()["sessionId"]
    user1_token = user1_headers["Authorization"].split(" ")[1]

    # Setup user 2 (observer)
    user2_headers = get_jwt_auth_headers(server, "user2")
    response = requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=user2_headers)
    assert response.status_code == 200
    user2_token = user2_headers["Authorization"].split(" ")[1]

    # User 2 connects via socket and tracks frame_update events
    sio_client = socketio.Client()
    received_events = []

    def on_frame_update(data):
        received_events.append(data)

    sio_client.on("frame_update", on_frame_update)
    sio_client.connect(server, auth={"token": user2_token})
    sio_client.emit("join:room", {"roomId": room})
    time.sleep(0.5)

    # User 1 acquires lock and updates step
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "updating step"},
        headers={**user1_headers, "X-Session-ID": user1_session},
    )
    assert response.status_code == 200

    # User 1 updates step
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 99},
        headers={**user1_headers, "X-Session-ID": user1_session},
    )
    assert response.status_code == 200

    # Wait for socket event
    time.sleep(0.5)

    # Verify user 2 received the frame_update event
    assert len(received_events) == 1
    event = received_events[0]
    assert event["frame"] == 99

    sio_client.disconnect()


def test_put_step_out_of_bounds(server, s22_xyz):
    """Test that PUT /step rejects step values beyond frame count."""
    import ase
    import ase.io
    import zndraw

    room = "test-step-out-of-bounds"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create room with known number of frames
    vis = zndraw.ZnDraw(url=server, room=room, user="test-user")
    structures = ase.io.read(s22_xyz, index=":")  # Read all structures (list of Atoms)
    for atoms in structures:
        if isinstance(atoms, ase.Atoms):  # Type guard for type checker
            vis.append(atoms)  # s22 has 22 structures (indices 0-21)

    # Join room
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    session_id = response.json()["sessionId"]

    # Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "testing out of bounds"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Test step at upper boundary (should succeed)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 21},  # Last valid index
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["step"] == 21

    # Test step exactly at frame count (should fail)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 22},  # One beyond last index
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 400
    data = response.json()
    assert "out of range" in data["error"]
    assert "22" in data["error"]  # Should mention the invalid step
    assert "0-21" in data["error"]  # Should mention valid range

    # Test step way beyond frame count (should fail)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 999999999},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 400
    data = response.json()
    assert "out of range" in data["error"]


def test_put_step_empty_room(server):
    """Test that PUT /step in empty room (0 frames) allows any step value (no validation)."""
    room = "test-step-empty-room"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Join empty room
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    session_id = response.json()["sessionId"]

    # Verify room is empty
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=auth_headers)
    assert response.status_code == 200
    assert response.json()["totalFrames"] == 0

    # Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "testing empty room"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Set step 0 in empty room (should succeed - no validation when empty)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 0},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["step"] == 0

    # Should also allow arbitrary steps in empty room
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 999},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["step"] == 999
