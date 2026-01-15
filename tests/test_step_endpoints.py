"""Tests for step/frame REST endpoints with lock-based access control.

The @check_lock decorator allows operations to proceed when no lock exists (FIFO ordering).
Only blocks when ANOTHER session holds a lock on the target.
"""

import pytest
import redis
import requests

from zndraw.app.constants import SocketEvents


@pytest.fixture
def room_with_step_lock(server, connect_room):
    """Create a room, join it, acquire step lock, and yield connection info.

    Uses connect_room to keep socket alive during test (prevents lock cleanup).
    """
    conn = connect_room("test-step-lock-room")

    # Acquire lock for step
    response = requests.post(
        f"{server}/api/rooms/{conn.room_id}/locks/step/acquire",
        json={"msg": "testing step lock"},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    lock_data = response.json()
    assert lock_data["success"] is True
    lock_token = lock_data["lockToken"]

    yield server, conn.room_id, conn.session_id, conn.headers, lock_token


def test_get_step_without_auth(server, get_jwt_auth_headers, create_and_join_room):
    """Test that GET /step returns data for existing room."""
    room = "test-get-step"
    auth_headers = get_jwt_auth_headers(server)

    # Create room first
    create_and_join_room(server, room, auth_headers)

    # Get step with authentication
    response = requests.get(
        f"{server}/api/rooms/{room}/step", headers=auth_headers, timeout=10
    )
    assert response.status_code == 200
    data = response.json()
    assert "step" in data
    assert "totalFrames" in data
    # Step might be None for empty room
    assert data["step"] is None or isinstance(data["step"], int)
    assert isinstance(data["totalFrames"], int)


def test_put_step_without_lock(server, get_jwt_auth_headers, create_and_join_room):
    """Test that PUT /step succeeds when no lock is held (FIFO ordering).

    With @check_lock, operations proceed if no lock exists.
    """
    room = "test-put-step-no-lock"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create and join room
    session_id = create_and_join_room(server, room, auth_headers)

    # Set step without acquiring lock first - should succeed with @check_lock
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 42},
        headers={**auth_headers, "X-Session-ID": session_id},
        timeout=10,
    )

    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert data["step"] == 42


def test_put_step_with_lock(room_with_step_lock):
    """Test that PUT /step succeeds with valid lock."""
    server, room, session_id, auth_headers, lock_token = room_with_step_lock

    # Set step with valid lock
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 42},
        headers={**auth_headers, "X-Session-ID": session_id},
        timeout=10,
    )

    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert data["step"] == 42

    # Verify step was updated
    response = requests.get(
        f"{server}/api/rooms/{room}/step", headers=auth_headers, timeout=10
    )
    assert response.status_code == 200
    data = response.json()
    assert data["step"] == 42


def test_put_step_negative_value(server, get_jwt_auth_headers, create_and_join_room):
    """Test that PUT /step rejects negative step values."""
    room = "test-step-negative"
    auth_headers = get_jwt_auth_headers(server, "test-user")
    session_id = create_and_join_room(server, room, auth_headers)

    # Try to set negative step
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": -1},
        headers={**auth_headers, "X-Session-ID": session_id},
        timeout=10,
    )

    assert response.status_code == 400
    data = response.json()
    assert "non-negative" in data["error"]


def test_put_step_invalid_type(server, get_jwt_auth_headers, create_and_join_room):
    """Test that PUT /step rejects non-integer step values."""
    room = "test-step-invalid-type"
    auth_headers = get_jwt_auth_headers(server, "test-user")
    session_id = create_and_join_room(server, room, auth_headers)

    # Try to set step with invalid type (string)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": "not-a-number"},
        headers={**auth_headers, "X-Session-ID": session_id},
        timeout=10,
    )

    assert response.status_code == 400
    data = response.json()
    assert "Invalid step value" in data["error"]


def test_put_step_missing_session_id(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that PUT /step rejects requests without X-Session-ID header."""
    room = "test-missing-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create room
    create_and_join_room(server, room, auth_headers)

    # Try to set step without session ID header
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers=auth_headers,  # Missing X-Session-ID
        timeout=10,
    )

    assert response.status_code == 400
    data = response.json()
    assert "X-Session-ID" in data["error"]


def test_put_step_invalid_session_id(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that PUT /step rejects requests with invalid session ID."""
    room = "test-invalid-session"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create room
    create_and_join_room(server, room, auth_headers)

    # Try to set step with invalid session ID
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers={**auth_headers, "X-Session-ID": "invalid-session-id-12345"},
        timeout=10,
    )

    assert response.status_code == 401
    data = response.json()
    assert "Invalid or expired session" in data["error"]


def test_put_step_continuous_updates(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test continuous step updates without lock (FIFO ordering)."""
    room = "test-step-continuous"
    auth_headers = get_jwt_auth_headers(server, "test-user")
    session_id = create_and_join_room(server, room, auth_headers)

    # Simulate continuous updates (like dragging a slider)
    for step in [0, 10, 20, 30, 40]:
        response = requests.put(
            f"{server}/api/rooms/{room}/step",
            json={"step": step},
            headers={**auth_headers, "X-Session-ID": session_id},
            timeout=10,
        )
        assert response.status_code == 200
        assert response.json()["step"] == step

    # Verify final step
    response = requests.get(
        f"{server}/api/rooms/{room}/step", headers=auth_headers, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["step"] == 40


def test_put_step_after_lock_expiry(server, get_jwt_auth_headers, create_and_join_room):
    """Test that PUT /step succeeds after lock expires.

    With @check_lock, if no lock exists, operations proceed (FIFO ordering).
    """
    room = "test-lock-expiry"
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Create and join room
    session_id = create_and_join_room(server, room, auth_headers)

    # Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "testing expiry"},
        headers={**auth_headers, "X-Session-ID": session_id},
        timeout=10,
    )
    assert response.status_code == 200

    # Manually expire the lock in Redis (simulate TTL expiry)
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    lock_key = f"room:{room}:lock:step"
    r.delete(lock_key)

    # Set step after lock expired - should succeed with @check_lock
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers={**auth_headers, "X-Session-ID": session_id},
        timeout=10,
    )

    assert response.status_code == 200
    assert response.json()["step"] == 10


def test_put_step_session_user_mismatch(
    server, get_jwt_auth_headers, create_and_join_room
):
    """Test that PUT /step rejects requests where session doesn't match JWT user."""
    room = "test-session-mismatch"

    # User 1 creates room and joins
    user1_headers = get_jwt_auth_headers(server, "user1")
    session_id = create_and_join_room(server, room, user1_headers)

    # User 2 gets auth token
    user2_headers = get_jwt_auth_headers(server, "user2")

    # Try to use user1's session with user2's JWT token
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 10},
        headers={**user2_headers, "X-Session-ID": session_id},
        timeout=10,
    )

    assert response.status_code == 403
    data = response.json()
    assert "Session/user mismatch" in data["error"]


def test_put_step_blocked_by_other_session_lock(server, connect_room):
    """Test that different session cannot update step when another holds the lock."""
    # Session 1: Create room, join and acquire lock
    conn1 = connect_room("test-different-session")

    response = requests.post(
        f"{server}/api/rooms/{conn1.room_id}/locks/step/acquire",
        json={"msg": "session 1 lock"},
        headers=conn1.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Session 2: Same user joins in a different session (simulating another browser tab)
    conn2 = connect_room("test-different-session")

    # Session 2 tries to set step - should be blocked by session 1's lock
    response = requests.put(
        f"{server}/api/rooms/{conn1.room_id}/step",
        json={"step": 10},
        headers=conn2.headers,
        timeout=10,
    )

    assert response.status_code == 423
    data = response.json()
    assert "locked" in data["error"].lower() or "test-user" in data["error"]


def test_atomic_pattern_acquire_set_release(server, connect_room):
    """Test atomic update pattern: acquire lock -> set step -> release lock."""
    conn = connect_room("test-atomic-pattern")

    # 1. Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{conn.room_id}/locks/step/acquire",
        json={"msg": "atomic update"},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # 2. Set step
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/step",
        json={"step": 99},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["step"] == 99

    # 3. Release lock
    response = requests.post(
        f"{server}/api/rooms/{conn.room_id}/locks/step/release",
        json={"lockToken": lock_token},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200

    # Verify step was set
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/step", headers=conn.headers, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["step"] == 99


def test_step_lock_independence_from_trajectory_lock(server, connect_room):
    """Test that step operations are independent of trajectory:meta lock.

    With @check_lock, step operations check for 'step' lock, not trajectory:meta.
    """
    conn = connect_room("test-lock-independence")

    # Acquire trajectory:meta lock (but NOT step lock)
    response = requests.post(
        f"{server}/api/rooms/{conn.room_id}/locks/trajectory:meta/acquire",
        json={"msg": "trajectory lock"},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200

    # Setting step should succeed (no step lock held, FIFO ordering)
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/step",
        json={"step": 50},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["step"] == 50


def test_step_update_emits_frame_update_event(server, connect_room):
    """Test that updating step emits frame_update socket event to other clients."""
    import time

    room = "test-frame-update-event"

    # Setup user 1 - create room and join (socket stays connected for lock)
    conn1 = connect_room(room, user="user1")

    # Setup user 2 (observer) - join existing room via socket
    conn2 = connect_room(room, user="user2")

    # User 2 tracks frame_update events
    received_events = []

    def on_frame_update(data):
        received_events.append(data)

    conn2.sio.on(SocketEvents.FRAME_UPDATE, on_frame_update)
    time.sleep(0.3)

    # User 1 updates step (no lock needed with @check_lock)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 99},
        headers=conn1.headers,
        timeout=10,
    )
    assert response.status_code == 200

    # Wait for socket event
    time.sleep(0.5)

    # Verify user 2 received the frame_update event
    assert len(received_events) == 1
    event = received_events[0]
    assert event["frame"] == 99


def test_put_step_out_of_bounds(server, s22_xyz, connect_room):
    """Test that PUT /step rejects step values beyond frame count."""
    import ase
    import ase.io

    import zndraw

    room = "test-step-out-of-bounds"

    # Create room with known number of frames via ZnDraw (handles create+join)
    vis = zndraw.ZnDraw(url=server, room=room, user="test-user")
    structures = ase.io.read(s22_xyz, index=":")  # Read all structures (list of Atoms)
    for atoms in structures:
        if isinstance(atoms, ase.Atoms):  # Type guard for type checker
            vis.append(atoms)  # s22 has 22 structures (indices 0-21)

    # Join room with socket that stays connected (for lock)
    conn = connect_room(room)

    # Test step at upper boundary (should succeed)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 21},  # Last valid index
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["step"] == 21

    # Test step exactly at frame count (should fail)
    response = requests.put(
        f"{server}/api/rooms/{room}/step",
        json={"step": 22},  # One beyond last index
        headers=conn.headers,
        timeout=10,
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
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 400
    data = response.json()
    assert "out of range" in data["error"]


def test_put_step_empty_room(server, connect_room):
    """Test that PUT /step in empty room (0 frames) allows any step value (no validation)."""
    conn = connect_room("test-step-empty-room")

    # Verify room is empty
    response = requests.get(
        f"{server}/api/rooms/{conn.room_id}/step", headers=conn.headers, timeout=10
    )
    assert response.status_code == 200
    assert response.json()["totalFrames"] == 0

    # Set step 0 in empty room (should succeed - no validation when empty)
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/step",
        json={"step": 0},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["step"] == 0

    # Should also allow arbitrary steps in empty room
    response = requests.put(
        f"{server}/api/rooms/{conn.room_id}/step",
        json={"step": 999},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["step"] == 999


def test_acquire_lock_idempotent_same_session(server, connect_room):
    """Test that acquiring a lock twice from the same session returns success (idempotent)."""
    conn = connect_room("test-idempotent-acquire")

    # First acquire
    response1 = requests.post(
        f"{server}/api/rooms/{conn.room_id}/locks/step/acquire",
        json={"msg": "first acquire"},
        headers=conn.headers,
        timeout=10,
    )
    assert response1.status_code == 200
    data1 = response1.json()
    assert data1["success"] is True
    assert data1.get("refreshed", False) is False  # New lock
    token1 = data1["lockToken"]

    # Second acquire (same session) - should succeed and return same token
    response2 = requests.post(
        f"{server}/api/rooms/{conn.room_id}/locks/step/acquire",
        json={"msg": "second acquire"},
        headers=conn.headers,
        timeout=10,
    )
    assert response2.status_code == 200
    data2 = response2.json()
    assert data2["success"] is True
    assert data2["refreshed"] is True  # Refreshed existing lock
    assert data2["lockToken"] == token1  # Same token returned


def test_acquire_lock_different_session_same_user_fails(server, connect_room):
    """Test that different session from same user cannot acquire lock held by another session."""
    # Session 1: Create room, join and acquire lock
    conn1 = connect_room("test-diff-session-acquire")

    response = requests.post(
        f"{server}/api/rooms/{conn1.room_id}/locks/step/acquire",
        json={"msg": "session 1 lock"},
        headers=conn1.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Session 2: Same user joins in a different session (simulating another browser tab)
    conn2 = connect_room("test-diff-session-acquire")

    # Session 2 tries to acquire lock (should fail - different session holds it)
    response = requests.post(
        f"{server}/api/rooms/{conn1.room_id}/locks/step/acquire",
        json={"msg": "session 2 trying"},
        headers=conn2.headers,
        timeout=10,
    )
    assert response.status_code == 423  # Locked
    data = response.json()
    assert data["success"] is False
    assert "Lock already held" in data["error"]
