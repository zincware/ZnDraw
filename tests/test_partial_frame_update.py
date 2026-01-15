"""Tests for the partial frame update endpoint.

This endpoint allows updating specific keys (e.g., arrays.positions) within a frame
without having to send the entire frame data.

With @check_lock, operations proceed if no lock exists (FIFO ordering).
Only blocks when ANOTHER session holds a lock on the target.
"""

import msgpack
import numpy as np
import pytest
import requests

from zndraw import ZnDraw


@pytest.fixture
def room_with_frames_and_lock(server, s22, connect_room):
    """Create a room with frames and acquire a lock for editing.

    Uses connect_room to keep socket alive during test (prevents lock cleanup).
    """
    room = "partial-update-test"

    # Create room and add frames via ZnDraw
    vis = ZnDraw(url=server, room=room, user="test-user")
    vis.extend(s22[:3])  # Add first 3 frames

    # Join room with socket that stays connected (for lock)
    conn = connect_room(room)

    # Acquire lock for trajectory:meta
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "testing partial update"},
        headers=conn.headers,
        timeout=10,
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    yield server, room, conn.session_id, conn.headers, vis


@pytest.fixture
def room_with_frames_no_lock(server, s22, connect_room):
    """Create a room with frames but NO lock (for testing FIFO behavior).

    Uses connect_room to keep socket alive during test.
    """
    room = "partial-update-no-lock-test"

    # Create room and add frames via ZnDraw
    vis = ZnDraw(url=server, room=room, user="test-user")
    vis.extend(s22[:3])  # Add first 3 frames

    # Join room with socket that stays connected
    conn = connect_room(room)

    yield server, room, conn.session_id, conn.headers, vis


def _get_numpy_dtype_string(arr: np.ndarray) -> str:
    """Get numpy dtype string in the format msgpack-numpy expects."""
    dtype = arr.dtype
    # Map common dtypes to msgpack-numpy format
    dtype_map = {
        np.float64: "<f8",
        np.float32: "<f4",
        np.int64: "<i8",
        np.int32: "<i4",
        np.int16: "<i2",
        np.int8: "|i1",
        np.uint64: "<u8",
        np.uint32: "<u4",
        np.uint16: "<u2",
        np.uint8: "|u1",
        np.bool_: "|b1",
    }
    return dtype_map.get(dtype.type, str(dtype))


def _encode_numpy_for_msgpack(data: dict) -> bytes:
    """Encode a dictionary with numpy arrays to msgpack format.

    This mimics how the frontend encodes data using msgpack-numpy.
    The frontend sends string keys like 'nd', 'type', etc. which the
    backend's decode_with_string_keys hook converts to byte keys.
    """

    def encode_value(obj):
        if isinstance(obj, np.ndarray):
            # Use msgpack-numpy format with string keys (as frontend sends)
            return {
                "nd": True,
                "type": _get_numpy_dtype_string(obj),
                "kind": "",
                "shape": list(obj.shape),
                "data": obj.tobytes(),
            }
        if isinstance(obj, (np.floating, np.integer)):
            # Scalar numpy types - convert to ndarray
            arr = np.asarray(obj)
            return {
                "nd": True,
                "type": _get_numpy_dtype_string(arr),
                "kind": "",
                "shape": list(arr.shape),
                "data": arr.tobytes(),
            }
        return obj

    encoded = {k: encode_value(v) for k, v in data.items()}
    return msgpack.packb(encoded)


def test_partial_update_positions(room_with_frames_and_lock):
    """Test that partial update successfully modifies positions."""
    server, room, session_id, auth_headers, vis = room_with_frames_and_lock

    # Get original positions for frame 0
    original_atoms = vis[0]
    original_positions = original_atoms.positions.copy()

    # Create new positions (shift all atoms by [1, 2, 3])
    new_positions = original_positions + np.array([1.0, 2.0, 3.0])
    new_positions = new_positions.astype(np.float64)

    # Encode the update
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    # Send partial update
    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 200
    result = response.json()
    assert result["success"] is True
    assert result["frame_id"] == 0
    assert "arrays.positions" in result["updated_keys"]

    # Verify the positions were actually updated
    # Need to reconnect to get fresh data
    vis2 = ZnDraw(url=server, room=room, user="verify-user")
    updated_atoms = vis2[0]

    np.testing.assert_array_almost_equal(
        updated_atoms.positions,
        new_positions,
        decimal=5,
        err_msg="Positions were not updated correctly",
    )


def test_partial_update_without_lock_proceeds(room_with_frames_no_lock):
    """Test that partial update succeeds without lock (FIFO ordering).

    With @check_lock, operations proceed if no lock exists.
    """
    server, room, session_id, auth_headers, vis = room_with_frames_no_lock

    # Get original positions for frame 0
    original_atoms = vis[0]
    original_positions = original_atoms.positions.copy()

    # Create new positions
    new_positions = original_positions + np.array([1.0, 2.0, 3.0])
    new_positions = new_positions.astype(np.float64)

    # Encode the update
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    # Send partial update WITHOUT lock - should succeed with @check_lock
    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 200
    result = response.json()
    assert result["success"] is True


def test_partial_update_info(room_with_frames_and_lock):
    """Test that partial update successfully modifies info data."""
    server, room, session_id, auth_headers, vis = room_with_frames_and_lock

    # Create new info data
    energy_value = np.float64(42.5)

    # Encode the update
    update_data = _encode_numpy_for_msgpack({"info.test_energy": energy_value})

    # Send partial update
    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 200
    result = response.json()
    assert result["success"] is True
    assert "info.test_energy" in result["updated_keys"]

    # Verify the info was actually updated
    vis2 = ZnDraw(url=server, room=room, user="verify-user")
    updated_atoms = vis2[0]

    assert "test_energy" in updated_atoms.info
    np.testing.assert_almost_equal(
        updated_atoms.info["test_energy"],
        42.5,
        decimal=5,
    )


def test_partial_update_multiple_keys(room_with_frames_and_lock):
    """Test that partial update can modify multiple keys at once."""
    server, room, session_id, auth_headers, vis = room_with_frames_and_lock

    original_atoms = vis[0]
    num_atoms = len(original_atoms)

    # Create new data for multiple keys
    new_positions = original_atoms.positions + np.array([0.5, 0.5, 0.5])
    new_positions = new_positions.astype(np.float64)
    custom_array = np.ones(num_atoms, dtype=np.float64) * 3.14

    # Encode multiple updates
    update_data = _encode_numpy_for_msgpack(
        {
            "arrays.positions": new_positions,
            "arrays.custom_data": custom_array,
        }
    )

    # Send partial update
    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 200
    result = response.json()
    assert result["success"] is True
    assert set(result["updated_keys"]) == {"arrays.positions", "arrays.custom_data"}


def test_partial_update_invalid_frame_index(room_with_frames_and_lock):
    """Test that partial update returns 404 for non-existing frame."""
    server, room, session_id, auth_headers, _ = room_with_frames_and_lock

    new_positions = np.array([[0, 0, 0]], dtype=np.float64)
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    # Try to update frame 999 which doesn't exist
    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/999/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 404
    result = response.json()
    assert "Invalid frame index" in result["error"]


def test_partial_update_negative_frame_index(room_with_frames_and_lock):
    """Test that partial update rejects negative frame index.

    Flask's <int:frame_id> route pattern doesn't match negative numbers,
    so this returns 405 Method Not Allowed (no matching route).
    """
    server, room, session_id, auth_headers, _ = room_with_frames_and_lock

    new_positions = np.array([[0, 0, 0]], dtype=np.float64)
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    # Try to update frame -1
    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/-1/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    # Flask's <int:frame_id> doesn't match negative numbers, returns 405
    assert response.status_code == 405


def test_partial_update_empty_body(room_with_frames_and_lock):
    """Test that partial update returns 400 for empty update body."""
    server, room, session_id, auth_headers, _ = room_with_frames_and_lock

    # Send empty dict
    update_data = msgpack.packb({})

    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 400
    result = response.json()
    assert "No keys provided" in result["error"]


def test_partial_update_blocked_by_other_session_lock(server, s22, connect_room):
    """Test that partial update is blocked when another session holds the lock."""
    room = "partial-update-blocked"

    # Create room and add frames
    vis = ZnDraw(url=server, room=room, user="test-user")
    vis.extend(s22[:1])

    # Session 1 joins and acquires lock
    conn1 = connect_room(room, user="user1")
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "user1 lock"},
        headers=conn1.headers,
        timeout=10,
    )
    assert response.status_code == 200

    # Session 2 joins (different user)
    conn2 = connect_room(room, user="user2")

    # Session 2 tries to update - should be blocked
    new_positions = np.array([[0, 0, 0]], dtype=np.float64)
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **conn2.headers,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 423  # Locked by another session
    result = response.json()
    assert "locked" in result["error"].lower() or "user1" in result["error"]


def test_partial_update_requires_session_id(server, s22, get_jwt_auth_headers):
    """Test that partial update requires X-Session-ID header."""
    room = "partial-update-no-session"

    # Create room and add frames
    vis = ZnDraw(url=server, room=room, user="test-user")
    vis.extend(s22[:1])

    # Get auth headers but don't include session ID
    auth_headers = get_jwt_auth_headers(server, "test-user")

    new_positions = np.array([[0, 0, 0]], dtype=np.float64)
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **auth_headers,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 400
    result = response.json()
    assert "X-Session-ID" in result["error"]


def test_partial_update_different_frame(room_with_frames_and_lock):
    """Test that partial update works on frames other than frame 0."""
    server, room, session_id, auth_headers, vis = room_with_frames_and_lock

    # Get original positions for frame 1
    original_atoms = vis[1]
    original_positions = original_atoms.positions.copy()

    # Create new positions
    new_positions = original_positions * 2.0
    new_positions = new_positions.astype(np.float64)

    # Encode the update
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    # Send partial update to frame 1
    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/1/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 200
    result = response.json()
    assert result["success"] is True
    assert result["frame_id"] == 1

    # Verify frame 1 was updated
    vis2 = ZnDraw(url=server, room=room, user="verify-user")
    updated_atoms = vis2[1]

    np.testing.assert_array_almost_equal(
        updated_atoms.positions,
        new_positions,
        decimal=5,
    )

    # Verify frame 0 was NOT modified
    frame0_atoms = vis2[0]
    original_frame0 = vis[0]
    np.testing.assert_array_almost_equal(
        frame0_atoms.positions,
        original_frame0.positions,
        decimal=5,
        err_msg="Frame 0 should not have been modified",
    )


def test_partial_update_preserves_other_data(room_with_frames_and_lock):
    """Test that partial update preserves data not being updated."""
    server, room, session_id, auth_headers, vis = room_with_frames_and_lock

    # Get original data for frame 0
    original_atoms = vis[0]
    original_numbers = original_atoms.numbers.copy()
    original_cell = original_atoms.cell.array.copy()

    # Update only positions
    new_positions = original_atoms.positions + np.array([5.0, 5.0, 5.0])
    new_positions = new_positions.astype(np.float64)

    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=update_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 200

    # Verify other data was preserved
    vis2 = ZnDraw(url=server, room=room, user="verify-user")
    updated_atoms = vis2[0]

    # Positions should be updated
    np.testing.assert_array_almost_equal(
        updated_atoms.positions,
        new_positions,
        decimal=5,
    )

    # Numbers should be preserved
    np.testing.assert_array_equal(
        updated_atoms.numbers,
        original_numbers,
        err_msg="Atomic numbers should be preserved",
    )

    # Cell should be preserved
    np.testing.assert_array_almost_equal(
        updated_atoms.cell.array,
        original_cell,
        decimal=5,
        err_msg="Cell should be preserved",
    )


def test_partial_update_empty_room(server, connect_room):
    """Test that partial update returns 404 for empty room (no frames)."""
    conn = connect_room("partial-update-empty-room")

    # Try to update frame 0 in empty room (no lock needed, FIFO)
    new_positions = np.array([[0, 0, 0]], dtype=np.float64)
    update_data = _encode_numpy_for_msgpack({"arrays.positions": new_positions})

    response = requests.patch(
        f"{server}/api/rooms/{conn.room_id}/frames/0/partial",
        data=update_data,
        headers={
            **conn.headers,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 404
    result = response.json()
    assert "No frames found" in result["error"]


def test_partial_update_malformed_msgpack(room_with_frames_and_lock):
    """Test that partial update returns 400 for malformed msgpack data."""
    server, room, session_id, auth_headers, _ = room_with_frames_and_lock

    # Send invalid msgpack bytes (random garbage data)
    malformed_data = b"\x91\x92\x93\xff\xfe\xfd\xfc"

    response = requests.patch(
        f"{server}/api/rooms/{room}/frames/0/partial",
        data=malformed_data,
        headers={
            **auth_headers,
            "X-Session-ID": session_id,
            "Content-Type": "application/msgpack",
        },
        timeout=10,
    )

    assert response.status_code == 400
    result = response.json()
    assert "error" in result
