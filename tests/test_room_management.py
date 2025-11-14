"""Tests for room management API endpoints.

This module tests the new room management architecture including:
- Room metadata (description, locked, hidden)
- Default room management
- Room duplication
- Lock enforcement
"""

import pytest
import redis
import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw


def test_list_rooms_includes_metadata(server, s22):
    """Test that GET /api/rooms returns all metadata fields."""
    # Create room with data and set metadata
    vis = ZnDraw(url=server, room="test-room-1", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set("room:test-room-1:description", "Test room description")
    r.set("room:test-room-1:locked", "1")
    r.set("room:test-room-1:hidden", "0")
    r.set("default_room", "test-room-1")

    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200

    rooms = response.json()
    test_room = [room for room in rooms if room["id"] == "test-room-1"][0]

    assert test_room["description"] == "Test room description"
    assert test_room["frameCount"] == 1
    assert test_room["locked"] is True
    assert test_room["hidden"] is False
    assert test_room["isDefault"] is True


def test_list_rooms_without_description(server, s22):
    """Test that rooms without description return null."""
    vis = ZnDraw(url=server, room="test-room-2", user="user1")
    vis.append(s22[0])

    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200

    rooms = response.json()
    test_room = [room for room in rooms if room["id"] == "test-room-2"][0]

    assert test_room["description"] is None
    assert test_room["locked"] is False  # Default
    assert test_room["hidden"] is False  # Default
    assert test_room["isDefault"] is False


def test_list_rooms_metadata_locked(server, s22):
    """Test that GET /api/rooms returns metadataLocked when lock is held."""
    vis = ZnDraw(url=server, room="test-room-metalock", user="user1")
    vis.append(s22[0])

    # Initially, metadata lock should be False
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    test_room = [room for room in rooms if room["id"] == "test-room-metalock"][0]
    assert test_room["metadataLocked"] is False

    # Acquire metadata lock - need to store the lock instance
    lock = vis.get_lock()
    lock.acquire()
    try:
        # Now metadataLocked should be True
        response = requests.get(f"{server}/api/rooms")
        assert response.status_code == 200
        rooms = response.json()
        test_room = [room for room in rooms if room["id"] == "test-room-metalock"][0]
        assert test_room["metadataLocked"] is True
    finally:
        # Release lock - use the same instance
        lock.release()

    # metadataLocked should be False again
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    test_room = [room for room in rooms if room["id"] == "test-room-metalock"][0]
    assert test_room["metadataLocked"] is False


def test_list_rooms_metadata_locked_context_manager(server, s22):
    """Test that GET /api/rooms returns metadataLocked when using with vis.lock."""
    vis = ZnDraw(url=server, room="test-room-metalock-ctx", user="user1")
    vis.append(s22[0])

    # Initially, metadata lock should be False
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    test_room = [room for room in rooms if room["id"] == "test-room-metalock-ctx"][0]
    assert test_room["metadataLocked"] is False

    # Acquire metadata lock using context manager
    with vis.get_lock():
        # Now metadataLocked should be True
        response = requests.get(f"{server}/api/rooms")
        assert response.status_code == 200
        rooms = response.json()
        test_room = [room for room in rooms if room["id"] == "test-room-metalock-ctx"][
            0
        ]
        assert test_room["metadataLocked"] is True

    # metadataLocked should be False again after context exit
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    test_room = [room for room in rooms if room["id"] == "test-room-metalock-ctx"][0]
    assert test_room["metadataLocked"] is False


def test_metadata_lock_refresh_long_operation(server, s22):
    """Test that lock refresh keeps the lock active during long operations."""
    import time

    vis = ZnDraw(url=server, room="test-room-metalock-long", user="user1")
    vis.append(s22[0])

    # Use the default server-controlled TTL (60s with refresh every 30s)
    lock = vis.get_lock(msg="Long running operation")

    # Acquire lock and hold it for a few seconds
    lock.acquire()
    try:
        # Verify lock is held
        response = requests.get(f"{server}/api/rooms")
        rooms = response.json()
        test_room = [room for room in rooms if room["id"] == "test-room-metalock-long"][
            0
        ]
        assert test_room["metadataLocked"] is True

        # Wait 3 seconds - verify lock is still held
        time.sleep(3)

        # Lock should still be held (background refresh thread keeps it alive)
        response = requests.get(f"{server}/api/rooms")
        rooms = response.json()
        test_room = [room for room in rooms if room["id"] == "test-room-metalock-long"][
            0
        ]
        assert test_room["metadataLocked"] is True
    finally:
        lock.release()

    # After release, lock should be gone
    response = requests.get(f"{server}/api/rooms")
    rooms = response.json()
    test_room = [room for room in rooms if room["id"] == "test-room-metalock-long"][0]
    assert test_room["metadataLocked"] is False


def test_metadata_lock_ttl_validation(server, s22):
    """Test that lock TTL is controlled by server configuration."""
    vis = ZnDraw(url=server, room="test-room-ttl-validation", user="user1")
    vis.append(s22[0])

    # TTL is now server-controlled via LockConfig.DEFAULT_TTL (60 seconds)
    # Client cannot specify custom TTL
    with vis.get_lock(msg="Testing server-controlled TTL") as lock:
        # Verify lock was acquired successfully
        # Server will use DEFAULT_TTL=60 and DEFAULT_REFRESH_INTERVAL=30
        assert lock._is_held
        assert lock._ttl == 60  # Server-provided TTL
        assert lock._refresh_interval == 30  # Server-provided refresh interval


def test_get_room_details(server, s22):
    """Test that GET /api/rooms/{room_id} returns detailed metadata."""
    vis = ZnDraw(url=server, room="test-room-3", user="user1")
    vis.extend(s22[:3])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set("room:test-room-3:description", "Detailed room")
    r.set("room:test-room-3:locked", "0")
    r.set("room:test-room-3:hidden", "1")

    response = requests.get(f"{server}/api/rooms/test-room-3", headers=get_jwt_auth_headers(server))
    assert response.status_code == 200

    room = response.json()
    assert room["id"] == "test-room-3"
    assert room["description"] == "Detailed room"
    assert room["frameCount"] == 3
    assert room["locked"] is False
    assert room["hidden"] is True


def test_get_nonexistent_room(server):
    """Test that getting a nonexistent room returns 404."""
    response = requests.get(f"{server}/api/rooms/nonexistent-room", headers=get_jwt_auth_headers(server))
    assert response.status_code == 404


def test_update_room_description(server, s22):
    """Test updating room description."""
    vis = ZnDraw(url=server, room="test-room-4", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)

    # Set description
    response = requests.patch(
        f"{server}/api/rooms/test-room-4", json={"description": "New description"}
    )
    assert response.status_code == 200
    assert r.get("room:test-room-4:description") == "New description"

    # Clear description
    response = requests.patch(
        f"{server}/api/rooms/test-room-4", json={"description": None}
    )
    assert response.status_code == 200
    assert r.get("room:test-room-4:description") is None


def test_update_room_locked_flag(server, s22):
    """Test updating room locked flag."""
    vis = ZnDraw(url=server, room="test-room-5", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)

    # Lock room
    response = requests.patch(f"{server}/api/rooms/test-room-5", json={"locked": True})
    assert response.status_code == 200
    assert r.get("room:test-room-5:locked") == "1"

    # Unlock room
    response = requests.patch(f"{server}/api/rooms/test-room-5", json={"locked": False})
    assert response.status_code == 200
    assert r.get("room:test-room-5:locked") == "0"


def test_update_room_hidden_flag(server, s22):
    """Test updating room hidden flag."""
    vis = ZnDraw(url=server, room="test-room-6", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)

    # Hide room
    response = requests.patch(f"{server}/api/rooms/test-room-6", json={"hidden": True})
    assert response.status_code == 200
    assert r.get("room:test-room-6:hidden") == "1"

    # Unhide room
    response = requests.patch(f"{server}/api/rooms/test-room-6", json={"hidden": False})
    assert response.status_code == 200
    assert r.get("room:test-room-6:hidden") == "0"


def test_update_multiple_fields(server, s22):
    """Test updating multiple room fields at once."""
    vis = ZnDraw(url=server, room="test-room-7", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)

    response = requests.patch(
        f"{server}/api/rooms/test-room-7",
        json={"description": "Multi-field update", "locked": True, "hidden": False},
    )
    assert response.status_code == 200
    assert r.get("room:test-room-7:description") == "Multi-field update"
    assert r.get("room:test-room-7:locked") == "1"
    assert r.get("room:test-room-7:hidden") == "0"


def test_update_nonexistent_room_fails(server):
    """Test that updating a nonexistent room returns 404."""
    response = requests.patch(
        f"{server}/api/rooms/nonexistent-room", json={"description": "Should fail"}
    )
    assert response.status_code == 404


def test_get_default_room(server, s22):
    """Test getting the default room."""
    vis = ZnDraw(url=server, room="default-room", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set("default_room", "default-room")

    response = requests.get(f"{server}/api/rooms/default")
    assert response.status_code == 200

    data = response.json()
    assert data["roomId"] == "default-room"


def test_get_default_room_when_none_set(server):
    """Test getting default room when none is set."""
    response = requests.get(f"{server}/api/rooms/default")
    assert response.status_code == 200

    data = response.json()
    assert data["roomId"] is None


def test_set_default_room(server, s22):
    """Test setting the default room."""
    vis = ZnDraw(url=server, room="new-default-room", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)

    response = requests.put(
        f"{server}/api/rooms/default", json={"roomId": "new-default-room"}
    )
    assert response.status_code == 200
    assert r.get("default_room") == "new-default-room"


def test_unset_default_room(server, s22):
    """Test unsetting the default room."""
    vis = ZnDraw(url=server, room="some-room", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set("default_room", "some-room")

    response = requests.put(f"{server}/api/rooms/default", json={"roomId": None})
    assert response.status_code == 200
    assert r.get("default_room") is None


def test_set_nonexistent_default_room_fails(server):
    """Test that setting a nonexistent room as default fails."""
    response = requests.put(
        f"{server}/api/rooms/default", json={"roomId": "nonexistent-room"}
    )
    assert response.status_code == 404


def test_duplicate_room_basic(server, s22):
    """Test basic room duplication."""
    vis = ZnDraw(url=server, room="source-room", user="user1")
    vis.extend(s22[:2])

    response = requests.post(f"{server}/api/rooms/source-room/duplicate", json={})
    assert response.status_code == 200

    data = response.json()
    assert data["status"] == "ok"
    assert "roomId" in data
    assert data["frameCount"] == 2

    new_room_id = data["roomId"]

    # Verify new room has same frames
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    new_indices = r.zrange(f"room:{new_room_id}:trajectory:indices", 0, -1)
    assert len(new_indices) == 2


def test_duplicate_room_with_custom_id(server, s22):
    """Test room duplication with custom room ID."""
    vis = ZnDraw(url=server, room="source-room-2", user="user1")
    vis.append(s22[0])

    new_room_id = "custom-new-room"
    response = requests.post(
        f"{server}/api/rooms/source-room-2/duplicate", json={"newRoomId": new_room_id}
    )
    assert response.status_code == 200

    data = response.json()
    assert data["roomId"] == new_room_id

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    assert r.exists(f"room:{new_room_id}:trajectory:indices")


def test_duplicate_room_with_description(server, s22):
    """Test room duplication with description."""
    vis = ZnDraw(url=server, room="source-room-3", user="user1")
    vis.append(s22[0])

    response = requests.post(
        f"{server}/api/rooms/source-room-3/duplicate",
        json={"description": "Copied room"},
    )
    assert response.status_code == 200

    data = response.json()
    new_room_id = data["roomId"]

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    assert r.get(f"room:{new_room_id}:description") == "Copied room"


def test_duplicate_room_copies_geometries(server, s22):
    """Test that room duplication copies geometries."""
    vis = ZnDraw(url=server, room="source-room-4", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.hset("room:source-room-4:geometries", "particles", '{"type": "Sphere"}')

    response = requests.post(f"{server}/api/rooms/source-room-4/duplicate", json={})
    assert response.status_code == 200

    new_room_id = response.json()["roomId"]
    geometries = r.hgetall(f"room:{new_room_id}:geometries")
    assert "particles" in geometries


def test_duplicate_room_copies_bookmarks(server, s22):
    """Test that room duplication copies bookmarks."""
    vis = ZnDraw(url=server, room="source-room-5", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.hset("room:source-room-5:bookmarks", "source-room-5:0", "First frame")

    response = requests.post(f"{server}/api/rooms/source-room-5/duplicate", json={})
    assert response.status_code == 200

    new_room_id = response.json()["roomId"]
    bookmarks = r.hgetall(f"room:{new_room_id}:bookmarks")
    assert "source-room-5:0" in bookmarks  # Physical key remains same


def test_duplicate_room_initializes_flags(server, s22):
    """Test that duplicated room has correct initial flags."""
    vis = ZnDraw(url=server, room="source-room-6", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set("room:source-room-6:locked", "1")
    r.set("room:source-room-6:hidden", "1")

    response = requests.post(f"{server}/api/rooms/source-room-6/duplicate", json={})
    assert response.status_code == 200

    new_room_id = response.json()["roomId"]

    # New room should be unlocked and visible
    assert r.get(f"room:{new_room_id}:locked") == "0"
    assert r.get(f"room:{new_room_id}:hidden") == "0"
    assert r.get(f"room:{new_room_id}:current_frame") == "0"


def test_duplicate_nonexistent_room_fails(server):
    """Test that duplicating a nonexistent room fails."""
    response = requests.post(f"{server}/api/rooms/nonexistent-room/duplicate", json={})
    assert response.status_code == 404


def test_duplicate_to_existing_room_fails(server, s22):
    """Test that duplicating to an existing room ID fails."""
    vis1 = ZnDraw(url=server, room="source-room-7", user="user1")
    vis1.append(s22[0])

    vis2 = ZnDraw(url=server, room="existing-room", user="user1")
    vis2.append(s22[0])

    response = requests.post(
        f"{server}/api/rooms/source-room-7/duplicate",
        json={"newRoomId": "existing-room"},
    )
    assert response.status_code == 409


def test_unlocked_room_allows_mutations(server, s22):
    """Test that unlocked rooms allow mutations."""
    vis = ZnDraw(url=server, room="unlocked-room", user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set("room:unlocked-room:locked", "0")

    # Should be able to append another frame
    vis.append(s22[1])
    assert len(vis) == 2


@pytest.mark.parametrize("locked_value", ["1", "0"])
def test_locked_room_allows_reads(server, s22, locked_value):
    """Test that locked rooms still allow read operations."""
    room_id = f"read-test-room-{locked_value}"
    vis = ZnDraw(url=server, room=room_id, user="user1")
    vis.append(s22[0])

    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set(f"room:{room_id}:locked", locked_value)

    # GET requests should work regardless of lock status
    response = requests.get(f"{server}/api/rooms/{room_id}/frames")
    # May fail for other reasons, but not due to lock (should be 200 with data)
    assert response.status_code != 403


def test_lock_with_message(server, s22):
    """Test that lock can be acquired with a custom message."""
    vis = ZnDraw(url=server, room="test-lock-msg", user="user1")
    vis.append(s22[0])

    # Acquire lock with message
    with vis.get_lock(msg="Uploading trajectory data"):
        # Check lock status via REST API
        response = requests.get(
            f"{server}/api/rooms/test-lock-msg/locks/trajectory:meta"
        )
        assert response.status_code == 200
        lock_status = response.json()

        assert lock_status["locked"] is True
        assert lock_status["target"] == "trajectory:meta"
        assert "metadata" in lock_status
        assert lock_status["metadata"]["msg"] == "Uploading trajectory data"
        assert lock_status["metadata"]["userName"] == "user1"
        assert "timestamp" in lock_status["metadata"]

    # After release, lock should be gone
    response = requests.get(f"{server}/api/rooms/test-lock-msg/locks/trajectory:meta")
    assert response.status_code == 200
    lock_status = response.json()
    assert lock_status["locked"] is False


def test_lock_with_custom_metadata(server, s22):
    """Test that lock can be acquired with a message."""
    vis = ZnDraw(url=server, room="test-lock-metadata", user="user1")
    vis.append(s22[0])

    # Acquire lock with message
    with vis.get_lock(msg="Processing simulation step 42/100"):
        # Check lock status
        response = requests.get(
            f"{server}/api/rooms/test-lock-metadata/locks/trajectory:meta"
        )
        assert response.status_code == 200
        lock_status = response.json()

        assert lock_status["locked"] is True
        assert lock_status["metadata"]["msg"] == "Processing simulation step 42/100"

    # After release, lock should be gone
    response = requests.get(
        f"{server}/api/rooms/test-lock-metadata/locks/trajectory:meta"
    )
    lock_status = response.json()
    assert lock_status["locked"] is False


def test_lock_with_message_and_metadata(server, s22):
    """Test that lock can be acquired with a message."""
    vis = ZnDraw(url=server, room="test-lock-both", user="user1")
    vis.append(s22[0])

    # Acquire lock with message
    with vis.get_lock(msg="Processing batch 1/10"):
        # Check lock status
        response = requests.get(
            f"{server}/api/rooms/test-lock-both/locks/trajectory:meta"
        )
        assert response.status_code == 200
        lock_status = response.json()

        assert lock_status["locked"] is True
        assert lock_status["metadata"]["msg"] == "Processing batch 1/10"


def test_lock_reentrant_with_different_messages(server, s22):
    """Test that nested locking is not supported and raises RuntimeError."""
    vis = ZnDraw(url=server, room="test-lock-reentrant", user="user1")
    vis.append(s22[0])

    with vis.get_lock(msg="Outer operation"):
        # Check outer message
        response = requests.get(
            f"{server}/api/rooms/test-lock-reentrant/locks/trajectory:meta"
        )
        lock_status = response.json()
        assert lock_status["metadata"]["msg"] == "Outer operation"

        # Nested locking should fail to acquire (server rejects it)
        with pytest.raises(RuntimeError, match="Failed to acquire lock"):
            with vis.get_lock(msg="Inner operation"):
                pass

    # After release, lock should be gone
    response = requests.get(
        f"{server}/api/rooms/test-lock-reentrant/locks/trajectory:meta"
    )
    lock_status = response.json()
    assert lock_status["locked"] is False


def test_lock_without_metadata_still_works(server, s22):
    """Test that lock without metadata still works (backward compatibility)."""
    vis = ZnDraw(url=server, room="test-lock-no-meta", user="user1")
    vis.append(s22[0])

    # Acquire lock without metadata (old style)
    with vis.get_lock():
        # Check lock is acquired
        response = requests.get(
            f"{server}/api/rooms/test-lock-no-meta/locks/trajectory:meta"
        )
        assert response.status_code == 200
        lock_status = response.json()

        assert lock_status["locked"] is True
        # Metadata should be empty since no message was provided
        assert lock_status["metadata"] == {}

    # After release, lock should be gone
    response = requests.get(
        f"{server}/api/rooms/test-lock-no-meta/locks/trajectory:meta"
    )
    lock_status = response.json()
    assert lock_status["locked"] is False


def test_lock_status_unlocked_target(server):
    """Test getting lock status for an unlocked target."""
    response = requests.get(f"{server}/api/rooms/nonexistent/locks/trajectory:meta")
    assert response.status_code == 200
    lock_status = response.json()

    assert lock_status["locked"] is False
    assert lock_status["target"] == "trajectory:meta"
    assert "holder" not in lock_status
    assert "metadata" not in lock_status
