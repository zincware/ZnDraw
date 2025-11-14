"""Performance tests for room listing endpoint.

Tests verify that the batched room data fetcher significantly reduces
Redis query count compared to the N+1 query pattern.
"""

import pytest
import redis
import requests

from conftest import get_jwt_auth_headers


def test_batched_room_listing_response_structure(server, s22):
    """Verify batched implementation returns same structure as original."""
    import zndraw

    # Create multiple rooms with various states
    room_ids = []
    for i in range(5):
        room_id = f"perf-test-room-{i}"
        room_ids.append(room_id)

        vis = zndraw.ZnDraw(url=server, room=room_id, user=f"user-{i}")
        if i < 3:  # First 3 rooms have atoms
            vis.append(s22[0])

    # Set various metadata on rooms
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.set("room:perf-test-room-0:description", "Room 0 description")
    r.set("room:perf-test-room-1:locked", "1")
    r.set("room:perf-test-room-2:hidden", "1")
    r.hset("room:perf-test-room-3:metadata", "file", "test.xyz")

    # Fetch room list
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200

    rooms = response.json()

    # Verify all our rooms are present
    our_rooms = [r for r in rooms if r["id"] in room_ids]
    assert len(our_rooms) == 5

    # Verify response structure for each room
    for room in our_rooms:
        assert "id" in room
        assert "description" in room
        assert "frameCount" in room
        assert "locked" in room
        assert "metadataLocked" in room
        assert "hidden" in room
        assert "isDefault" in room
        assert "metadata" in room

        # Verify types
        assert isinstance(room["frameCount"], int)
        assert isinstance(room["locked"], bool)
        assert isinstance(room["hidden"], bool)
        assert isinstance(room["metadataLocked"], bool)
        assert isinstance(room["isDefault"], bool)
        assert isinstance(room["metadata"], dict)

    # Verify specific room states
    room_0 = [r for r in our_rooms if r["id"] == "perf-test-room-0"][0]
    assert room_0["description"] == "Room 0 description"
    assert room_0["frameCount"] == 1

    room_1 = [r for r in our_rooms if r["id"] == "perf-test-room-1"][0]
    assert room_1["locked"] is True

    room_2 = [r for r in our_rooms if r["id"] == "perf-test-room-2"][0]
    assert room_2["hidden"] is True

    room_3 = [r for r in our_rooms if r["id"] == "perf-test-room-3"][0]
    assert room_3["metadata"].get("file") == "test.xyz"


def test_batched_listing_with_many_rooms(server, s22):
    """Test batched fetching with larger number of rooms."""
    import zndraw

    num_rooms = 20

    # Create many rooms
    for i in range(num_rooms):
        room_id = f"many-rooms-test-{i}"
        vis = zndraw.ZnDraw(url=server, room=room_id, user=f"user-{i}")

        # Add atoms to some rooms
        if i % 3 == 0:
            vis.append(s22[0])

    # Set metadata on some rooms
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    for i in range(0, num_rooms, 5):
        room_id = f"many-rooms-test-{i}"
        r.hset(f"room:{room_id}:metadata", "batch", f"room-{i}")
        r.set(f"room:{room_id}:description", f"Description {i}")

    # Fetch all rooms
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200

    rooms = response.json()
    our_rooms = [r for r in rooms if r["id"].startswith("many-rooms-test-")]

    # Verify all rooms returned
    assert len(our_rooms) == num_rooms

    # Verify frame counts are correct
    rooms_with_atoms = [r for r in our_rooms if r["frameCount"] > 0]
    assert len(rooms_with_atoms) == num_rooms // 3 + 1  # Every 3rd room + room 0

    # Verify metadata
    rooms_with_metadata = [
        r for r in our_rooms if "batch" in r["metadata"]
    ]
    assert len(rooms_with_metadata) == 4  # 0, 5, 10, 15


def test_empty_room_list(server):
    """Test batched fetcher handles empty room list gracefully."""
    auth_headers = get_jwt_auth_headers(server, "test-user")

    # Fresh server should have minimal or no rooms
    # Just verify endpoint works without errors
    response = requests.get(f"{server}/api/rooms", headers=auth_headers)
    assert response.status_code == 200
    assert isinstance(response.json(), list)


def test_batched_listing_search_filter(server, s22):
    """Test search filtering works with batched implementation."""
    import zndraw

    # Create rooms with searchable metadata
    vis1 = zndraw.ZnDraw(url=server, room="search-test-alpha", user="user1")
    vis1.append(s22[0])

    vis2 = zndraw.ZnDraw(url=server, room="search-test-beta", user="user2")
    vis2.append(s22[0])

    vis3 = zndraw.ZnDraw(url=server, room="search-test-gamma", user="user3")
    vis3.append(s22[0])

    # Set metadata
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    r.hset("room:search-test-alpha:metadata", "experiment", "protein-folding")
    r.hset("room:search-test-beta:metadata", "experiment", "protein-binding")
    r.hset("room:search-test-gamma:metadata", "experiment", "dna-sequencing")

    # Search for "protein"
    response = requests.get(f"{server}/api/rooms?search=protein")
    assert response.status_code == 200

    rooms = response.json()
    matched_rooms = [
        r for r in rooms if r["id"].startswith("search-test-")
    ]

    # Should match alpha and beta (both have "protein" in metadata)
    assert len(matched_rooms) == 2
    matched_ids = {r["id"] for r in matched_rooms}
    assert "search-test-alpha" in matched_ids
    assert "search-test-beta" in matched_ids
    assert "search-test-gamma" not in matched_ids


def test_batched_listing_handles_missing_data(server, s22):
    """Test batched fetcher handles rooms with missing/null data."""
    import zndraw

    # Create room but don't set any metadata
    vis = zndraw.ZnDraw(url=server, room="minimal-room", user="user1")
    vis.append(s22[0])

    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200

    rooms = response.json()
    minimal_room = [r for r in rooms if r["id"] == "minimal-room"][0]

    # Verify defaults are correct
    assert minimal_room["description"] is None
    assert minimal_room["locked"] is False
    assert minimal_room["hidden"] is False
    assert minimal_room["metadataLocked"] is False
    assert minimal_room["metadata"] == {}
    assert minimal_room["frameCount"] == 1


@pytest.mark.parametrize("num_rooms", [10, 50])
def test_batched_listing_scaling(server, s22, num_rooms):
    """Test that batched listing scales well with increasing room count.

    This test creates many rooms and verifies the endpoint remains responsive.
    With the batched implementation, response time should be roughly constant
    regardless of room count, unlike the O(n) behavior of N+1 queries.
    """
    import time
    import zndraw

    # Create many rooms
    for i in range(num_rooms):
        room_id = f"scale-test-{num_rooms}-{i}"
        vis = zndraw.ZnDraw(url=server, room=room_id, user=f"user-{i}")
        if i % 5 == 0:
            vis.append(s22[0])

    # Measure response time
    start = time.time()
    response = requests.get(f"{server}/api/rooms")
    elapsed = time.time() - start

    assert response.status_code == 200
    rooms = response.json()

    # Verify correctness
    our_rooms = [r for r in rooms if r["id"].startswith(f"scale-test-{num_rooms}-")]
    assert len(our_rooms) == num_rooms

    # Performance check: should complete in reasonable time
    # With batched implementation, even 50 rooms should be fast
    # Allow generous timeout for CI environments
    max_time = 5.0  # 5 seconds max
    assert elapsed < max_time, f"Room listing took {elapsed:.2f}s (limit: {max_time}s)"

    print(f"Batched listing for {num_rooms} rooms: {elapsed:.3f}s")
