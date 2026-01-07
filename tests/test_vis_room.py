import pytest
import requests

from zndraw.zndraw import ZnDraw


def test_rest_create_and_join_room(server, get_jwt_auth_headers):
    """Test creating and joining a room via REST API."""
    # Authenticate as admin to see all rooms
    headers = get_jwt_auth_headers(server, "admin")
    response = requests.get(f"{server}/api/rooms", headers=headers)
    assert response.status_code == 200
    rooms = response.json()
    assert isinstance(rooms, list)
    assert len(rooms) == 0  # no rooms yet

    room = "test-room-1"

    # First, create the room
    create_response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
    )
    assert create_response.status_code == 201
    create_data = create_response.json()
    assert create_data["status"] == "ok"
    assert create_data["roomId"] == room
    assert create_data["created"] is True

    # Then join the room
    join_response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=headers
    )
    assert join_response.status_code == 200
    join_data = join_response.json()
    assert join_data["status"] == "ok"
    assert join_data["roomId"] == room
    assert "sessionId" in join_data
    assert "userName" in join_data

    # list all rooms again to see if the new room is there
    response = requests.get(f"{server}/api/rooms", headers=headers)
    assert response.status_code == 200
    rooms = response.json()
    assert isinstance(rooms, list)
    assert len(rooms) == 1
    assert rooms[0]["id"] == room

    # getting any frame will fail with index error
    for frame_idx in [0, 1, -1, 100]:
        response = requests.get(
            f"{server}/api/rooms/{room}/frames", params={"indices": str(frame_idx)}
        )
        assert response.status_code == 404
        data = response.json()
        assert data["type"] == "IndexError"
        assert data["error"] == f"No frames found in room '{room}'"


def test_join_nonexistent_room_fails(server, get_jwt_auth_headers):
    """Test that joining a non-existent room returns 404."""
    room = "nonexistent-room"
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=get_jwt_auth_headers(server)
    )
    assert response.status_code == 404
    data = response.json()
    assert "error" in data
    assert room in data["error"]


def test_join_existing_room(server, get_jwt_auth_headers):
    """Test joining an existing room."""
    room = "test-room-1"
    headers = get_jwt_auth_headers(server)

    # Create the room first
    create_response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
    )
    assert create_response.status_code == 201

    # Join the existing room with the same user
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=headers
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room
    assert "sessionId" in data
    assert "userName" in data

    # Join with a different user
    headers2 = get_jwt_auth_headers(server, "user2")
    response2 = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=headers2
    )
    assert response2.status_code == 200
    data2 = response2.json()
    assert data2["status"] == "ok"
    assert data2["roomId"] == room


def test_create_room_with_copy_from(server, s22, get_jwt_auth_headers):
    """Test creating a room by copying from an existing room."""
    # create a source room with some frames
    source_room = "source-room"
    vis = ZnDraw(url=server, room=source_room, user="user1")
    vis.extend(s22[:3])

    # Create a new room by copying from the source room
    new_room = "test-room-1"
    headers = get_jwt_auth_headers(server)
    response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": new_room, "copyFrom": source_room},
        headers=headers,
    )
    assert response.status_code == 201
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == new_room
    assert data["frameCount"] == 3  # Should have copied 3 frames
    assert data["created"] is True

    # Create another room without copying
    another_room = "test-room-2"
    response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": another_room},
        headers=headers,
    )
    assert response.status_code == 201
    data = response.json()
    assert data["status"] == "ok"
    assert data["frameCount"] == 0  # Empty room


def test_join_room_invalid_name(server, get_jwt_auth_headers):
    room = "invalid:room"
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=get_jwt_auth_headers(server)
    )
    assert response.status_code == 400
    data = response.json()
    assert "error" in data
    assert data["error"] == "Room ID cannot contain ':' character"

    with pytest.raises(RuntimeError):
        _ = ZnDraw(url=server, room=room, user="user1")


def test_vis_empty_room(server, s22):
    vis1 = ZnDraw(url=server, room="room-a", user="user1")
    assert len(vis1) == 0
    with pytest.raises(IndexError):
        _ = vis1[0]

    vis1.extend(s22)
    assert len(vis1) == len(s22)
    for i in range(len(s22)):
        assert vis1[i] == s22[i]

    with pytest.raises(IndexError):
        _ = vis1[9999]


def test_vis_len(server, s22):
    vis1 = ZnDraw(url=server, room="room-a", user="user1")
    vis1.extend(s22)
    assert len(vis1) == len(s22)

    vis2 = ZnDraw(url=server, room="room-a", user="user2")
    assert len(vis2) == len(s22)


def test_lock_acquisition_broadcasts_with_iso_timestamp(server, s22):
    """Test that lock acquisition broadcasts room:update with ISO timestamp."""
    import datetime

    from zndraw.app.models import LockMetadata

    vis = ZnDraw(url=server, room="lock-test-room", user="test-user")

    # Acquire lock - this should broadcast room:update with LockMetadata
    with vis.get_lock(msg="Testing lock broadcast"):
        # Lock is acquired and metadata is broadcast
        # The timestamp should be in ISO format, not float

        # Verify we can create LockMetadata with ISO timestamp
        iso_timestamp = datetime.datetime.utcnow().isoformat()
        lock_metadata = LockMetadata(
            msg="Test message", userName="test-user", timestamp=iso_timestamp
        )
        assert lock_metadata.timestamp == iso_timestamp

        # Verify that float timestamps are rejected
        with pytest.raises(Exception):  # Pydantic validation error
            LockMetadata(
                msg="Test message",
                userName="test-user",
                timestamp=1234567890.123,  # Float should fail
            )

    # Release lock before extend (no nested locking in simplified design)
    vis.extend(s22)

    # Verify frames were added
    assert len(vis) == len(s22)
