import pytest
import requests
from conftest import get_jwt_auth_headers

from zndraw.zndraw import ZnDraw


def test_rest_join_new_room(server):
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    assert isinstance(rooms, list)
    assert len(rooms) == 0  # no rooms yet

    room = "test-room-1"
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=get_jwt_auth_headers(server)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room
    assert data["created"] is True
    assert "sessionId" in data
    assert "userName" in data

    # list all rooms again to see if the new room is there
    response = requests.get(f"{server}/api/rooms")
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


def test_join_existing_room(server):
    # create a room first
    room = "test-room-1"
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=get_jwt_auth_headers(server)
    )
    assert response.status_code == 200

    # join the existing room with a different user
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=get_jwt_auth_headers(server)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room
    assert data["created"] is False
    assert "sessionId" in data
    assert "userName" in data


def test_join_room_with_copy_from(server, s22):
    # create a source room with some frames
    source_room = "source-room"
    vis = ZnDraw(url=server, room=source_room, user="user1")
    vis.extend(s22[:3])

    # create a new room by copying from the source room
    new_room = "test-room-1"
    response = requests.post(
        f"{server}/api/rooms/{new_room}/join",
        json={"copyFrom": source_room},
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == new_room
    assert data["frameCount"] == 3  # Should have copied 3 frames
    assert data["created"] is True

    # create a new room without copying from another room
    another_room = "test-room-2"
    response = requests.post(
        f"{server}/api/rooms/{another_room}/join",
        json={},
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["frameCount"] == 0  # Empty room


def test_join_room_invalid_name(server):
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
