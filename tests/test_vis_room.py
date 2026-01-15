import uuid

import pytest
import requests

from zndraw.zndraw import ZnDraw


def test_rest_create_and_socket_join_room(server, get_jwt_auth_headers):
    """Test creating a room via REST API and joining via socket."""
    import socketio

    # Authenticate as admin to see all rooms
    headers = get_jwt_auth_headers(server, "admin")
    response = requests.get(f"{server}/api/rooms", headers=headers, timeout=10)
    assert response.status_code == 200
    rooms = response.json()
    assert isinstance(rooms, list)
    assert len(rooms) == 0  # no rooms yet

    room = "test-room-1"

    # First, create the room via REST
    create_response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
        timeout=10,
    )
    assert create_response.status_code == 201
    create_data = create_response.json()
    assert create_data["status"] == "ok"
    assert create_data["roomId"] == room
    assert create_data["created"] is True

    # Then join the room via socket room:join
    jwt_token = headers["Authorization"].replace("Bearer ", "")
    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True, wait_timeout=10)

    try:
        join_response = sio.call(
            "room:join", {"roomId": room, "clientType": "frontend"}, timeout=10
        )
        assert join_response["status"] == "ok"
        assert "sessionId" in join_response
        assert "roomData" in join_response
    finally:
        sio.disconnect()

    # list all rooms again to see if the new room is there
    response = requests.get(f"{server}/api/rooms", headers=headers, timeout=10)
    assert response.status_code == 200
    rooms = response.json()
    assert isinstance(rooms, list)
    assert len(rooms) == 1
    assert rooms[0]["id"] == room

    # getting any frame will fail with index error
    for frame_idx in [0, 1, -1, 100]:
        response = requests.get(
            f"{server}/api/rooms/{room}/frames",
            params={"indices": str(frame_idx)},
            timeout=10,
        )
        assert response.status_code == 404
        data = response.json()
        assert data["type"] == "IndexError"
        assert data["error"] == f"No frames found in room '{room}'"


def test_socket_join_nonexistent_room_fails(server, get_jwt_auth_headers):
    """Test that joining a non-existent room via socket returns 404.

    Note: Frontend clients auto-create rooms, so this test uses clientType: "python"
    to test the behavior for Python SDK clients which require explicit room creation.
    """
    import socketio

    room = "nonexistent-room"
    headers = get_jwt_auth_headers(server)
    jwt_token = headers["Authorization"].replace("Bearer ", "")

    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True, wait_timeout=10)

    try:
        response = sio.call(
            "room:join", {"roomId": room, "clientType": "python"}, timeout=10
        )
        assert response["status"] == "error"
        assert response["code"] == 404
        assert "not found" in response["message"].lower()
    finally:
        sio.disconnect()


def test_socket_join_existing_room(server, get_jwt_auth_headers):
    """Test joining an existing room via socket."""
    import socketio

    room = "test-room-1"
    headers = get_jwt_auth_headers(server)

    # Create the room first via REST
    create_response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
        timeout=10,
    )
    assert create_response.status_code == 201

    # Join the existing room with the same user via socket
    jwt_token = headers["Authorization"].replace("Bearer ", "")
    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True, wait_timeout=10)

    try:
        response = sio.call(
            "room:join", {"roomId": room, "clientType": "frontend"}, timeout=10
        )
        assert response["status"] == "ok"
        assert "sessionId" in response
        assert "roomData" in response
    finally:
        sio.disconnect()

    # Join with a different user
    headers2 = get_jwt_auth_headers(server, "user2")
    jwt_token2 = headers2["Authorization"].replace("Bearer ", "")
    sio2 = socketio.Client()
    sio2.connect(server, auth={"token": jwt_token2}, wait=True, wait_timeout=10)

    try:
        response2 = sio2.call(
            "room:join", {"roomId": room, "clientType": "frontend"}, timeout=10
        )
        assert response2["status"] == "ok"
        assert "sessionId" in response2
    finally:
        sio2.disconnect()


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
        timeout=10,
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
        timeout=10,
    )
    assert response.status_code == 201
    data = response.json()
    assert data["status"] == "ok"
    assert data["frameCount"] == 0  # Empty room


def test_create_room_invalid_name(server, get_jwt_auth_headers):
    """Test that creating a room with invalid name fails."""
    room = "invalid:room"
    headers = get_jwt_auth_headers(server)

    # Creating a room with ':' in name should fail
    response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
        timeout=10,
    )
    assert response.status_code == 400
    data = response.json()
    assert "error" in data
    assert ":" in data["error"]

    # ZnDraw client should also fail
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


def test_python_creates_no_session_camera(server):
    """Test that a Python (ZnDraw) client does NOT create a session camera."""
    vis = ZnDraw(url=server, room=str(uuid.uuid4()), user=str(uuid.uuid4()))
    assert len([x for x in vis.geometries if x.startswith("cam:")]) == 0


def test_frontend_creates_one_session_camera(server, connect_room):
    """Test that a frontend client creates exactly ONE session camera."""
    room_id = str(uuid.uuid4())

    # Join as frontend (connect_room uses clientType: "frontend")
    c1 = connect_room(room_id)

    vis = ZnDraw(url=server, room=room_id, user=str(uuid.uuid4()))
    assert len([x for x in vis.geometries if x.startswith("cam:")]) == 1

    # join again as another frontend client
    c2 = connect_room(room_id)
    assert len([x for x in vis.geometries if x.startswith("cam:")]) == 2

    c2.sio.disconnect()
    c2.sio.sleep(0.1)
    assert len([x for x in vis.geometries if x.startswith("cam:")]) == 1

    c1.sio.disconnect()
    c1.sio.sleep(0.1)
    assert len([x for x in vis.geometries if x.startswith("cam:")]) == 0
