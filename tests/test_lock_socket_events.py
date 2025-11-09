"""Tests for lock socket events (lock:update)."""
import json
import time

import pytest
import redis
import requests
import socketio

from conftest import get_jwt_auth_headers


@pytest.fixture
def socketio_client():
    """Create a Socket.IO client for testing."""
    client = socketio.Client()
    yield client
    if client.connected:
        client.disconnect()


@pytest.fixture
def authenticated_session(server):
    """Join a room and return server, room, session_id, auth_headers, socketio_client."""
    room = "test-lock-socket"
    user_name = "test-user"

    # Get auth headers
    auth_headers = get_jwt_auth_headers(server, user_name)

    # Join the room to get session ID
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers
    )
    assert response.status_code == 200
    join_data = response.json()
    session_id = join_data["sessionId"]

    # Create socketio client
    client = socketio.Client()

    return server, room, session_id, auth_headers, client, user_name


def test_lock_acquire_emits_lock_update(authenticated_session):
    """Test that acquiring a lock emits lock:update event with action=acquired."""
    server, room, session_id, auth_headers, sio_client, user_name = authenticated_session

    # Track received events
    received_events = []

    def on_lock_update(data):
        received_events.append(data)

    sio_client.on("lock:update", on_lock_update)

    # Connect to socket
    sio_client.connect(server, auth={"token": auth_headers["Authorization"].split(" ")[1]})
    sio_client.emit("join:room", {"roomId": room})

    # Give socket time to connect and join room
    time.sleep(0.5)

    # Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "Creating test geometry"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Wait for socket event
    time.sleep(0.5)

    # Verify lock:update event was received
    assert len(received_events) == 1
    event = received_events[0]
    assert event["roomId"] == room
    assert event["target"] == "trajectory:meta"
    assert event["action"] == "acquired"
    assert event["holder"] == user_name
    assert event["message"] == "Creating test geometry"
    assert event["timestamp"] is not None
    assert event["sessionId"] == session_id  # Session ID should be included

    sio_client.disconnect()


def test_lock_refresh_emits_lock_update(authenticated_session):
    """Test that refreshing a lock emits lock:update event with action=refreshed."""
    server, room, session_id, auth_headers, sio_client, user_name = authenticated_session

    # Acquire lock first
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "Initial message"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # Track received events (after initial acquire)
    received_events = []

    def on_lock_update(data):
        received_events.append(data)

    sio_client.on("lock:update", on_lock_update)

    # Connect to socket
    sio_client.connect(server, auth={"token": auth_headers["Authorization"].split(" ")[1]})
    sio_client.emit("join:room", {"roomId": room})

    # Give socket time to connect
    time.sleep(0.5)

    # Refresh lock with new message
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/refresh",
        json={"lockToken": lock_token, "msg": "Updated message"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Wait for socket event
    time.sleep(0.5)

    # Verify lock:update event was received
    assert len(received_events) == 1
    event = received_events[0]
    assert event["roomId"] == room
    assert event["target"] == "trajectory:meta"
    assert event["action"] == "refreshed"
    assert event["holder"] == user_name
    assert event["message"] == "Updated message"
    assert event["timestamp"] is not None
    assert event["sessionId"] == session_id

    sio_client.disconnect()


def test_lock_release_emits_lock_update(authenticated_session):
    """Test that releasing a lock emits lock:update event with action=released."""
    server, room, session_id, auth_headers, sio_client, user_name = authenticated_session

    # Acquire lock first
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "Test lock"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # Track received events (after initial acquire)
    received_events = []

    def on_lock_update(data):
        received_events.append(data)

    sio_client.on("lock:update", on_lock_update)

    # Connect to socket
    sio_client.connect(server, auth={"token": auth_headers["Authorization"].split(" ")[1]})
    sio_client.emit("join:room", {"roomId": room})

    # Give socket time to connect
    time.sleep(0.5)

    # Release lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/release",
        json={"lockToken": lock_token},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Wait for socket event
    time.sleep(0.5)

    # Verify lock:update event was received
    assert len(received_events) == 1
    event = received_events[0]
    assert event["roomId"] == room
    assert event["target"] == "trajectory:meta"
    assert event["action"] == "released"
    assert event["holder"] is None
    assert event["message"] is None
    assert event["timestamp"] is None
    assert event["sessionId"] == session_id

    sio_client.disconnect()


def test_lock_acquire_without_message(authenticated_session):
    """Test that acquiring a lock without a message still emits lock:update."""
    server, room, session_id, auth_headers, sio_client, user_name = authenticated_session

    # Track received events
    received_events = []

    def on_lock_update(data):
        received_events.append(data)

    sio_client.on("lock:update", on_lock_update)

    # Connect to socket
    sio_client.connect(server, auth={"token": auth_headers["Authorization"].split(" ")[1]})
    sio_client.emit("join:room", {"roomId": room})

    # Give socket time to connect
    time.sleep(0.5)

    # Acquire lock without message
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    assert response.json()["success"] is True

    # Wait for socket event
    time.sleep(0.5)

    # Verify lock:update event was received
    assert len(received_events) == 1
    event = received_events[0]
    assert event["roomId"] == room
    assert event["target"] == "trajectory:meta"
    assert event["action"] == "acquired"
    assert event["holder"] == user_name
    assert event["message"] is None
    assert event["timestamp"] is None

    sio_client.disconnect()


def test_lock_events_for_different_targets(authenticated_session):
    """Test that lock events work for different lock targets."""
    server, room, session_id, auth_headers, sio_client, user_name = authenticated_session

    # Track received events
    received_events = []

    def on_lock_update(data):
        received_events.append(data)

    sio_client.on("lock:update", on_lock_update)

    # Connect to socket
    sio_client.connect(server, auth={"token": auth_headers["Authorization"].split(" ")[1]})
    sio_client.emit("join:room", {"roomId": room})

    # Give socket time to connect
    time.sleep(0.5)

    # Acquire lock for custom target
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/geometry:selection/acquire",
        json={"msg": "Selecting geometry"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200

    # Wait for socket event
    time.sleep(0.5)

    # Verify lock:update event was received for custom target
    assert len(received_events) == 1
    event = received_events[0]
    assert event["target"] == "geometry:selection"
    assert event["action"] == "acquired"
    assert event["holder"] == user_name

    sio_client.disconnect()


def test_multiple_clients_receive_lock_events(server):
    """Test that multiple clients in the same room receive lock:update events."""
    room = "test-multi-client"
    user1_name = "user1"
    user2_name = "user2"

    # User 1 setup
    user1_headers = get_jwt_auth_headers(server, user1_name)
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=user1_headers
    )
    user1_session = response.json()["sessionId"]

    # User 2 setup (listener)
    user2_headers = get_jwt_auth_headers(server, user2_name)
    requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=user2_headers)

    # Create socketio client for user 2 (listener)
    user2_events = []

    def on_lock_update(data):
        user2_events.append(data)

    sio_client = socketio.Client()
    sio_client.on("lock:update", on_lock_update)
    sio_client.connect(server, auth={"token": user2_headers["Authorization"].split(" ")[1]})
    sio_client.emit("join:room", {"roomId": room})

    # Give socket time to connect
    time.sleep(0.5)

    # User 1 acquires lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/trajectory:meta/acquire",
        json={"msg": "User 1 lock"},
        headers={**user1_headers, "X-Session-ID": user1_session},
    )
    assert response.status_code == 200

    # Wait for socket event
    time.sleep(0.5)

    # Verify user 2 received the lock event from user 1
    assert len(user2_events) == 1
    event = user2_events[0]
    assert event["holder"] == user1_name
    assert event["message"] == "User 1 lock"
    assert event["sessionId"] == user1_session  # Should include user1's sessionId

    sio_client.disconnect()


def test_lock_event_includes_session_id_for_filtering(authenticated_session):
    """Test that lock:update events include sessionId so clients can filter their own events."""
    server, room, session_id, auth_headers, sio_client, user_name = authenticated_session

    # Track all received events
    received_events = []

    def on_lock_update(data):
        received_events.append(data)

    sio_client.on("lock:update", on_lock_update)

    # Connect to socket
    sio_client.connect(server, auth={"token": auth_headers["Authorization"].split(" ")[1]})
    sio_client.emit("join:room", {"roomId": room})
    time.sleep(0.5)

    # Acquire lock
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/acquire",
        json={"msg": "Testing step lock"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    lock_token = response.json()["lockToken"]

    # Wait for socket event
    time.sleep(0.5)

    # Verify event was received with sessionId
    assert len(received_events) == 1
    acquire_event = received_events[0]
    assert acquire_event["target"] == "step"
    assert acquire_event["action"] == "acquired"
    assert acquire_event["sessionId"] == session_id
    assert acquire_event["holder"] == user_name

    # Client can now filter: if event.sessionId == mySessionId, ignore it
    my_session_id = session_id
    should_ignore = acquire_event["sessionId"] == my_session_id
    assert should_ignore is True  # This is our own event, should be ignored

    # Refresh the lock
    received_events.clear()
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/refresh",
        json={"lockToken": lock_token, "msg": "Still updating"},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    time.sleep(0.5)

    # Verify refresh event also includes sessionId
    assert len(received_events) == 1
    refresh_event = received_events[0]
    assert refresh_event["action"] == "refreshed"
    assert refresh_event["sessionId"] == session_id

    # Release the lock
    received_events.clear()
    response = requests.post(
        f"{server}/api/rooms/{room}/locks/step/release",
        json={"lockToken": lock_token},
        headers={**auth_headers, "X-Session-ID": session_id},
    )
    assert response.status_code == 200
    time.sleep(0.5)

    # Verify release event includes sessionId
    assert len(received_events) == 1
    release_event = received_events[0]
    assert release_event["action"] == "released"
    assert release_event["sessionId"] == session_id

    sio_client.disconnect()
