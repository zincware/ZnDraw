"""Integration tests for JWT authentication flow."""

import pytest
import requests
from conftest import get_jwt_auth_headers
from socketio.exceptions import ConnectionError as SocketIOConnectionError

from zndraw import ZnDraw


def test_login_flow(server):
    """Test full login flow - POST /api/login returns JWT token."""
    response = requests.post(f"{server}/api/login", json={"userName": "test-user"})

    assert response.status_code == 200
    data = response.json()

    # Verify response structure
    assert data["status"] == "ok"
    assert "token" in data
    assert "userName" in data
    assert "role" in data
    assert isinstance(data["token"], str)
    assert len(data["token"]) > 0
    assert data["userName"] == "test-user"


def test_login_without_username_creates_anonymous_guest(server):
    """Test that login without username creates anonymous guest user."""
    response = requests.post(f"{server}/api/login", json={})

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert "token" in data
    assert "userName" in data
    # Anonymous guest gets generated username like "user-abc123"
    assert data["userName"].startswith("user-")
    # In local mode (default), all users are admin
    assert data["role"] == "admin"


def test_login_with_empty_username_creates_anonymous_guest(server):
    """Test that login with empty/whitespace username creates anonymous guest."""
    response = requests.post(f"{server}/api/login", json={"userName": "   "})

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    # Empty username is treated same as no username - generates guest user
    assert data["userName"].startswith("user-")
    # In local mode (default), all users are admin
    assert data["role"] == "admin"


def test_join_room_with_jwt_succeeds(server):
    """Test joining a room with valid JWT token."""
    # Login to get JWT
    headers = get_jwt_auth_headers(server, "test-user")

    # Join room with JWT
    room = "test-room"
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=headers
    )

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room


def test_join_room_without_jwt_fails(server):
    """Test that joining a room without JWT token fails with 401."""
    room = "test-room"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})

    assert response.status_code == 401
    data = response.json()
    assert "error" in data


def test_join_room_with_invalid_jwt_fails(server):
    """Test that joining a room with invalid JWT token fails with 401."""
    room = "test-room"
    response = requests.post(
        f"{server}/api/rooms/{room}/join",
        json={},
        headers={"Authorization": "Bearer invalid-token-here"},
    )

    assert response.status_code == 401
    data = response.json()
    assert "error" in data


def test_join_room_with_malformed_auth_header_fails(server):
    """Test that malformed Authorization header fails."""
    room = "test-room"

    # Missing 'Bearer' prefix
    response = requests.post(
        f"{server}/api/rooms/{room}/join",
        json={},
        headers={"Authorization": "just-a-token"},
    )
    assert response.status_code == 401


def test_websocket_connection_with_jwt_succeeds(server):
    """Test WebSocket connection with valid JWT token."""
    # ZnDraw should automatically login and connect
    vis = ZnDraw(url=server, room="test-room", user="test-user")

    assert vis.socket.connected
    vis.socket.disconnect()


def test_websocket_connection_without_jwt_fails(server):
    """Test that WebSocket connection without JWT token fails."""
    import socketio

    sio = socketio.Client()

    # Try to connect without auth token
    with pytest.raises(SocketIOConnectionError):
        sio.connect(server, wait=True)


def test_websocket_connection_with_invalid_jwt_fails(server):
    """Test that WebSocket connection with invalid JWT fails."""
    import socketio

    sio = socketio.Client()

    # Try to connect with invalid token
    with pytest.raises(SocketIOConnectionError):
        sio.connect(server, auth={"token": "invalid-token"}, wait=True)


def test_multiple_logins_same_username_get_same_user(server):
    """Test that multiple logins with same username get the same user identity."""
    username = "duplicate-user"

    # First login
    response1 = requests.post(f"{server}/api/login", json={"userName": username})
    assert response1.status_code == 200
    data1 = response1.json()

    # Second login with same username
    response2 = requests.post(f"{server}/api/login", json={"userName": username})
    assert response2.status_code == 200
    data2 = response2.json()

    # Both should succeed with same username but different tokens
    assert data1["userName"] == data2["userName"] == username
    assert data1["token"] != data2["token"]  # Different tokens (new JTI)


def test_jwt_token_contains_correct_claims(server):
    """Test that JWT token contains expected claims."""
    import jwt as pyjwt

    username = "test-user"
    response = requests.post(f"{server}/api/login", json={"userName": username})
    assert response.status_code == 200
    data = response.json()
    token = data["token"]
    user_name = data["userName"]

    # Decode without verification to check claims
    payload = pyjwt.decode(token, options={"verify_signature": False})

    assert payload["sub"] == user_name
    assert "jti" in payload  # JWT ID should be present


def test_python_client_auto_login(server):
    """Test that ZnDraw Python client automatically logs in."""
    # ZnDraw should call login() automatically in __post_init__
    vis = ZnDraw(url=server, room="test-room", user="auto-login-user")

    # Verify client has JWT token
    assert vis.api.jwt_token is not None
    assert len(vis.api.jwt_token) > 0

    # Verify client can join room and perform operations
    assert vis.socket.connected
    assert len(vis) == 0  # Empty room

    vis.socket.disconnect()


def test_user_session_persists_in_redis(server, redis_client):
    """Test that user session data is stored in Redis after login."""
    username = "redis-test-user"

    # Login
    response = requests.post(f"{server}/api/login", json={"userName": username})
    assert response.status_code == 200
    data = response.json()
    user_name = data["userName"]

    # Check Redis for user session
    user_key = f"user:{user_name}"
    assert redis_client.exists(user_key) == 1

    # Verify stored data
    user_data = redis_client.hgetall(user_key)
    assert user_data["userName"] == username
    assert "createdAt" in user_data
    assert "lastLogin" in user_data


def test_join_room_updates_redis_room_users(server, redis_client):
    """Test that joining a room adds user to room's user set in Redis."""
    # Login and get headers
    headers = get_jwt_auth_headers(server, "room-joiner")

    # Extract username from token (decode without verification)
    import jwt as pyjwt

    token = headers["Authorization"].split(" ")[1]
    payload = pyjwt.decode(token, options={"verify_signature": False})
    user_name = payload["sub"]

    # Join room
    room = "redis-room-test"
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=headers
    )
    assert response.status_code == 200

    # Verify user is in room's user set
    room_users_key = f"room:{room}:users"
    assert redis_client.sismember(room_users_key, user_name) == 1
