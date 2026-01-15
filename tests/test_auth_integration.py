"""Integration tests for JWT authentication flow."""

import pytest
import requests
from socketio.exceptions import ConnectionError as SocketIOConnectionError

from zndraw import ZnDraw


def test_login_flow(server):
    """Test full auth flow - register user, then login to get JWT token."""
    # Step 1: Register user
    register_response = requests.post(
        f"{server}/api/user/register", json={"userName": "test-user"}
    )
    assert register_response.status_code == 201
    reg_data = register_response.json()
    assert reg_data["status"] == "ok"
    assert reg_data["userName"] == "test-user"

    # Step 2: Login to get JWT
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


def test_register_without_username_creates_anonymous_guest(server):
    """Test that register without username creates anonymous guest user."""
    # Step 1: Register without username (server generates one)
    register_response = requests.post(f"{server}/api/user/register", json={})
    assert register_response.status_code == 201
    reg_data = register_response.json()
    assert reg_data["status"] == "ok"
    # Anonymous guest gets generated username like "user-abc123"
    assert reg_data["userName"].startswith("user-")

    # Step 2: Login with generated username
    response = requests.post(
        f"{server}/api/login", json={"userName": reg_data["userName"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert "token" in data
    # In local mode (default), all users are admin
    assert data["role"] == "admin"


def test_register_with_empty_username_creates_anonymous_guest(server):
    """Test that register with empty/whitespace username creates anonymous guest."""
    # Step 1: Register with empty username (server generates one)
    register_response = requests.post(
        f"{server}/api/user/register", json={"userName": "   "}
    )
    assert register_response.status_code == 201
    reg_data = register_response.json()
    assert reg_data["status"] == "ok"
    # Empty username is treated same as no username - generates guest user
    assert reg_data["userName"].startswith("user-")

    # Step 2: Login with generated username
    response = requests.post(
        f"{server}/api/login", json={"userName": reg_data["userName"]}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    # In local mode (default), all users are admin
    assert data["role"] == "admin"


def test_join_room_with_jwt_succeeds(server, get_jwt_auth_headers):
    """Test joining a room via socket with valid JWT token."""
    import socketio

    # Login to get JWT
    headers = get_jwt_auth_headers(server, "test-user")
    jwt_token = headers["Authorization"].replace("Bearer ", "")

    # Create room first
    room = "test-room"
    create_response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
    )
    assert create_response.status_code == 201

    # Connect socket and join room
    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    try:
        response = sio.call("room:join", {"roomId": room, "clientType": "frontend"})
        assert response["status"] == "ok"
        assert "sessionId" in response
        assert "frameCount" in response
        assert "step" in response
        assert "locked" in response
    finally:
        sio.disconnect()


def test_join_room_without_jwt_fails(server):
    """Test that socket connection without JWT token fails."""
    import socketio
    from socketio.exceptions import ConnectionError as SocketIOConnectionError

    sio = socketio.Client()

    # Try to connect without auth token
    with pytest.raises(SocketIOConnectionError):
        sio.connect(server, wait=True)


def test_join_room_with_invalid_jwt_fails(server):
    """Test that socket connection with invalid JWT token fails."""
    import socketio
    from socketio.exceptions import ConnectionError as SocketIOConnectionError

    sio = socketio.Client()

    # Try to connect with invalid token
    with pytest.raises(SocketIOConnectionError):
        sio.connect(server, auth={"token": "invalid-token-here"}, wait=True)


def test_join_room_with_malformed_auth_header_fails(server):
    """Test that socket connection without proper token fails."""
    import socketio
    from socketio.exceptions import ConnectionError as SocketIOConnectionError

    sio = socketio.Client()

    # Try to connect with empty auth (no token key)
    with pytest.raises(SocketIOConnectionError):
        sio.connect(server, auth={}, wait=True)


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

    # Register user first
    register_response = requests.post(
        f"{server}/api/user/register", json={"userName": username}
    )
    assert register_response.status_code == 201

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

    username = "test-user-claims"

    # Register user first
    register_response = requests.post(
        f"{server}/api/user/register", json={"userName": username}
    )
    assert register_response.status_code == 201

    # Login to get JWT
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


def test_user_created_on_register(server, redis_client):
    """Test that user is created in Redis when /api/user/register is called.

    In the new flow, users are created via /api/user/register (the ONLY place
    for user creation). Login and socket.connect only authenticate/validate.
    """
    from zndraw.app.redis_keys import UserKeys

    username = "redis-register-test-user"

    # Check user doesn't exist yet
    keys = UserKeys(username)
    assert redis_client.exists(keys.hash_key()) == 0

    # Register - should create Redis entry
    response = requests.post(f"{server}/api/user/register", json={"userName": username})
    assert response.status_code == 201

    # Check Redis - user should NOW exist
    assert redis_client.exists(keys.hash_key()) == 1

    # Verify stored data
    user_data = redis_client.hgetall(keys.hash_key())
    assert user_data["userName"] == username
    assert "createdAt" in user_data


def test_login_without_registration_fails(server):
    """Test that login fails if user doesn't exist (wasn't registered).

    Note: Error message is generic to prevent username enumeration attacks.
    """
    username = "unregistered-user"

    # Login without registering first - should fail
    response = requests.post(f"{server}/api/login", json={"userName": username})

    assert response.status_code == 401
    data = response.json()
    assert "error" in data
    # Generic error to prevent username enumeration
    assert data["error"] == "Authentication failed"


def test_join_room_updates_redis_room_users(server, redis_client, get_jwt_auth_headers):
    """Test that joining a room adds user to room's user set in Redis."""
    import socketio

    # Login and get headers
    headers = get_jwt_auth_headers(server, "room-joiner")
    jwt_token = headers["Authorization"].replace("Bearer ", "")

    # Extract username from token (decode without verification)
    import jwt as pyjwt

    payload = pyjwt.decode(jwt_token, options={"verify_signature": False})
    user_name = payload["sub"]

    # Create room first
    room = "redis-room-test"
    create_response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
    )
    assert create_response.status_code == 201

    # Join room via socket
    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    try:
        response = sio.call("room:join", {"roomId": room, "clientType": "frontend"})
        assert response["status"] == "ok"

        # Verify user is in room's user set
        room_users_key = f"room:{room}:users"
        assert redis_client.sismember(room_users_key, user_name) == 1
    finally:
        sio.disconnect()


def test_join_nonexistent_room_fails(server, get_jwt_auth_headers):
    """Test that joining a non-existent room via socket returns 404.

    Room creation and room joining are separate operations.
    Use POST /api/rooms to create a room first.
    """
    import socketio

    headers = get_jwt_auth_headers(server, "test-user")
    jwt_token = headers["Authorization"].replace("Bearer ", "")

    # Connect socket
    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    try:
        room = "nonexistent-room-12345"
        response = sio.call("room:join", {"roomId": room, "clientType": "frontend"})

        assert response["status"] == "error"
        assert response["code"] == 404
        assert "not found" in response["message"].lower()
    finally:
        sio.disconnect()


def test_explicit_room_creation(server, get_jwt_auth_headers):
    """Test creating a room explicitly via POST /api/rooms."""
    headers = get_jwt_auth_headers(server, "room-creator")

    room = "explicit-room-test"
    response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room, "description": "Test room"},
        headers=headers,
    )

    assert response.status_code == 201
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room
    assert data["created"] is True


def test_explicit_room_creation_duplicate_fails(server, get_jwt_auth_headers):
    """Test that creating a room that already exists fails."""
    headers = get_jwt_auth_headers(server, "room-creator")

    room = "duplicate-room-test"

    # Create room first
    response1 = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
    )
    assert response1.status_code == 201

    # Try to create again - should fail
    response2 = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
    )
    assert response2.status_code == 409
    data = response2.json()
    assert "error" in data


def test_join_after_explicit_creation(server, get_jwt_auth_headers):
    """Test joining a room via socket after explicit creation."""
    import socketio

    headers = get_jwt_auth_headers(server, "test-user")
    jwt_token = headers["Authorization"].replace("Bearer ", "")

    room = "join-after-create-test"

    # Create room explicitly
    response1 = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=headers,
    )
    assert response1.status_code == 201

    # Join via socket - should succeed since room exists
    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    try:
        response = sio.call("room:join", {"roomId": room, "clientType": "frontend"})
        assert response["status"] == "ok"
        assert "sessionId" in response
    finally:
        sio.disconnect()
