"""Tests for room lock enforcement on modification operations."""

import json
import uuid

import pytest


@pytest.fixture
def test_app(redis_client):
    """Create Flask app with external Redis for testing."""
    from zndraw.config import ZnDrawConfig
    from zndraw.server import create_app

    # Use external Redis for testing so we can manipulate state
    config = ZnDrawConfig(redis_url="redis://localhost:6379")
    app = create_app(config=config)
    app.config["TESTING"] = True

    yield app


@pytest.fixture
def client(test_app):
    """Create Flask test client."""
    return test_app.test_client()


@pytest.fixture
def test_user(client):
    """Create a test user and return username and token."""
    user_name = f"test-user-{uuid.uuid4().hex[:8]}"
    # Register user first
    register_response = client.post("/api/user/register", json={"userName": user_name})
    assert register_response.status_code == 201
    # Then login
    response = client.post("/api/login", json={"userName": user_name})
    assert response.status_code == 200
    token = response.get_json()["token"]
    return {"userName": user_name, "token": token}


@pytest.fixture
def auth_headers(test_user):
    """Get authentication headers from test user."""
    return {"Authorization": f"Bearer {test_user['token']}"}


@pytest.fixture
def session_headers(redis_client, test_user):
    """Create session headers with valid session ID."""
    from zndraw.app.redis_keys import SessionKeys

    session_id = f"test-session-{uuid.uuid4().hex}"
    session_key = SessionKeys.session_data(session_id)

    # Create valid session data with the same user
    session_data = {
        "userId": test_user["userName"],
        "roomId": "test-room",
        "createdAt": "2025-11-27T00:00:00Z",
    }
    redis_client.set(session_key, json.dumps(session_data), ex=3600)

    return {"X-Session-ID": session_id}


@pytest.fixture
def acquire_lock(redis_client, test_user):
    """Fixture to acquire a lock for testing."""

    def _acquire(room_id: str, target: str, session_headers_dict: dict) -> str:
        from zndraw.app.redis_keys import RoomKeys

        session_id = session_headers_dict["X-Session-ID"]
        lock_token = f"token-{uuid.uuid4().hex}"

        keys = RoomKeys(room_id)
        lock_key = keys.lock(target)

        # Set the lock
        lock_data = {
            "sessionId": session_id,
            "userId": test_user["userName"],
            "token": lock_token,
        }
        redis_client.set(lock_key, json.dumps(lock_data), ex=60)

        return lock_token

    return _acquire


def test_locked_room_prevents_frame_deletion(
    client, redis_client, auth_headers, session_headers
):
    """Test that deleting frames is blocked when room is locked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)

    # Lock the room
    redis_client.set(keys.locked(), "1")

    # Try to delete frames - should fail with 403
    response = client.delete(
        f"/api/rooms/{room_id}/frames?frame_id=0",
        headers={**auth_headers, **session_headers},
    )

    assert response.status_code == 403
    assert "locked" in response.json["error"].lower()


def test_unlocked_room_allows_frame_deletion(
    client, redis_client, auth_headers, session_headers, acquire_lock
):
    """Test that deleting frames works when room is unlocked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)

    # Ensure room is unlocked
    redis_client.set(keys.locked(), "0")

    # Acquire trajectory lock
    acquire_lock(room_id, "trajectory:meta", session_headers)

    # Delete should succeed (will fail with 404 if no frames, but not 403)
    response = client.delete(
        f"/api/rooms/{room_id}/frames?frame_id=0",
        headers={**auth_headers, **session_headers},
    )

    assert response.status_code != 403


@pytest.mark.parametrize(
    "endpoint,method,data",
    [
        ("/api/rooms/{room_id}/frames?frame_id=0", "DELETE", None),
        ("/api/rooms/{room_id}/frames?action=append", "POST", b"msgpack_data"),
        ("/api/rooms/{room_id}/frames/bulk?start=0&stop=1", "PATCH", b"msgpack_data"),
    ],
)
def test_locked_room_blocks_frame_operations(
    client, redis_client, auth_headers, session_headers, endpoint, method, data
):
    """Test that all frame modification operations are blocked when room is locked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)
    redis_client.set(keys.locked(), "1")

    url = endpoint.format(room_id=room_id)

    if method == "DELETE":
        response = client.delete(url, headers={**auth_headers, **session_headers})
    elif method == "POST":
        response = client.post(
            url, data=data, headers={**auth_headers, **session_headers}
        )
    elif method == "PATCH":
        response = client.patch(
            url, data=data, headers={**auth_headers, **session_headers}
        )

    assert response.status_code == 403
    assert "locked" in response.json["error"].lower()


@pytest.mark.parametrize(
    "endpoint,method,data",
    [
        (
            "/api/rooms/{room_id}/geometries",
            "POST",
            {"key": "test", "type": "Sphere", "data": {}},
        ),
        ("/api/rooms/{room_id}/geometries/test", "DELETE", None),
        ("/api/rooms/{room_id}/figures", "POST", {"key": "test", "figure": {}}),
        ("/api/rooms/{room_id}/figures/test", "DELETE", None),
    ],
)
def test_locked_room_blocks_geometry_operations(
    client, redis_client, auth_headers, session_headers, endpoint, method, data
):
    """Test that geometry and figure operations are blocked when room is locked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)
    redis_client.set(keys.locked(), "1")

    url = endpoint.format(room_id=room_id)
    headers = {**auth_headers, **session_headers}

    if method == "POST":
        response = client.post(url, json=data, headers=headers)
    elif method == "DELETE":
        response = client.delete(url, headers=headers)

    assert response.status_code == 403
    assert "locked" in response.json["error"].lower()


@pytest.mark.parametrize(
    "endpoint,method,data",
    [
        (
            "/api/rooms/{room_id}/geometries",
            "POST",
            {"key": "test", "type": "Sphere", "data": {}},
        ),
        ("/api/rooms/{room_id}/geometries/test", "DELETE", None),
        ("/api/rooms/{room_id}/figures", "POST", {"key": "test", "figure": {}}),
        ("/api/rooms/{room_id}/figures/test", "DELETE", None),
    ],
)
def test_unlocked_room_allows_geometry_operations(
    client,
    redis_client,
    auth_headers,
    session_headers,
    acquire_lock,
    endpoint,
    method,
    data,
):
    """Test that geometry and figure operations work when room is unlocked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)

    # Ensure room is unlocked
    redis_client.set(keys.locked(), "0")

    # Acquire trajectory lock
    acquire_lock(room_id, "trajectory:meta", session_headers)

    url = endpoint.format(room_id=room_id)
    headers = {**auth_headers, **session_headers}

    if method == "POST":
        response = client.post(url, json=data, headers=headers)
    elif method == "DELETE":
        response = client.delete(url, headers=headers)

    # Should NOT return 403 (room locked error)
    # Will return 404 if trying to delete non-existent items, 200 for success
    assert response.status_code != 403


@pytest.mark.parametrize(
    "endpoint,method,data",
    [
        ("/api/rooms/{room_id}/bookmarks/0", "PUT", {"label": "Test Bookmark"}),
        ("/api/rooms/{room_id}/bookmarks/0", "DELETE", None),
    ],
)
def test_locked_room_blocks_bookmark_operations(
    client, redis_client, auth_headers, session_headers, endpoint, method, data
):
    """Test that bookmark operations are blocked when room is locked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)
    redis_client.set(keys.locked(), "1")

    url = endpoint.format(room_id=room_id)
    headers = {**auth_headers, **session_headers}

    if method == "PUT":
        response = client.put(url, json=data, headers=headers)
    elif method == "DELETE":
        response = client.delete(url, headers=headers)

    assert response.status_code == 403
    assert "locked" in response.json["error"].lower()


@pytest.mark.parametrize(
    "endpoint,method,data",
    [
        ("/api/rooms/{room_id}/bookmarks/0", "PUT", {"label": "Test Bookmark"}),
        ("/api/rooms/{room_id}/bookmarks/0", "DELETE", None),
    ],
)
def test_unlocked_room_allows_bookmark_operations(
    client,
    redis_client,
    auth_headers,
    session_headers,
    acquire_lock,
    endpoint,
    method,
    data,
):
    """Test that bookmark operations work when room is unlocked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)

    # Ensure room is unlocked
    redis_client.set(keys.locked(), "0")

    # Acquire trajectory lock
    acquire_lock(room_id, "trajectory:meta", session_headers)

    url = endpoint.format(room_id=room_id)
    headers = {**auth_headers, **session_headers}

    if method == "PUT":
        response = client.put(url, json=data, headers=headers)
    elif method == "DELETE":
        response = client.delete(url, headers=headers)

    # Should NOT return 403 (room locked error)
    # Will return 404 if trying to delete non-existent bookmarks, 200/400 for PUT
    assert response.status_code != 403


def test_locked_room_blocks_renormalize(
    client, redis_client, auth_headers, session_headers
):
    """Test that renormalize operation is blocked when room is locked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)
    redis_client.set(keys.locked(), "1")

    response = client.post(
        f"/api/rooms/{room_id}/renormalize", headers={**auth_headers, **session_headers}
    )

    assert response.status_code == 403
    assert "locked" in response.json["error"].lower()


def test_unlocked_room_allows_renormalize(
    client, redis_client, auth_headers, session_headers, acquire_lock
):
    """Test that renormalize operation works when room is unlocked."""
    from zndraw.app.redis_keys import RoomKeys

    room_id = "test-room"
    keys = RoomKeys(room_id)

    # Ensure room is unlocked
    redis_client.set(keys.locked(), "0")

    # Acquire trajectory lock
    acquire_lock(room_id, "trajectory:meta", session_headers)

    response = client.post(
        f"/api/rooms/{room_id}/renormalize", headers={**auth_headers, **session_headers}
    )

    # Should NOT return 403 (room locked error)
    # Will return 404 if room doesn't exist, 200 for success
    assert response.status_code != 403
