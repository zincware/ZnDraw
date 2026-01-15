import shutil
import signal
import socket
import subprocess
import time
import typing as t
import uuid
from dataclasses import dataclass

import ase.collections
import ase.io
import eventlet  # noqa - eventlet must be installed for flask-socketio to start a production server
import pytest
import redis
import requests
import socketio
import znh5md

from zndraw.start_celery import run_celery_worker


@dataclass
class RoomConnection:
    """Active connection to a room with authentication.

    Attributes
    ----------
    sio
        The connected socket.io client.
    headers
        HTTP headers with Authorization and X-Session-ID for REST calls.
    room_id
        The room identifier.
    session_id
        The session identifier assigned by the server.
    """

    sio: socketio.Client
    headers: dict
    room_id: str
    session_id: str


def _get_free_port() -> int:
    """Get a free port using socket binding.

    This is safer than random port selection as it guarantees the port
    is actually available at the time of binding.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("", 0))
        return sock.getsockname()[1]


def _wait_for_server(port: int, timeout: float = 10.0) -> bool:
    """Wait for server to be ready on the given port.

    Parameters
    ----------
    port
        The port to check for server availability.
    timeout
        Maximum time to wait in seconds.

    Returns
    -------
    bool
        True if server is ready, False if timeout reached.
    """
    start = time.time()
    while time.time() - start < timeout:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                sock.settimeout(0.1)
                sock.connect(("127.0.0.1", port))
                return True
        except (ConnectionRefusedError, OSError):
            time.sleep(0.1)
    return False


def _get_jwt_auth_headers(server_url: str, user_name: str | None = None) -> dict:
    """Register user and get JWT authentication headers for API requests.

    Follows the proper auth flow: register first, then login.

    Args:
        server_url: The base URL of the server (e.g., "http://127.0.0.1:5000")
        user_name: Optional username. If not provided, generates a random one.

    Returns:
        Dictionary with Authorization header containing JWT token
    """
    import requests

    if user_name is None:
        user_name = f"test-user-{uuid.uuid4().hex[:8]}"

    # Step 1: Register user (creates in backend)
    register_response = requests.post(
        f"{server_url}/api/user/register", json={"userName": user_name}, timeout=10
    )
    # 409 = already exists, which is fine (allows reusing usernames in tests)
    if register_response.status_code not in (200, 201, 409):
        raise RuntimeError(f"Registration failed: {register_response.text}")

    # Step 2: Login (get JWT)
    response = requests.post(
        f"{server_url}/api/login", json={"userName": user_name}, timeout=10
    )

    if response.status_code != 200:
        raise RuntimeError(f"Login failed: {response.text}")

    data = response.json()
    return {"Authorization": f"Bearer {data['token']}"}


def _create_room_via_rest(
    server: str, room_id: str, auth_headers: dict, template: str = "none"
) -> None:
    """Create a room via REST API.

    Parameters
    ----------
    server
        The base URL of the server.
    room_id
        The room ID to create.
    auth_headers
        Headers dict containing Authorization header.
    template
        Template name. Defaults to "none" for 0 frames (test isolation).
        Use "empty" for 1 frame with empty ase.Atoms.

    Raises
    ------
    RuntimeError
        If room creation fails (except 409 Conflict which is ignored).
    """
    payload = {"roomId": room_id, "template": template}

    response = requests.post(
        f"{server}/api/rooms",
        json=payload,
        headers=auth_headers,
        timeout=10,
    )
    if response.status_code not in (200, 201, 409):
        raise RuntimeError(f"Failed to create room {room_id}: {response.text}")


def _create_and_join_room(server: str, room: str, auth_headers: dict) -> str:
    """Create a room and join it via socket, returning session_id.

    Creates the room via REST API first (ensures 0 frames), then joins via socket.
    This prevents frontend auto-creation which would apply the "empty" template.

    Parameters
    ----------
    server
        The base URL of the server.
    room
        The room ID to create and join.
    auth_headers
        Headers dict containing Authorization header.

    Returns
    -------
    str
        The session_id from the room:join response.
    """
    import socketio

    jwt_token = auth_headers["Authorization"].replace("Bearer ", "")

    _create_room_via_rest(server, room, auth_headers)

    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    try:
        response = sio.call("room:join", {"roomId": room, "clientType": "frontend"})

        if response.get("status") != "ok":
            raise RuntimeError(f"Failed to join room {room}: {response.get('message')}")

        return response["sessionId"]
    finally:
        sio.disconnect()


def _join_room_and_get_headers(
    server: str, room_id: str, user: str = "test-user"
) -> dict:
    """Join a room and return auth headers with session ID."""
    auth_headers = _get_jwt_auth_headers(server, user)
    session_id = _create_and_join_room(server, room_id, auth_headers)
    return {**auth_headers, "X-Session-ID": session_id}


def _create_room_connection(
    server: str, room_id: str, user: str = "test-user"
) -> RoomConnection:
    """Create a room connection that stays connected.

    Creates the room via REST API first (ensures 0 frames), then joins via socket.
    This prevents frontend auto-creation which would apply the "empty" template.
    The caller is responsible for disconnecting the socket.

    Parameters
    ----------
    server
        The base URL of the server.
    room_id
        The room identifier.
    user
        Username for authentication.

    Returns
    -------
    RoomConnection
        Active connection with socket, headers, and session info.
    """
    auth_headers = _get_jwt_auth_headers(server, user)
    jwt_token = auth_headers["Authorization"].replace("Bearer ", "")

    _create_room_via_rest(server, room_id, auth_headers)

    sio = socketio.Client()
    sio.connect(server, auth={"token": jwt_token}, wait=True)

    response = sio.call("room:join", {"roomId": room_id, "clientType": "frontend"})

    if response.get("status") != "ok":
        sio.disconnect()
        raise RuntimeError(f"Failed to join room {room_id}: {response.get('message')}")

    session_id = response["sessionId"]
    headers = {**auth_headers, "X-Session-ID": session_id}

    return RoomConnection(
        sio=sio,
        headers=headers,
        room_id=room_id,
        session_id=session_id,
    )


@pytest.fixture
def get_free_port():
    """Fixture that provides the get_free_port function."""
    return _get_free_port


@pytest.fixture
def wait_for_server():
    """Fixture that provides the wait_for_server function."""
    return _wait_for_server


@pytest.fixture
def get_jwt_auth_headers():
    """Fixture that provides the get_jwt_auth_headers function."""
    return _get_jwt_auth_headers


@pytest.fixture
def create_and_join_room():
    """Fixture that provides the create_and_join_room function."""
    return _create_and_join_room


@pytest.fixture
def join_room_and_get_headers():
    """Fixture that provides the join_room_and_get_headers function."""
    return _join_room_and_get_headers


@pytest.fixture
def connect_room(server):
    """Factory for room connections with automatic cleanup.

    Creates room connections that stay alive for the test duration.
    All connections are automatically disconnected after the test.

    Usage
    -----
    def test_something(connect_room):
        conn = connect_room("my-room")
        # conn.headers, conn.sio, conn.room_id, conn.session_id available
        # No cleanup needed - fixture handles it
    """
    connections: list[RoomConnection] = []

    def create(room_id: str, user: str = "test-user") -> RoomConnection:
        conn = _create_room_connection(server, room_id, user)
        connections.append(conn)
        return conn

    yield create

    # Automatic cleanup after test
    for conn in connections:
        if conn.sio.connected:
            conn.sio.disconnect()


@pytest.fixture
def app():
    """Create a Flask app for unit testing."""
    from zndraw.config import ZnDrawConfig
    from zndraw.server import create_app

    # Create app with in-memory storage for testing
    config = ZnDrawConfig(redis_url=None)
    test_app = create_app(config=config)
    test_app.config["TESTING"] = True

    yield test_app


@pytest.fixture
def redis_client():
    """Create a Redis client and clean up after test."""
    client = redis.Redis(host="localhost", port=6379, decode_responses=True)
    yield client
    client.flushall()


# Test admin credentials (used in admin mode)
TEST_ADMIN_USERNAME = "test-admin"
TEST_ADMIN_PASSWORD = "test-admin-password"


def _create_server_process(
    tmp_path, admin_mode: bool = False
) -> t.Generator[str, None, None]:
    """Create and manage a zndraw server subprocess.

    Parameters
    ----------
    tmp_path
        Temporary directory for server data and logs.
    admin_mode
        If True, start server with admin credentials configured.
        If False, start in local mode (all users are admin).

    Yields
    ------
    str
        Server URL (e.g., "http://127.0.0.1:5000")
    """
    import os

    port = _get_free_port()
    storage_path = tmp_path / "zndraw-data"
    redis_url = "redis://localhost:6379"

    # Log files for debugging
    server_log = tmp_path / "server.log"
    server_err = tmp_path / "server_err.log"

    # Create environment for subprocess
    env = os.environ.copy()

    if admin_mode:
        # Admin mode: set admin credentials
        env["ZNDRAW_ADMIN_USERNAME"] = TEST_ADMIN_USERNAME
        env["ZNDRAW_ADMIN_PASSWORD"] = TEST_ADMIN_PASSWORD
    else:
        # Local mode: remove any admin credentials
        env.pop("ZNDRAW_ADMIN_USERNAME", None)
        env.pop("ZNDRAW_ADMIN_PASSWORD", None)

    # Start zndraw-server subprocess
    # Using unique port ensures new server starts (no --force-new-server needed)
    with open(server_log, "w") as stdout_f, open(server_err, "w") as stderr_f:
        proc = subprocess.Popen(
            [
                "zndraw",
                "--port",
                str(port),
                "--no-celery",
                "--storage-path",
                str(storage_path),
                "--redis-url",
                redis_url,
                "--no-browser",
            ],
            stdout=stdout_f,
            stderr=stderr_f,
            env=env,
        )

    # Wait for the server to be ready
    if not _wait_for_server(port):
        proc.kill()
        raise TimeoutError("Server did not start in time")

    try:
        yield f"http://127.0.0.1:{port}"
    finally:
        proc.send_signal(signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            raise RuntimeError("Server did not shut down in time")
        finally:
            # Print logs for debugging
            if server_log.exists():
                print("\n=== Server stdout ===")
                print(server_log.read_text())
            if server_err.exists():
                print("\n=== Server stderr ===")
                print(server_err.read_text())

            # Clean up storage and Redis
            shutil.rmtree(storage_path, ignore_errors=True)
            r = redis.Redis.from_url(redis_url, decode_responses=True)
            r.flushall()

            # Clean up server info file
            from zndraw.server_manager import remove_server_info

            remove_server_info(port)


@pytest.fixture
def server(tmp_path) -> t.Generator[str, None, None]:
    """Server in LOCAL mode (all users are admin).

    Use this fixture for tests that don't need to differentiate
    between admin and non-admin users.
    """
    yield from _create_server_process(tmp_path, admin_mode=False)


@pytest.fixture
def server_admin_mode(tmp_path) -> t.Generator[str, None, None]:
    """Server in ADMIN mode (only granted users are admin).

    Use this fixture for tests that need to verify admin vs non-admin behavior.
    Admin credentials: TEST_ADMIN_USERNAME / TEST_ADMIN_PASSWORD

    In admin mode:
    - Regular users are NOT admins by default
    - Only the configured admin user has admin privileges
    - Users can be granted admin via AdminService.grant_admin()
    """
    yield from _create_server_process(tmp_path, admin_mode=True)


@pytest.fixture
def celery_worker():
    from zndraw.config import ZnDrawConfig

    config = ZnDrawConfig(redis_url="redis://localhost:6379")
    worker = run_celery_worker(config)
    try:
        yield worker
    finally:
        worker.terminate()
        try:
            worker.wait(timeout=10)
        except subprocess.TimeoutExpired:
            # Eventlet workers on Python 3.13 may not respond to SIGTERM gracefully
            # Force kill and continue (test cleanup should still succeed)
            worker.kill()
            worker.wait(timeout=5)


@pytest.fixture
def s22() -> list[ase.Atom]:
    """Return a list of 22 atoms."""
    return list(ase.collections.s22)


@pytest.fixture
def s22_xyz(s22, tmp_path) -> str:
    """Return the S22 trajectory as an Atoms object with multiple frames."""
    traj_path = tmp_path / "s22.xyz"
    ase.io.write(traj_path, s22)
    return traj_path.as_posix()


@pytest.fixture
def s22_h5(s22, tmp_path) -> str:
    """Return the S22 trajectory as an Atoms object with multiple frames."""
    traj_path = tmp_path / "s22.h5"
    znh5md.write(traj_path, s22)
    return traj_path.as_posix()


@pytest.fixture
def s22_db(s22, tmp_path) -> str:
    """Return the S22 trajectory as an ASE database with multiple structures."""
    import ase.db

    db_path = tmp_path / "s22.db"
    db = ase.db.connect(str(db_path))

    # Write all S22 structures to the database
    for i, atoms in enumerate(s22):
        db.write(atoms, name=f"s22_{i}", index=i)

    return db_path.as_posix()


@pytest.fixture
def s22_json_db(s22, tmp_path) -> str:
    """Return the S22 trajectory as a JSON-format ASE database."""
    import ase.db

    db_path = tmp_path / "s22.json"
    db = ase.db.connect(str(db_path))

    # Write all S22 structures to the database
    for i, atoms in enumerate(s22):
        db.write(atoms, name=f"s22_{i}", index=i)

    return db_path.as_posix()


@pytest.fixture
def joined_room(server, request):
    """Create a room and return tuple of (server_url, room_name).

    Room name is automatically generated from test function name.
    If test name is "test_foo_bar", room will be "test-foo-bar".

    Note: SessionId is created via socket room:join when needed.
    This fixture just ensures the room exists.

    Example:
        def test_my_feature(joined_room):
            server, room = joined_room
            # Room exists, use _join_room_and_get_headers for sessionId
            headers = _join_room_and_get_headers(server, room)
    """
    import requests

    # Generate room name from test function name
    test_name = request.node.name
    room = test_name.replace("_", "-")

    # Create the room
    auth_headers = _get_jwt_auth_headers(server)
    create_response = requests.post(
        f"{server}/api/rooms",
        json={"roomId": room},
        headers=auth_headers,
        timeout=10,
    )
    assert create_response.status_code in (200, 201, 409), (
        f"Failed to create room {room}"
    )

    return server, room
