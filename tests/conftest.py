import random
import shutil
import signal
import socket
import subprocess
import time
import typing as t
import uuid

import ase.collections
import ase.io
import eventlet  # noqa - eventlet must be installed for flask-socketio to start a production server
import pytest
import redis
import znh5md

from zndraw.start_celery import run_celery_worker


def get_jwt_auth_headers(server_url: str, user_name: str | None = None) -> dict:
    """Login and get JWT authentication headers for API requests.

    Args:
        server_url: The base URL of the server (e.g., "http://127.0.0.1:5000")
        user_name: Optional username. If not provided, generates a random one.

    Returns:
        Dictionary with Authorization header containing JWT token
    """
    import requests

    if user_name is None:
        user_name = f"test-user-{uuid.uuid4().hex[:8]}"

    response = requests.post(f"{server_url}/api/login", json={"userName": user_name})

    if response.status_code != 200:
        raise RuntimeError(f"Login failed: {response.text}")

    data = response.json()
    return {"Authorization": f"Bearer {data['token']}"}


@pytest.fixture
def app():
    """Create a Flask app for unit testing."""
    from zndraw.server import create_app

    # Create app with in-memory storage for testing
    test_app = create_app(storage_path=":memory:", redis_url=None)
    test_app.config["TESTING"] = True

    yield test_app


@pytest.fixture
def redis_client():
    """Create a Redis client and clean up after test."""
    client = redis.Redis(host="localhost", port=6379, decode_responses=True)
    yield client
    client.flushall()


@pytest.fixture
def server(tmp_path) -> t.Generator[str, None, None]:
    port = random.randint(10000, 20000)
    storage_path = tmp_path / "zndraw-data"
    redis_url = "redis://localhost:6379"

    # Log files for debugging
    server_log = tmp_path / "server.log"
    server_err = tmp_path / "server_err.log"

    # Start zndraw-server subprocess
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
        )

    # Wait for the server to be ready
    for _ in range(100):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                sock.settimeout(0.1)
                sock.connect(("127.0.0.1", port))
                break
        except (ConnectionRefusedError, OSError):
            time.sleep(0.1)
    else:
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
            remove_server_info()


@pytest.fixture
def celery_worker():
    redis_url = "redis://localhost:6379"
    worker = run_celery_worker(redis_url)
    try:
        yield worker
    finally:
        worker.terminate()
        try:
            worker.wait(timeout=5)
        except subprocess.TimeoutExpired:
            worker.kill()
            raise RuntimeError("Celery worker did not shut down in time")


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
    """Join a room and return tuple of (server_url, room_name).

    Room name is automatically generated from test function name.
    If test name is "test_foo_bar", room will be "test-foo-bar".

    Example:
        def test_my_feature(joined_room):
            server, room = joined_room
            # Room is already joined, ready to use
            response = requests.get(f"{server}/api/rooms/{room}/geometries")
    """
    import requests

    # Generate room name from test function name
    test_name = request.node.name
    room = test_name.replace("_", "-")

    # Join the room with JWT authentication
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=get_jwt_auth_headers(server)
    )
    assert response.status_code == 200, f"Failed to join room {room}"

    return server, room
