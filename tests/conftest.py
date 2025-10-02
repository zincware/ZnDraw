import random
import shutil
import signal
import socket
import subprocess
import time

import ase.collections
import eventlet  # noqa - eventlet must be installed for flask-socketio to start a production server
import pytest
import redis

from zndraw.start_celery import run_celery_worker


@pytest.fixture
def server(tmp_path):
    port = random.randint(10000, 20000)
    storage_path = tmp_path / "zndraw-data.zarr"
    redis_url = "redis://localhost:6379"

    # Start zndraw-server subprocess
    proc = subprocess.Popen(
        [
            "zndraw-server",
            "--port",
            str(port),
            "--no-celery",
            "--storage-path",
            str(storage_path),
            "--redis-url",
            redis_url,
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
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
            # Clean up storage and Redis
            shutil.rmtree(storage_path, ignore_errors=True)
            r = redis.Redis.from_url(redis_url, decode_responses=True)
            r.flushall()


@pytest.fixture
def celery_worker():
    worker = run_celery_worker()
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
