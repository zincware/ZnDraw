import random
import shutil
import signal
import socket
import subprocess
import time
import typing as t
from pathlib import Path

import ase.collections
import eventlet  # noqa - eventlet must be installed for flask-socketio to start a production server
import pytest
import redis
import ase.io
import znh5md

from zndraw.start_celery import run_celery_worker


@pytest.fixture
def server(tmp_path) -> t.Generator[str, None, None]:
    port = random.randint(10000, 20000)
    storage_path = tmp_path / "zndraw-data.zarr"
    redis_url = "redis://localhost:6379"

    # Log files for debugging
    server_log = tmp_path / "server.log"
    server_err = tmp_path / "server_err.log"

    # Start zndraw-server subprocess
    with open(server_log, "w") as stdout_f, open(server_err, "w") as stderr_f:
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