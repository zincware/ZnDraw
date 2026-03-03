"""Server lifecycle management for restart scenarios in tests.

This module provides a ServerProvider class that allows starting, stopping,
and restarting ZnDraw servers mid-test. Use this ONLY when you need restart
capability. For normal tests, use the `server` or `server_admin_mode` fixtures.
"""

import os
import socket
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ServerProvider:
    """Manages a ZnDraw server with restart capability.

    Use this ONLY when you need to restart the server mid-test.
    For normal tests, use the `server` or `server_admin_mode` fixtures.

    Parameters
    ----------
    port
        Server port number.
    storage_path
        Path to LMDB storage directory.
    redis_url
        Redis connection URL.
    admin_user
        Admin username (enables deployment mode if set with password).
    admin_password
        Admin password (enables deployment mode if set with username).

    Example
    -------
    def test_extension_persistence(server_provider):
        server = server_provider
        vis = ZnDraw(url=server.url, room="test")
        vis.register_extension(MyExt, public=True)

        server.restart()

        # Verify extension re-registers after restart
        schema = requests.get(f"{server.url}/api/schema/modifiers").json()
        assert "MyExt" in schema
    """

    port: int
    storage_path: Path
    redis_url: str
    admin_user: str | None = None
    admin_password: str | None = None

    # Internal state (not part of dataclass comparison)
    _process: subprocess.Popen | None = field(default=None, repr=False, compare=False)

    @property
    def url(self) -> str:
        """Server URL."""
        return f"http://127.0.0.1:{self.port}"

    @property
    def is_running(self) -> bool:
        """Check if server process is running and responsive."""
        if self._process is None or self._process.poll() is not None:
            return False
        try:
            import requests

            resp = requests.get(f"{self.url}/health", timeout=1)
            return resp.status_code == 200
        except Exception:
            return False

    def start(self) -> str:
        """Start the server.

        Returns
        -------
        str
            Server URL.

        Raises
        ------
        RuntimeError
            If server is already running or fails to start.
        """
        if self._process is not None and self._process.poll() is None:
            raise RuntimeError("Server already running")

        env = os.environ.copy()
        if self.admin_user and self.admin_password:
            env["ZNDRAW_ADMIN_USERNAME"] = self.admin_user
            env["ZNDRAW_ADMIN_PASSWORD"] = self.admin_password
        else:
            env.pop("ZNDRAW_ADMIN_USERNAME", None)
            env.pop("ZNDRAW_ADMIN_PASSWORD", None)

        self._process = subprocess.Popen(
            [
                "zndraw",
                "--port",
                str(self.port),
                "--no-celery",
                "--storage-path",
                str(self.storage_path),
                "--redis-url",
                self.redis_url,
                "--no-browser",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
        )

        if not _wait_for_server(self.port):
            self.kill()
            raise RuntimeError(f"Server failed to start on port {self.port}")

        return self.url

    def stop(self) -> None:
        """Graceful shutdown (SIGTERM)."""
        if self._process is None:
            return
        self._process.terminate()
        try:
            self._process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            self._process.kill()
            self._process.wait()
        self._process = None

    def kill(self) -> None:
        """Force kill (SIGKILL)."""
        if self._process is None:
            return
        self._process.kill()
        self._process.wait()
        self._process = None

    def restart(self) -> str:
        """Stop and start server.

        Returns
        -------
        str
            Server URL.
        """
        self.stop()
        return self.start()

    def flush_redis(self) -> None:
        """Flush all data from Redis.

        Clears all keys in the Redis database, including:
        - Public extensions
        - Room data
        - Session data
        """
        import redis

        client = redis.Redis.from_url(self.redis_url, decode_responses=True)
        client.flushall()

    def fresh_restart(self) -> str:
        """Stop server, flush Redis, and start fresh.

        Use this when you need a completely clean slate with no
        persisted data (extensions, rooms, etc.).

        Returns
        -------
        str
            Server URL.
        """
        self.stop()
        self.flush_redis()
        return self.start()


def _wait_for_server(port: int, timeout: float = 30.0) -> bool:
    """Wait for server to become ready.

    Parameters
    ----------
    port
        Port to check.
    timeout
        Maximum wait time in seconds.

    Returns
    -------
    bool
        True if server is ready, False if timeout.
    """
    start = time.time()
    while time.time() - start < timeout:
        try:
            with socket.socket() as sock:
                sock.settimeout(0.1)
                sock.connect(("127.0.0.1", port))
                return True
        except (ConnectionRefusedError, OSError):
            time.sleep(0.1)
    return False
