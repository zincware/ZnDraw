"""Server management utilities for ZnDraw CLI.

This module handles server process discovery, lifecycle management,
PID file management, and token storage for multi-instance support.
"""

from __future__ import annotations

import json
import logging
import os
import stat
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import httpx
from pydantic import BaseModel, TypeAdapter

log = logging.getLogger(__name__)

DEFAULT_PORT = 8000  # FastAPI default


@dataclass
class ServerInfo:
    """Information about a running ZnDraw server.

    Attributes
    ----------
    pid
        Process ID of the server.
    port
        Port number the server is running on.
    version
        Version of ZnDraw the server is running.
    shutdown_token
        Random token for secure shutdown from CLI.
    """

    pid: int
    port: int
    version: str
    shutdown_token: str | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> ServerInfo:
        """Create ServerInfo from dictionary."""
        return cls(**data)


def get_zndraw_dir() -> Path:
    """Get the ZnDraw configuration directory.

    Returns
    -------
    Path
        Path to ~/.zndraw directory (created if not exists).
    """
    zndraw_dir = Path.home() / ".zndraw"
    zndraw_dir.mkdir(exist_ok=True)
    return zndraw_dir


def get_pid_file_path(port: int) -> Path:
    """Get the path to a port-specific server PID file.

    Parameters
    ----------
    port
        Port number for the server.

    Returns
    -------
    Path
        Path to ~/.zndraw/server-{port}.pid
    """
    return get_zndraw_dir() / f"server-{port}.pid"


def read_server_info(port: int) -> ServerInfo | None:
    """Read server information from a port-specific PID file.

    Parameters
    ----------
    port
        Port number to read server info for.

    Returns
    -------
    ServerInfo | None
        Server information if the file exists and is valid, None otherwise.
    """
    pid_file = get_pid_file_path(port)
    if not pid_file.exists():
        return None

    try:
        with open(pid_file) as f:
            data = json.load(f)
        return ServerInfo.from_dict(data)
    except (json.JSONDecodeError, KeyError, ValueError, TypeError):
        return None


def write_server_info(server_info: ServerInfo) -> None:
    """Write server information to the port-specific PID file.

    Parameters
    ----------
    server_info
        Server information to write. Port is derived from server_info.port.
    """
    pid_file = get_pid_file_path(server_info.port)
    with open(pid_file, "w") as f:
        json.dump(server_info.to_dict(), f)


def remove_server_info(port: int) -> None:
    """Remove the port-specific server PID file.

    Parameters
    ----------
    port
        Port number of the server to remove info for.
    """
    pid_file = get_pid_file_path(port)
    if pid_file.exists():
        pid_file.unlink()


def list_all_pid_files() -> list[Path]:
    """List all server PID files in the ZnDraw directory.

    Returns
    -------
    list[Path]
        List of paths to server-{port}.pid files.
    """
    zndraw_dir = get_zndraw_dir()
    return sorted(zndraw_dir.glob("server-*.pid"))


def is_process_running(pid: int) -> bool:
    """Check if a process with the given PID is running.

    Parameters
    ----------
    pid
        Process ID to check.

    Returns
    -------
    bool
        True if the process is running, False otherwise.
    """
    try:
        os.kill(pid, 0)
        return True
    except (OSError, ProcessLookupError):
        return False


def is_server_responsive(port: int, timeout: float = 2.0) -> bool:
    """Check if a server is responsive on the given port.

    Parameters
    ----------
    port
        Port number to check.
    timeout
        Request timeout in seconds.

    Returns
    -------
    bool
        True if the server responds, False otherwise.
    """
    try:
        with httpx.Client(timeout=timeout) as client:
            response = client.get(f"http://localhost:{port}/v1/health")
            return response.status_code == 200
    except httpx.RequestError:
        return False


def wait_for_server_ready(
    url: str, timeout: float = 30.0, poll_interval: float = 0.2
) -> bool:
    """Wait for a server to become ready by polling the health endpoint.

    Uses exponential backoff starting from poll_interval, capped at 2 seconds.

    Parameters
    ----------
    url
        Base URL of the server (e.g., 'http://localhost:8000').
    timeout
        Maximum time to wait in seconds.
    poll_interval
        Initial interval between health checks in seconds.

    Returns
    -------
    bool
        True if the server became ready within the timeout, False otherwise.
    """
    start_time = time.time()
    current_interval = poll_interval
    max_interval = 2.0
    attempt = 0

    with httpx.Client(timeout=2.0) as client:
        while time.time() - start_time < timeout:
            attempt += 1
            try:
                response = client.get(f"{url}/v1/health")
                if response.status_code == 200:
                    elapsed = time.time() - start_time
                    log.debug(f"Server ready after {elapsed:.2f}s ({attempt} attempts)")
                    return True
            except httpx.RequestError:
                pass

            time.sleep(current_interval)
            current_interval = min(current_interval * 1.5, max_interval)

    log.warning(f"Server not ready after {timeout}s ({attempt} attempts)")
    return False


def get_server_status(port: int) -> tuple[bool, ServerInfo | None, str]:
    """Get the status of a server on a specific port.

    Parameters
    ----------
    port
        Port number to check.

    Returns
    -------
    tuple[bool, ServerInfo | None, str]
        A tuple of (is_running, server_info, status_message).
    """
    server_info = read_server_info(port)

    if server_info is None:
        return False, None, f"No PID file for port {port}"

    if not is_process_running(server_info.pid):
        remove_server_info(port)
        return False, None, f"Process {server_info.pid} is not running (cleaned up)"

    if not is_server_responsive(server_info.port):
        return (
            False,
            server_info,
            f"Process {server_info.pid} running but not responsive on port {port}",
        )

    return (
        True,
        server_info,
        f"Server running (PID: {server_info.pid}, Port: {port}, "
        f"Version: {server_info.version})",
    )


def find_running_server(port: int | None = None) -> ServerInfo | None:
    """Find a running server, optionally for a specific port.

    If port is specified, only checks that port.
    If port is None, checks default port first, then smallest port number.

    Parameters
    ----------
    port
        Specific port to check, or None to auto-discover.

    Returns
    -------
    ServerInfo | None
        Server info if a running server is found, None otherwise.
    """
    if port is not None:
        is_running, server_info, _ = get_server_status(port)
        return server_info if is_running else None

    # Check default port first
    is_running, server_info, _ = get_server_status(DEFAULT_PORT)
    if is_running and server_info is not None:
        return server_info

    # Scan all PID files and find running servers
    running_servers: list[ServerInfo] = []
    for pid_file in list_all_pid_files():
        try:
            port_str = pid_file.stem.split("-")[1]
            file_port = int(port_str)
        except (IndexError, ValueError):
            continue

        # Skip default port; it was already checked explicitly above
        if file_port == DEFAULT_PORT:
            continue

        is_running, server_info, _ = get_server_status(file_port)
        if is_running and server_info is not None:
            running_servers.append(server_info)

    if not running_servers:
        return None

    # Return server with smallest port number
    return min(running_servers, key=lambda s: s.port)


def shutdown_server(server_info: ServerInfo) -> bool:
    """Shutdown the server gracefully via API endpoint.

    Parameters
    ----------
    server_info
        Information about the server to shutdown.

    Returns
    -------
    bool
        True if the server was successfully shut down, False otherwise.
    """
    if not is_process_running(server_info.pid):
        return False

    try:
        headers = {}
        if server_info.shutdown_token:
            headers["X-Shutdown-Token"] = server_info.shutdown_token

        with httpx.Client(timeout=5.0) as client:
            client.post(
                f"http://localhost:{server_info.port}/v1/admin/shutdown",
                headers=headers,
            )
    except httpx.RequestError:
        pass

    # Wait for process to exit
    for _ in range(25):
        if not is_process_running(server_info.pid):
            remove_server_info(server_info.port)
            return True
        time.sleep(0.2)
    return False


# --- Token Storage ---

_tokens_adapter = TypeAdapter(dict[str, "TokenEntry"])


class TokenEntry(BaseModel):
    """A stored authentication token for a ZnDraw server."""

    access_token: str
    email: str
    stored_at: datetime


class TokenStore:
    """Persistent token storage in ``tokens.json``.

    Stores JWT tokens keyed by server URL, analogous to browser localStorage.
    File permissions are set to 0600 (owner read/write only).

    Parameters
    ----------
    directory
        Directory to store ``tokens.json`` in. Defaults to ``~/.zndraw``.
    """

    def __init__(self, directory: Path | None = None) -> None:
        self.directory = directory or get_zndraw_dir()

    @property
    def path(self) -> Path:
        """Path to the tokens.json file."""
        return self.directory / "tokens.json"

    def _read_all(self) -> dict[str, TokenEntry]:
        if not self.path.exists():
            return {}
        try:
            raw = self.path.read_text()
            return _tokens_adapter.validate_json(raw)
        except (json.JSONDecodeError, ValueError):
            return {}

    def _write_all(self, data: dict[str, TokenEntry]) -> None:
        raw = _tokens_adapter.dump_json(data, indent=2)
        self.path.write_bytes(raw)
        self.path.chmod(stat.S_IRUSR | stat.S_IWUSR)

    def get(self, server_url: str) -> TokenEntry | None:
        """Get stored token for a server URL."""
        return self._read_all().get(server_url)

    def set(self, server_url: str, entry: TokenEntry) -> None:
        """Store a token for a server URL."""
        data = self._read_all()
        data[server_url] = entry
        self._write_all(data)

    def delete(self, server_url: str) -> None:
        """Remove stored token for a server URL."""
        data = self._read_all()
        if server_url in data:
            del data[server_url]
            self._write_all(data)
