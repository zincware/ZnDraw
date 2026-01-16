"""Server management utilities for ZnDraw CLI."""

import json
import logging
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import requests

log = logging.getLogger(__name__)

DEFAULT_PORT = 5000


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
    def from_dict(cls, data: dict) -> "ServerInfo":
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
        response = requests.get(f"http://localhost:{port}/health", timeout=timeout)
        return response.status_code == 200
    except requests.RequestException:
        return False


def wait_for_server_ready(
    url: str, timeout: float = 30.0, poll_interval: float = 0.2
) -> bool:
    """Wait for a server to become ready by polling the health endpoint.

    Uses exponential backoff starting from poll_interval, capped at 2 seconds.

    Parameters
    ----------
    url
        Base URL of the server (e.g., 'http://localhost:5000').
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

    while time.time() - start_time < timeout:
        attempt += 1
        try:
            response = requests.get(f"{url}/health", timeout=2.0)
            if response.status_code == 200:
                elapsed = time.time() - start_time
                log.debug(f"Server ready after {elapsed:.2f}s ({attempt} attempts)")
                return True
        except requests.RequestException:
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
        f"Server running (PID: {server_info.pid}, Port: {port}, Version: {server_info.version})",
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

        requests.post(
            f"http://localhost:{server_info.port}/api/shutdown",
            headers=headers,
            timeout=5.0,
        )
    except requests.RequestException:
        pass

    # Wait for process to exit
    for _ in range(25):
        if not is_process_running(server_info.pid):
            remove_server_info(server_info.port)
            return True
        time.sleep(0.2)
    return False
