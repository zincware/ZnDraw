"""Server management utilities for ZnDraw CLI."""

import json
import os
import signal
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import requests


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

    """

    pid: int
    port: int
    version: str

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "ServerInfo":
        """Create ServerInfo from dictionary."""
        return cls(**data)


def get_pid_file_path() -> Path:
    """Get the path to the server PID file.

    Returns
    -------
    Path
        Path to ~/.zndraw/server.pid

    """
    zndraw_dir = Path.home() / ".zndraw"
    zndraw_dir.mkdir(exist_ok=True)
    return zndraw_dir / "server.pid"


def read_server_info() -> ServerInfo | None:
    """Read server information from the PID file.

    Returns
    -------
    ServerInfo | None
        Server information if the file exists and is valid, None otherwise.

    """
    pid_file = get_pid_file_path()
    if not pid_file.exists():
        return None

    try:
        with open(pid_file, "r") as f:
            data = json.load(f)
        return ServerInfo.from_dict(data)
    except (json.JSONDecodeError, KeyError, ValueError):
        return None


def write_server_info(server_info: ServerInfo) -> None:
    """Write server information to the PID file.

    Parameters
    ----------
    server_info
        Server information to write.

    """
    pid_file = get_pid_file_path()
    with open(pid_file, "w") as f:
        json.dump(server_info.to_dict(), f)


def remove_server_info() -> None:
    """Remove the server PID file."""
    pid_file = get_pid_file_path()
    if pid_file.exists():
        pid_file.unlink()


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
        # Signal 0 doesn't kill the process, just checks if it exists
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
    except (requests.RequestException, requests.ConnectionError):
        return False


def get_server_status() -> tuple[bool, ServerInfo | None, str]:
    """Get the status of the local server.

    Returns
    -------
    tuple[bool, ServerInfo | None, str]
        A tuple of (is_running, server_info, status_message).

    """
    server_info = read_server_info()

    if server_info is None:
        return False, None, "No server.pid file found"

    if not is_process_running(server_info.pid):
        return False, server_info, f"Process {server_info.pid} is not running"

    if not is_server_responsive(server_info.port):
        return (
            False,
            server_info,
            f"Process {server_info.pid} is running but server on port {server_info.port} is not responsive",
        )

    return (
        True,
        server_info,
        f"Server is running (PID: {server_info.pid}, Port: {server_info.port}, Version: {server_info.version})",
    )


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
        response = requests.post(
            f"http://localhost:{server_info.port}/api/shutdown", timeout=5.0
        )
        return response.status_code == 200
    except (requests.RequestException, requests.ConnectionError):
        # the server shutdown will happen instantly, so we do see an error here!
        for _ in range(25):  # wait for up to 5 seconds
            if not is_process_running(server_info.pid):
                return True
            time.sleep(0.2)
        return False
