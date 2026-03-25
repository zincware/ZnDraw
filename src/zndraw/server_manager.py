"""Server management utilities for ZnDraw CLI.

Retained utilities: process checking, server readiness polling, and shutdown.
Server registry and token storage moved to zndraw.state_file.
"""

from __future__ import annotations

import logging
import os
import time

import httpx

log = logging.getLogger(__name__)

DEFAULT_PORT = 8000


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
    except (OSError, ProcessLookupError):
        return False
    return True


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
                    log.debug(
                        "Server ready after %.2fs (%s attempts)",
                        elapsed,
                        attempt,
                    )
                    return True
            except httpx.RequestError:
                pass

            time.sleep(current_interval)
            current_interval = min(current_interval * 1.5, max_interval)

    log.warning("Server not ready after %ss (%s attempts)", timeout, attempt)
    return False


def shutdown_server(url: str, local_token: str | None = None) -> bool:
    """Shutdown the server gracefully via API endpoint.

    Parameters
    ----------
    url
        Server URL to shut down.
    local_token
        Bearer token for local admin auth.

    Returns
    -------
    bool
        True if the server was successfully shut down, False otherwise.
    """
    headers: dict[str, str] = {}
    if local_token:
        headers["Authorization"] = f"Bearer {local_token}"
    try:
        with httpx.Client(timeout=5.0) as client:
            client.post(f"{url}/v1/admin/shutdown", headers=headers)
    except httpx.RequestError:
        pass

    # Wait for server to stop responding
    for _ in range(25):
        try:
            with httpx.Client(timeout=2.0) as check:
                resp = check.get(f"{url}/v1/health")
                if resp.status_code != 200:
                    return True
        except httpx.RequestError:
            return True
        time.sleep(0.2)
    return False
