"""Connection management for zndraw-cli.

Handles server URL resolution (flag -> env -> PID file) and token
resolution (flag -> env -> guest auth). Provides a Connection dataclass
that wraps httpx.Client with base_url and auth header.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field
from typing import Any, NoReturn

import httpx

from zndraw.server_manager import find_running_server

# Exit codes
EXIT_OK = 0
EXIT_CLIENT_ERROR = 1
EXIT_SERVER_ERROR = 2
EXIT_CONNECTION_ERROR = 3


def _problem_json(title: str, detail: str, status: int) -> str:
    """Format an RFC 9457 problem JSON string."""
    return json.dumps(
        {
            "type": "about:blank",
            "title": title,
            "status": status,
            "detail": detail,
        }
    )


def die(title: str, detail: str, status: int, exit_code: int) -> NoReturn:
    """Print a problem JSON error to stderr and exit."""
    sys.stderr.write(_problem_json(title, detail, status) + "\n")
    raise SystemExit(exit_code)


@dataclass
class Connection:
    """REST-only connection to a zndraw server.

    Attributes
    ----------
    base_url
        Server base URL, e.g. ``http://localhost:8000``.
    token
        JWT bearer token.
    client
        Pre-configured httpx client.
    """

    base_url: str
    token: str
    client: httpx.Client = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.client = httpx.Client(
            base_url=self.base_url,
            headers={"Authorization": f"Bearer {self.token}"},
            timeout=30.0,
        )

    def request(
        self,
        method: str,
        path: str,
        **kwargs: Any,
    ) -> httpx.Response:
        """Make an HTTP request with structured error handling.

        Parameters
        ----------
        method
            HTTP method (GET, POST, PUT, PATCH, DELETE).
        path
            URL path relative to base_url.
        **kwargs
            Passed to ``httpx.Client.request``.

        Returns
        -------
        httpx.Response
            The response object on success (2xx/3xx).

        Raises
        ------
        SystemExit
            On 4xx/5xx errors, prints problem JSON to stderr and exits.
        """
        try:
            resp = self.client.request(method, path, **kwargs)  # type: ignore[arg-type]
        except httpx.RequestError as exc:
            die("Connection Error", str(exc), 503, EXIT_CONNECTION_ERROR)

        if resp.status_code >= 400:
            # Try to pass through server's problem JSON
            content_type = resp.headers.get("content-type", "")
            if "problem+json" in content_type or "application/json" in content_type:
                sys.stderr.write(resp.text + "\n")
            else:
                sys.stderr.write(
                    _problem_json(
                        resp.reason_phrase or "Error",
                        resp.text[:500],
                        resp.status_code,
                    )
                    + "\n"
                )
            exit_code = (
                EXIT_CLIENT_ERROR if resp.status_code < 500 else EXIT_SERVER_ERROR
            )
            raise SystemExit(exit_code)

        return resp

    def get(self, path: str, **kwargs: Any) -> httpx.Response:
        """GET request."""
        return self.request("GET", path, **kwargs)

    def post(self, path: str, **kwargs: Any) -> httpx.Response:
        """POST request."""
        return self.request("POST", path, **kwargs)

    def put(self, path: str, **kwargs: Any) -> httpx.Response:
        """PUT request."""
        return self.request("PUT", path, **kwargs)

    def patch(self, path: str, **kwargs: Any) -> httpx.Response:
        """PATCH request."""
        return self.request("PATCH", path, **kwargs)

    def delete(self, path: str, **kwargs: Any) -> httpx.Response:
        """DELETE request."""
        return self.request("DELETE", path, **kwargs)

    def close(self) -> None:
        """Close the underlying HTTP client."""
        self.client.close()


def resolve_url(url: str | None) -> str:
    """Resolve the server URL from flag, env, or PID file.

    Parameters
    ----------
    url
        Explicit URL from ``--url`` flag or ``ZNDRAW_URL`` env var.
        ``None`` triggers PID file auto-discovery.
    """
    if url is not None:
        return url.rstrip("/")

    server_info = find_running_server()
    if server_info is not None:
        return f"http://localhost:{server_info.port}"

    die(
        "No Server Found",
        "No running zndraw server found. Start one with `zndraw` or pass `--url`.",
        503,
        EXIT_CONNECTION_ERROR,
    )


def resolve_token(base_url: str, token: str | None) -> str:
    """Resolve the auth token from flag, env, or guest auth.

    Parameters
    ----------
    base_url
        Server URL for guest auth fallback.
    token
        Explicit token from ``--token`` flag or ``ZNDRAW_TOKEN`` env var.
        ``None`` triggers guest auth.
    """
    if token is not None:
        return token

    # Auto-create guest session
    try:
        with httpx.Client(base_url=base_url, timeout=10.0) as client:
            resp = client.post("/v1/auth/guest")
            resp.raise_for_status()
            data = resp.json()
            return data["access_token"]
    except (httpx.RequestError, httpx.HTTPStatusError, KeyError) as exc:
        die(
            "Authentication Failed",
            f"Failed to create guest session: {exc}",
            401,
            EXIT_CONNECTION_ERROR,
        )


def get_current_step(conn: Connection, room: str) -> int:
    """Fetch the current step (frame index) for a room.

    Parameters
    ----------
    conn
        Active connection.
    room
        Room ID.
    """
    resp = conn.get(f"/v1/rooms/{room}/step")
    return resp.json()["step"]


def get_connection(url: str | None, token: str | None) -> Connection:
    """Create a Connection from resolved URL and token.

    Parameters
    ----------
    url
        From ``--url`` flag / ``ZNDRAW_URL`` env var, or None.
    token
        From ``--token`` flag / ``ZNDRAW_TOKEN`` env var, or None.
    """
    base_url = resolve_url(url)
    resolved_token = resolve_token(base_url, token)
    return Connection(base_url=base_url, token=resolved_token)
