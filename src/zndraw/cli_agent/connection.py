"""Connection management for zndraw-cli.

Handles server URL resolution (flag -> env -> PID file) and token
resolution (flag -> env -> guest auth). Provides a Connection dataclass
that wraps httpx.Client with base_url and auth header.
"""

from __future__ import annotations

import contextlib
import json
import sys
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Annotated, Any, NoReturn

import httpx
import typer

if TYPE_CHECKING:
    from collections.abc import Generator

    from zndraw import ZnDraw

from zndraw.server_manager import find_running_server

# Shared type aliases for per-subcommand options
UrlOpt = Annotated[
    str | None,
    typer.Option("--url", envvar="ZNDRAW_URL", help="ZnDraw server URL"),
]
TokenOpt = Annotated[
    str | None,
    typer.Option("--token", envvar="ZNDRAW_TOKEN", help="Auth token"),
]
RoomOpt = Annotated[
    str | None,
    typer.Option("--room", envvar="ZNDRAW_ROOM", help="Room ID"),
]
UserOpt = Annotated[
    str | None,
    typer.Option("--user", envvar="ZNDRAW_USER", help="User email for authentication"),
]
PasswordOpt = Annotated[
    str | None,
    typer.Option(
        "--password", envvar="ZNDRAW_PASSWORD", help="Password for authentication"
    ),
]


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


def resolve_token(
    base_url: str,
    token: str | None,
    user: str | None = None,
    password: str | None = None,
) -> str:
    """Resolve auth token — CLI wrapper around shared auth_utils.

    Parameters
    ----------
    base_url
        Server URL.
    token
        Explicit token from ``--token`` / ``ZNDRAW_TOKEN``.
    user
        User email from ``--user`` / ``ZNDRAW_USER``.
    password
        Password from ``--password`` / ``ZNDRAW_PASSWORD``.
    """
    from zndraw.auth_utils import resolve_token as _resolve_token

    try:
        return _resolve_token(base_url, token=token, user=user, password=password)
    except ValueError as exc:
        die(str(exc), str(exc), 400, EXIT_CLIENT_ERROR)
    except httpx.HTTPStatusError as exc:
        status = exc.response.status_code
        title = "Authentication Failed" if status < 500 else "Server Error"
        exit_code = EXIT_CLIENT_ERROR if status < 500 else EXIT_SERVER_ERROR
        die(title, str(exc), status, exit_code)
    except (httpx.RequestError, KeyError) as exc:
        die(
            "Authentication Failed",
            f"Failed to authenticate: {exc}",
            401,
            EXIT_CONNECTION_ERROR,
        )


def resolve_room(room: str | None) -> str:
    """Resolve room from ``--room`` option or ``ZNDRAW_ROOM`` env var.

    Parameters
    ----------
    room
        Room value from the ``--room`` option (includes env var fallback
        via Typer's ``envvar``), or None.
    """
    if room is None:
        die(
            "Room Required",
            "Pass --room or set ZNDRAW_ROOM env var.",
            400,
            EXIT_CLIENT_ERROR,
        )
    return room


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


def get_connection(
    url: str | None,
    token: str | None,
    user: str | None = None,
    password: str | None = None,
) -> Connection:
    """Create a Connection from resolved URL and token.

    Parameters
    ----------
    url
        From ``--url`` flag / ``ZNDRAW_URL`` env var, or None.
    token
        From ``--token`` flag / ``ZNDRAW_TOKEN`` env var, or None.
    user
        From ``--user`` flag / ``ZNDRAW_USER`` env var, or None.
    password
        From ``--password`` flag / ``ZNDRAW_PASSWORD`` env var, or None.
    """
    base_url = resolve_url(url)
    resolved_token = resolve_token(base_url, token, user=user, password=password)
    return Connection(base_url=base_url, token=resolved_token)


def get_zndraw(
    url: str | None,
    token: str | None,
    room: str,
    user: str | None = None,
    password: str | None = None,
) -> ZnDraw:
    """Create a ZnDraw instance from CLI context.

    Parameters
    ----------
    url
        From ``--url`` flag / ``ZNDRAW_URL`` env var, or None.
    token
        From ``--token`` flag / ``ZNDRAW_TOKEN`` env var, or None.
    room
        Room ID.
    user
        From ``--user`` flag / ``ZNDRAW_USER`` env var, or None.
    password
        From ``--password`` flag / ``ZNDRAW_PASSWORD`` env var, or None.
    """
    from zndraw import ZnDraw

    base_url = resolve_url(url)
    resolved_token = resolve_token(base_url, token, user=user, password=password)
    return ZnDraw(
        url=base_url, room=room, token=resolved_token, create_if_missing=False
    )


@contextlib.contextmanager
def cli_error_handler() -> Generator[None, None, None]:
    """Context manager for CLI error handling.

    Catches httpx errors and Python exceptions raised by
    ``APIManager.raise_for_status`` (KeyError, PermissionError,
    ValueError, ZnDrawError, RoomLockedError) and prints RFC 9457
    problem JSON to stderr.
    """
    try:
        yield
    except httpx.HTTPStatusError as exc:
        content_type = exc.response.headers.get("content-type", "")
        if "problem+json" in content_type or "application/json" in content_type:
            sys.stderr.write(exc.response.text + "\n")
        else:
            sys.stderr.write(
                _problem_json(
                    exc.response.reason_phrase or "Error",
                    exc.response.text[:500],
                    exc.response.status_code,
                )
                + "\n"
            )
        exit_code = (
            EXIT_CLIENT_ERROR if exc.response.status_code < 500 else EXIT_SERVER_ERROR
        )
        raise SystemExit(exit_code) from exc
    except httpx.RequestError as exc:
        die("Connection Error", str(exc), 503, EXIT_CONNECTION_ERROR)
    except KeyError as exc:
        die("Not Found", str(exc), 404, EXIT_CLIENT_ERROR)
    except IndexError as exc:
        die("Not Found", str(exc), 404, EXIT_CLIENT_ERROR)
    except PermissionError as exc:
        die("Forbidden", str(exc), 403, EXIT_CLIENT_ERROR)
    except ValueError as exc:
        die("Unprocessable Entity", str(exc), 422, EXIT_CLIENT_ERROR)
    except _zndraw_exceptions() as exc:
        from zndraw.client import RoomLockedError

        if isinstance(exc, RoomLockedError):
            die("Room Locked", str(exc), 423, EXIT_CLIENT_ERROR)
        die("Server Error", str(exc), 500, EXIT_SERVER_ERROR)


def _zndraw_exceptions() -> tuple[type[Exception], ...]:
    """Lazily import ZnDraw exception types."""
    from zndraw.client import RoomLockedError, ZnDrawError

    return (ZnDrawError, RoomLockedError)
