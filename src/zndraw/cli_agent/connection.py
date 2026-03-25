"""Connection management for zndraw-cli.

Handles server URL resolution and token resolution via ClientSettings
(pydantic-settings source chain). Provides a Connection dataclass
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

# Shared type aliases for per-subcommand options
UrlOpt = Annotated[
    str | None,
    typer.Option("--url", help="ZnDraw server URL [env: ZNDRAW_URL]."),
]
TokenOpt = Annotated[
    str | None,
    typer.Option("--token", help="Auth token [env: ZNDRAW_TOKEN]."),
]
RoomOpt = Annotated[
    str | None,
    typer.Option("--room", help="Room ID [env: ZNDRAW_ROOM]."),
]
UserOpt = Annotated[
    str | None,
    typer.Option("--user", help="User email [env: ZNDRAW_USER]."),
]
PasswordOpt = Annotated[
    str | None,
    typer.Option("--password", help="Password [env: ZNDRAW_PASSWORD]."),
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


def resolve_room(room: str | None) -> str:
    """Resolve room from ``--room`` flag, env var, or ClientSettings.

    Parameters
    ----------
    room
        Room value from the ``--room`` flag, or None.
    """
    if room is not None:
        return room

    from zndraw.client.settings import ClientSettings

    try:
        settings = ClientSettings()
    except Exception:
        pass
    else:
        if settings.room is not None:
            return settings.room

    die(
        "Room Required",
        "Pass --room or set ZNDRAW_ROOM env var.",
        400,
        EXIT_CLIENT_ERROR,
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


def get_connection(
    url: str | None,
    token: str | None,
    user: str | None = None,
    password: str | None = None,
) -> Connection:
    """Create a Connection using ClientSettings for resolution.

    Parameters
    ----------
    url
        From ``--url`` flag, or None.
    token
        From ``--token`` flag, or None.
    user
        From ``--user`` flag, or None.
    password
        From ``--password`` flag, or None.
    """
    from zndraw.auth_utils import guest_login, login_with_credentials
    from zndraw.client.settings import ClientSettings

    overrides = {
        k: v
        for k, v in {
            "url": url,
            "token": token,
            "user": user,
            "password": password,
        }.items()
        if v is not None
    }

    try:
        settings = ClientSettings(**overrides)
    except Exception as exc:
        die("Configuration Error", str(exc), 400, EXIT_CLIENT_ERROR)

    if settings.url is None:
        die(
            "No Server Found",
            "No running zndraw server found. Start one with `uv run zndraw` or pass `--url`.",
            503,
            EXIT_CONNECTION_ERROR,
        )

    if settings.token is not None:
        resolved_token = settings.token
    elif settings.user and settings.password:
        try:
            resolved_token = login_with_credentials(
                settings.url, settings.user, settings.password
            )
        except (httpx.HTTPStatusError, httpx.RequestError, KeyError) as exc:
            die("Authentication Failed", str(exc), 401, EXIT_CONNECTION_ERROR)
    else:
        try:
            resolved_token = guest_login(settings.url)
        except (httpx.HTTPStatusError, httpx.RequestError, KeyError) as exc:
            die("Authentication Failed", str(exc), 401, EXIT_CONNECTION_ERROR)

    return Connection(base_url=settings.url, token=resolved_token)


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
        From ``--url`` flag, or None.
    token
        From ``--token`` flag, or None.
    room
        Room ID.
    user
        From ``--user`` flag, or None.
    password
        From ``--password`` flag, or None.
    """
    from zndraw import ZnDraw

    return ZnDraw(
        url=url,
        room=room,
        token=token,
        user=user,
        password=password,
        create_if_missing=False,
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
