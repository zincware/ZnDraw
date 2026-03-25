"""Custom pydantic-settings source backed by StateFile.

Resolves URL via health-check-based server discovery and token
via local_token (localhost) or stored access_token (remote).
"""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING, Any
from urllib.parse import urlparse

import httpx
from pydantic_settings import BaseSettings, PydanticBaseSettingsSource

from zndraw.state_file import StateFile

if TYPE_CHECKING:
    from pydantic.fields import FieldInfo

    from zndraw.state_file import ServerEntry

log = logging.getLogger(__name__)


def _is_pid_alive(pid: int) -> bool:
    """Check if a process is alive via os.kill(pid, 0).

    Parameters
    ----------
    pid : int
        Process ID to check.

    Returns
    -------
    bool
        True if the process is alive.
    """
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    return True


def _is_url_healthy(url: str, timeout: float = 2.0) -> bool:
    """Check if a server is healthy via GET /v1/health.

    Parameters
    ----------
    url : str
        Server URL to check.
    timeout : float
        Request timeout in seconds.

    Returns
    -------
    bool
        True if the server responds with 200.
    """
    try:
        with httpx.Client(timeout=timeout) as client:
            resp = client.get(f"{url}/v1/health")
            return resp.status_code == 200
    except httpx.RequestError:
        return False


def _is_localhost(url: str) -> bool:
    """Check if a URL points to localhost.

    Parameters
    ----------
    url : str
        URL to check.

    Returns
    -------
    bool
        True if the URL's hostname is localhost, 127.0.0.1, or ::1.
    """
    hostname = urlparse(url).hostname
    return hostname in {"localhost", "127.0.0.1", "::1"}


class StateFileSource(PydanticBaseSettingsSource):
    """Pydantic-settings source backed by ~/.zndraw/state.json.

    Resolves ``url`` via health-check-based server discovery
    (localhost preferred, sorted by last_used descending) and
    ``token`` via local_token (local servers) or stored access_token
    (remote servers).

    Parameters
    ----------
    settings_cls : type[BaseSettings]
        The settings class being configured.
    state_file : StateFile | None
        StateFile instance (defaults to ~/.zndraw).
    """

    def __init__(
        self,
        settings_cls: type[BaseSettings],
        state_file: StateFile | None = None,
    ) -> None:
        super().__init__(settings_cls)
        self._state_file = state_file or StateFile()
        self._current_state: dict[str, Any] = {}

    def get_field_value(
        self,
        field: FieldInfo,  # noqa: ARG002
        field_name: str,
    ) -> tuple[Any, str, bool]:
        """Not used — __call__ handles everything."""
        return None, field_name, False

    def __call__(self) -> dict[str, Any]:
        """Resolve URL and token from state.json.

        Returns
        -------
        dict[str, Any]
            Resolved settings (url, token) or empty dict.
        """
        result: dict[str, Any] = {}

        # Check if URL was already resolved by a higher-priority source
        url_from_above = self.current_state.get("url") or self._current_state.get("url")

        if url_from_above is None:
            # Discover URL from state.json
            url = self._discover_url()
            if url is not None:
                result["url"] = url
        else:
            url = url_from_above

        # Resolve token for the determined URL
        if url is not None:
            token = self._resolve_token(url)
            if token is not None:
                result["token"] = token

        return result

    def _discover_url(self) -> str | None:
        """Discover a healthy server URL from state.json.

        Returns
        -------
        str | None
            The first healthy server URL, or None.
        """
        data = self._state_file.read()
        if not data.servers:
            return None

        local: list[tuple[str, ServerEntry]] = []
        remote: list[tuple[str, ServerEntry]] = []

        for url, entry in data.servers.items():
            if _is_localhost(url):
                local.append((url, entry))
            else:
                remote.append((url, entry))

        local.sort(key=lambda x: x[1].last_used, reverse=True)
        remote.sort(key=lambda x: x[1].last_used, reverse=True)

        for url, entry in local:
            if entry.pid is not None and not _is_pid_alive(entry.pid):
                log.debug("Removing dead local server: %s (PID %d)", url, entry.pid)
                self._state_file.remove_server(url)
                continue
            if _is_url_healthy(url):
                return url

        for url, _entry in remote:
            if _is_url_healthy(url):
                return url
            log.warning(
                "Server %s is unreachable. To remove: zndraw-cli auth logout --url %s",
                url,
                url,
            )

        return None

    def _resolve_token(self, url: str) -> str | None:
        """Resolve token for a URL.

        Parameters
        ----------
        url : str
            Server URL to resolve token for.

        Returns
        -------
        str | None
            The resolved token, or None.
        """
        data = self._state_file.read()
        server = data.servers.get(url)
        if server is not None and _is_localhost(url) and server.local_token:
            return server.local_token
        token_entry = data.tokens.get(url)
        if token_entry is not None:
            return token_entry.access_token
        return None
