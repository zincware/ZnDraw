"""Unified state file for server registry and token storage.

Replaces separate PID files (~/.zndraw/server-{port}.pid) and
tokens.json with a single ~/.zndraw/state.json (mode 0600).

All writes use tempfile + os.rename for atomicity.
"""

from __future__ import annotations

import contextlib
import json
import logging
import os
import stat
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from pydantic import BaseModel, TypeAdapter

log = logging.getLogger(__name__)


class ServerEntry(BaseModel):
    """Registry entry for a known ZnDraw server.

    Attributes
    ----------
    added_at : datetime
        When this server was first registered.
    last_used : datetime
        When this server was last successfully connected to.
    pid : int | None
        Process ID (local servers only).
    version : str | None
        ZnDraw version (local servers only).
    local_token : str | None
        Per-start superuser token (local servers only).
    """

    added_at: datetime
    last_used: datetime
    pid: int | None = None
    version: str | None = None
    local_token: str | None = None
    access_token: str | None = None


class TokenEntry(BaseModel):
    """Stored authentication token for a ZnDraw server.

    Attributes
    ----------
    access_token : str
        JWT access token.
    email : str
        Email address of the authenticated user.
    stored_at : datetime
        When the token was stored.
    """

    access_token: str
    email: str
    stored_at: datetime


class StateData(BaseModel):
    """Root schema for state.json."""

    servers: dict[str, ServerEntry] = {}
    tokens: dict[str, TokenEntry] = {}


_state_adapter = TypeAdapter(StateData)


class StateFile:
    """Manages ~/.zndraw/state.json with atomic writes.

    Parameters
    ----------
    directory : Path | None
        Directory containing state.json. Defaults to ~/.zndraw.
    """

    def __init__(self, directory: Path | None = None) -> None:
        self.directory = directory or (Path.home() / ".zndraw")
        self.directory.mkdir(exist_ok=True)

    @property
    def path(self) -> Path:
        """Path to state.json."""
        return self.directory / "state.json"

    def read(self) -> StateData:
        """Read state.json, returning empty StateData if missing or corrupt."""
        if not self.path.exists():
            return StateData()
        try:
            raw = self.path.read_text()
            return _state_adapter.validate_json(raw)
        except (json.JSONDecodeError, ValueError):
            log.warning("Corrupt state.json at %s, starting fresh", self.path)
            return StateData()

    def _write(self, data: StateData) -> None:
        """Atomic write: tempfile + os.replace, mode 0600."""
        raw = data.model_dump_json(indent=2)
        fd, tmp_path = tempfile.mkstemp(dir=self.directory, suffix=".tmp")
        try:
            with os.fdopen(fd, "wb") as f:
                fd = -1
                f.write(raw.encode())
                f.flush()
                os.fchmod(f.fileno(), stat.S_IRUSR | stat.S_IWUSR)
                os.fsync(f.fileno())
            Path(tmp_path).replace(self.path)
        except BaseException:
            if fd >= 0:
                os.close(fd)
            with contextlib.suppress(OSError):
                Path(tmp_path).unlink()
            raise

    # --- Server registry ---

    def add_server(self, url: str, entry: ServerEntry) -> None:
        """Register or update a server entry.

        Parameters
        ----------
        url : str
            Server URL (used as key).
        entry : ServerEntry
            Server metadata to store.
        """
        data = self.read()
        data.servers[url] = entry
        self._write(data)

    def remove_server(self, url: str) -> None:
        """Remove a server entry (no-op if absent).

        Parameters
        ----------
        url : str
            Server URL to remove.
        """
        data = self.read()
        if url in data.servers:
            del data.servers[url]
            self._write(data)

    def get_server(self, url: str) -> ServerEntry | None:
        """Get a server entry by URL, or None.

        Parameters
        ----------
        url : str
            Server URL to look up.

        Returns
        -------
        ServerEntry | None
            The server entry, or None if not found.
        """
        return self.read().servers.get(url)

    def update_last_used(self, url: str) -> None:
        """Update the last_used timestamp for a server (no-op if absent).

        Parameters
        ----------
        url : str
            Server URL to update.
        """
        data = self.read()
        if url in data.servers:
            data.servers[url].last_used = datetime.now(UTC)
            self._write(data)

    # --- Token storage ---

    def add_token(self, url: str, entry: TokenEntry) -> None:
        """Store an authentication token for a server URL.

        Parameters
        ----------
        url : str
            Server URL (used as key).
        entry : TokenEntry
            Token data to store.
        """
        data = self.read()
        data.tokens[url] = entry
        self._write(data)

    def remove_token(self, url: str) -> None:
        """Remove a stored token (no-op if absent).

        Parameters
        ----------
        url : str
            Server URL whose token should be removed.
        """
        data = self.read()
        if url in data.tokens:
            del data.tokens[url]
            self._write(data)

    def get_token(self, url: str) -> TokenEntry | None:
        """Get stored token for a server URL, or None.

        Parameters
        ----------
        url : str
            Server URL to look up.

        Returns
        -------
        TokenEntry | None
            The token entry, or None if not found.
        """
        return self.read().tokens.get(url)
