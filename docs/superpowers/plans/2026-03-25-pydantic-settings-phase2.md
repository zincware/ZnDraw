# Pydantic-Settings Phase 2: Client/CLI Unification

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unify client-side settings resolution via `ClientSettings(BaseSettings)`, a unified `StateFile` replacing PID files + tokens.json, and local admin token support — so that `ZnDraw()` with zero args connects to a local server as admin.

**Architecture:** Build three new modules (`state_file.py`, `settings_sources.py`, `client/settings.py`) that compose into a single settings chain: init > env > pyproject.toml > state file > guest fallback. Then rewire `ZnDraw`, CLI, and auth to use the new chain. Old code (ServerInfo, PID files, TokenStore, manual resolve functions) is removed.

**Tech Stack:** pydantic-settings (custom `PydanticBaseSettingsSource`), pydantic `BaseModel`, httpx (health checks), `os.rename` (atomic writes).

**Spec:** `docs/superpowers/specs/2026-03-25-pydantic-settings-unification-design.md` (Phase 2 section)

---

## File Structure

### New Files

| File | Responsibility |
|---|---|
| `src/zndraw/state_file.py` | `StateFile` class — unified `~/.zndraw/state.json` with `ServerEntry`, `TokenEntry` models, atomic writes, migration from old format |
| `src/zndraw/settings_sources.py` | `StateFileSource(PydanticBaseSettingsSource)` — URL discovery via health checks, token resolution via local_token/stored token |
| `src/zndraw/client/settings.py` | `ClientSettings(BaseSettings)` — env prefix `ZNDRAW_`, pyproject `[tool.zndraw]`, full source chain |
| `tests/test_state_file.py` | StateFile unit tests |
| `tests/test_state_file_source.py` | StateFileSource unit tests |
| `tests/test_client_settings.py` | ClientSettings integration tests |

### Modified Files

| File | Change |
|---|---|
| `src/zndraw/dependencies.py` | Add `get_local_token_user` dependency, `LocalTokenOrAdminDep` |
| `src/zndraw/routes/admin.py` | Shutdown endpoint uses `LocalTokenOrAdminDep` |
| `src/zndraw/client/core.py` | `__post_init__` uses `ClientSettings`, remove `_resolve_url` |
| `src/zndraw/auth_utils.py` | Simplify to `validate_credentials`, `login_with_credentials`, `guest_login` |
| `src/zndraw/cli.py` | Server start writes `state.json`, stop removes entry, shutdown uses `local_token` |
| `src/zndraw/cli_agent/connection.py` | Remove `envvar=` from type aliases, replace `resolve_url`/`resolve_token` with `ClientSettings` |
| `src/zndraw/cli_agent/auth.py` | `login`/`logout` use `StateFile` instead of `TokenStore` |
| `src/zndraw/server_manager.py` | Remove `ServerInfo`, PID functions, `TokenStore`. Keep `wait_for_server_ready`, `is_process_running`, `is_server_responsive` |
| `tests/test_token_store.py` | Delete (replaced by `test_state_file.py`) |
| `tests/test_resolve_token.py` | Update for simplified `auth_utils` |

---

## Task 1: StateFile — Models, CRUD, Atomic Writes

**Files:**
- Create: `src/zndraw/state_file.py`
- Test: `tests/test_state_file.py`

- [ ] **Step 1: Write failing tests for StateFile**

```python
# tests/test_state_file.py
"""Tests for unified StateFile (replaces PID files + tokens.json)."""
from __future__ import annotations

import json
import os
import stat
from datetime import UTC, datetime

import pytest

from zndraw.state_file import ServerEntry, StateFile, TokenEntry


@pytest.fixture
def state_file(tmp_path):
    """StateFile backed by a temporary directory."""
    return StateFile(directory=tmp_path)


@pytest.fixture
def local_server_entry():
    """A local server entry with PID and local_token."""
    return ServerEntry(
        added_at=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
        last_used=datetime(2026, 3, 25, 14, 0, tzinfo=UTC),
        pid=12345,
        version="0.5.0",
        local_token="test-local-token",
    )


@pytest.fixture
def remote_server_entry():
    """A remote server entry (no PID, no local_token)."""
    return ServerEntry(
        added_at=datetime(2026, 3, 24, 9, 0, tzinfo=UTC),
        last_used=datetime(2026, 3, 25, 12, 0, tzinfo=UTC),
    )


@pytest.fixture
def sample_token():
    """A stored authentication token."""
    return TokenEntry(
        access_token="eyJ.test.token",
        email="user@example.com",
        stored_at=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
    )


# --- ServerEntry model ---


def test_server_entry_local_fields(local_server_entry):
    assert local_server_entry.pid == 12345
    assert local_server_entry.local_token == "test-local-token"
    assert local_server_entry.version == "0.5.0"


def test_server_entry_remote_defaults(remote_server_entry):
    assert remote_server_entry.pid is None
    assert remote_server_entry.local_token is None
    assert remote_server_entry.version is None


def test_server_entry_roundtrip(local_server_entry):
    data = local_server_entry.model_dump_json()
    restored = ServerEntry.model_validate_json(data)
    assert restored == local_server_entry


# --- TokenEntry model ---


def test_token_entry_fields(sample_token):
    assert sample_token.access_token == "eyJ.test.token"
    assert sample_token.email == "user@example.com"


def test_token_entry_roundtrip(sample_token):
    data = sample_token.model_dump_json()
    restored = TokenEntry.model_validate_json(data)
    assert restored == sample_token


# --- StateFile read/write ---


def test_read_empty(state_file):
    """Returns empty StateData when state.json doesn't exist."""
    data = state_file.read()
    assert data.servers == {}
    assert data.tokens == {}


def test_read_corrupt_json(state_file):
    """Returns empty StateData for corrupt JSON."""
    (state_file.directory / "state.json").write_text("not valid json{{{")
    data = state_file.read()
    assert data.servers == {}
    assert data.tokens == {}


def test_file_permissions(state_file, local_server_entry):
    """state.json must be mode 0600 (owner read/write only)."""
    state_file.add_server("http://localhost:8000", local_server_entry)
    mode = stat.S_IMODE(state_file.path.stat().st_mode)
    assert mode == 0o600


def test_atomic_write_no_partial(state_file, local_server_entry):
    """state.json should not contain partial data if write is interrupted."""
    state_file.add_server("http://localhost:8000", local_server_entry)
    # Verify file exists and is valid JSON
    raw = state_file.path.read_text()
    data = json.loads(raw)
    assert "servers" in data
    assert "http://localhost:8000" in data["servers"]


# --- Server CRUD ---


def test_add_and_get_server(state_file, local_server_entry):
    state_file.add_server("http://localhost:8000", local_server_entry)
    result = state_file.get_server("http://localhost:8000")
    assert result is not None
    assert result.pid == 12345
    assert result.local_token == "test-local-token"


def test_get_nonexistent_server(state_file):
    assert state_file.get_server("http://localhost:9999") is None


def test_remove_server(state_file, local_server_entry):
    state_file.add_server("http://localhost:8000", local_server_entry)
    state_file.remove_server("http://localhost:8000")
    assert state_file.get_server("http://localhost:8000") is None


def test_remove_nonexistent_server(state_file):
    """Removing a server that doesn't exist should not raise."""
    state_file.remove_server("http://localhost:9999")


def test_add_server_preserves_tokens(state_file, local_server_entry, sample_token):
    state_file.add_token("http://localhost:8000", sample_token)
    state_file.add_server("http://localhost:8000", local_server_entry)
    assert state_file.get_token("http://localhost:8000") is not None


def test_multiple_servers(state_file, local_server_entry, remote_server_entry):
    state_file.add_server("http://localhost:8000", local_server_entry)
    state_file.add_server("https://remote.example.com", remote_server_entry)
    assert state_file.get_server("http://localhost:8000") is not None
    assert state_file.get_server("https://remote.example.com") is not None


def test_update_last_used(state_file, local_server_entry):
    state_file.add_server("http://localhost:8000", local_server_entry)
    original = state_file.get_server("http://localhost:8000")
    assert original is not None
    old_last_used = original.last_used

    state_file.update_last_used("http://localhost:8000")
    updated = state_file.get_server("http://localhost:8000")
    assert updated is not None
    assert updated.last_used > old_last_used


def test_update_last_used_nonexistent(state_file):
    """Updating last_used for nonexistent server should not raise."""
    state_file.update_last_used("http://localhost:9999")


# --- Token CRUD ---


def test_add_and_get_token(state_file, sample_token):
    state_file.add_token("https://remote.example.com", sample_token)
    result = state_file.get_token("https://remote.example.com")
    assert result is not None
    assert result.access_token == "eyJ.test.token"
    assert result.email == "user@example.com"


def test_get_nonexistent_token(state_file):
    assert state_file.get_token("https://remote.example.com") is None


def test_remove_token(state_file, sample_token):
    state_file.add_token("https://remote.example.com", sample_token)
    state_file.remove_token("https://remote.example.com")
    assert state_file.get_token("https://remote.example.com") is None


def test_remove_nonexistent_token(state_file):
    """Removing a token that doesn't exist should not raise."""
    state_file.remove_token("https://remote.example.com")


def test_remove_token_preserves_servers(state_file, local_server_entry, sample_token):
    state_file.add_server("http://localhost:8000", local_server_entry)
    state_file.add_token("http://localhost:8000", sample_token)
    state_file.remove_token("http://localhost:8000")
    assert state_file.get_server("http://localhost:8000") is not None


def test_overwrite_token(state_file, sample_token):
    state_file.add_token("https://remote.example.com", sample_token)
    new_token = TokenEntry(
        access_token="new.token",
        email="new@example.com",
        stored_at=datetime(2026, 3, 26, tzinfo=UTC),
    )
    state_file.add_token("https://remote.example.com", new_token)
    result = state_file.get_token("https://remote.example.com")
    assert result is not None
    assert result.access_token == "new.token"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_state_file.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'zndraw.state_file'`

- [ ] **Step 3: Implement StateFile**

```python
# src/zndraw/state_file.py
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
    added_at
        When this server was first registered.
    last_used
        When this server was last successfully connected to.
    pid
        Process ID (local servers only).
    version
        ZnDraw version (local servers only).
    local_token
        Per-start superuser token (local servers only).
    """

    added_at: datetime
    last_used: datetime
    pid: int | None = None
    version: str | None = None
    local_token: str | None = None


class TokenEntry(BaseModel):
    """Stored authentication token for a ZnDraw server.

    Attributes
    ----------
    access_token
        JWT access token.
    email
        Email address of the authenticated user.
    stored_at
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
    directory
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
        """Atomic write: tempfile + os.rename, mode 0600."""
        raw = data.model_dump_json(indent=2)
        fd, tmp_path = tempfile.mkstemp(dir=self.directory, suffix=".tmp")
        try:
            os.write(fd, raw.encode())
            os.fchmod(fd, stat.S_IRUSR | stat.S_IWUSR)
            os.close(fd)
            fd = -1
            os.rename(tmp_path, str(self.path))
        except BaseException:
            if fd >= 0:
                os.close(fd)
            with contextlib.suppress(OSError):
                os.unlink(tmp_path)
            raise

    # --- Server registry ---

    def add_server(self, url: str, entry: ServerEntry) -> None:
        """Register or update a server entry."""
        data = self.read()
        data.servers[url] = entry
        self._write(data)

    def remove_server(self, url: str) -> None:
        """Remove a server entry (no-op if absent)."""
        data = self.read()
        if url in data.servers:
            del data.servers[url]
            self._write(data)

    def get_server(self, url: str) -> ServerEntry | None:
        """Get a server entry by URL, or None."""
        return self.read().servers.get(url)

    def update_last_used(self, url: str) -> None:
        """Update the last_used timestamp for a server (no-op if absent)."""
        data = self.read()
        if url in data.servers:
            data.servers[url].last_used = datetime.now(UTC)
            self._write(data)

    # --- Token storage ---

    def add_token(self, url: str, entry: TokenEntry) -> None:
        """Store an authentication token for a server URL."""
        data = self.read()
        data.tokens[url] = entry
        self._write(data)

    def remove_token(self, url: str) -> None:
        """Remove a stored token (no-op if absent)."""
        data = self.read()
        if url in data.tokens:
            del data.tokens[url]
            self._write(data)

    def get_token(self, url: str) -> TokenEntry | None:
        """Get stored token for a server URL, or None."""
        return self.read().tokens.get(url)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_state_file.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/state_file.py tests/test_state_file.py
git commit -m "feat: add StateFile for unified server registry and token storage"
```

---

## Task 2: StateFile — Migration from Old Format

**Files:**
- Modify: `src/zndraw/state_file.py`
- Modify: `tests/test_state_file.py`

- [ ] **Step 1: Write failing migration tests**

Append to `tests/test_state_file.py`:

```python
# --- Migration from old format ---


def test_migrate_pid_files(tmp_path):
    """Old server-{port}.pid files are migrated to state.json."""
    # Create old-format PID file
    pid_data = {"pid": 99999, "port": 8000, "version": "0.4.0", "shutdown_token": "old"}
    (tmp_path / "server-8000.pid").write_text(json.dumps(pid_data))

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()

    entry = sf.get_server("http://localhost:8000")
    assert entry is not None
    assert entry.pid == 99999
    assert entry.version == "0.4.0"
    # shutdown_token is dropped (never validated server-side)
    assert entry.local_token is None
    # Old file is deleted
    assert not (tmp_path / "server-8000.pid").exists()


def test_migrate_tokens_json(tmp_path):
    """Old tokens.json is migrated to state.json."""
    old_tokens = {
        "https://remote.example.com": {
            "access_token": "old.jwt",
            "email": "user@example.com",
            "stored_at": "2026-03-20T10:00:00Z",
        }
    }
    (tmp_path / "tokens.json").write_text(json.dumps(old_tokens))

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()

    token = sf.get_token("https://remote.example.com")
    assert token is not None
    assert token.access_token == "old.jwt"
    assert token.email == "user@example.com"
    # Old file is deleted
    assert not (tmp_path / "tokens.json").exists()


def test_migrate_tokens_json_adds_server_entry(tmp_path):
    """Migrated token URLs also get a server entry."""
    old_tokens = {
        "https://remote.example.com": {
            "access_token": "old.jwt",
            "email": "user@example.com",
            "stored_at": "2026-03-20T10:00:00Z",
        }
    }
    (tmp_path / "tokens.json").write_text(json.dumps(old_tokens))

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()

    server = sf.get_server("https://remote.example.com")
    assert server is not None
    assert server.pid is None  # remote server


def test_migrate_both_pid_and_tokens(tmp_path):
    """Both PID files and tokens.json are migrated together."""
    pid_data = {"pid": 99999, "port": 9000, "version": "0.4.0"}
    (tmp_path / "server-9000.pid").write_text(json.dumps(pid_data))
    old_tokens = {
        "https://remote.example.com": {
            "access_token": "old.jwt",
            "email": "user@example.com",
            "stored_at": "2026-03-20T10:00:00Z",
        }
    }
    (tmp_path / "tokens.json").write_text(json.dumps(old_tokens))

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()

    assert sf.get_server("http://localhost:9000") is not None
    assert sf.get_token("https://remote.example.com") is not None
    assert not (tmp_path / "server-9000.pid").exists()
    assert not (tmp_path / "tokens.json").exists()


def test_migrate_noop_when_no_old_files(tmp_path):
    """Migration is a no-op when there are no old files."""
    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()
    # state.json should NOT be created if there's nothing to migrate
    assert not sf.path.exists()


def test_migrate_idempotent(tmp_path):
    """Running migration twice is safe."""
    pid_data = {"pid": 99999, "port": 8000, "version": "0.4.0"}
    (tmp_path / "server-8000.pid").write_text(json.dumps(pid_data))

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()
    sf.migrate_if_needed()  # second call is a no-op

    assert sf.get_server("http://localhost:8000") is not None


def test_migrate_corrupt_pid_file_skipped(tmp_path):
    """Corrupt PID files are skipped during migration."""
    (tmp_path / "server-8000.pid").write_text("not valid json")

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()

    assert sf.get_server("http://localhost:8000") is None
    # Corrupt file is still deleted (cleanup)
    assert not (tmp_path / "server-8000.pid").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_state_file.py -k migrate -v`
Expected: FAIL — `AttributeError: 'StateFile' object has no attribute 'migrate_if_needed'`

- [ ] **Step 3: Implement migration**

Add to `src/zndraw/state_file.py` — `StateFile` class:

```python
    def migrate_if_needed(self) -> None:
        """Migrate old PID files and tokens.json to state.json.

        Called on first access. Old files are deleted after migration.
        Old shutdown_token values are dropped (never validated server-side).
        """
        pid_files = sorted(self.directory.glob("server-*.pid"))
        tokens_file = self.directory / "tokens.json"

        if not pid_files and not tokens_file.exists():
            return

        data = self.read()
        now = datetime.now(UTC)

        # Migrate PID files
        for pid_file in pid_files:
            try:
                raw = json.loads(pid_file.read_text())
                port = raw["port"]
                url = f"http://localhost:{port}"
                data.servers[url] = ServerEntry(
                    added_at=now,
                    last_used=now,
                    pid=raw.get("pid"),
                    version=raw.get("version"),
                    # shutdown_token dropped intentionally
                )
            except (json.JSONDecodeError, KeyError, ValueError):
                log.warning("Skipping corrupt PID file: %s", pid_file)
            pid_file.unlink()

        # Migrate tokens.json
        if tokens_file.exists():
            try:
                raw_tokens = json.loads(tokens_file.read_text())
                for url, entry_data in raw_tokens.items():
                    data.tokens[url] = TokenEntry(**entry_data)
                    # Also register as a server if not already known
                    if url not in data.servers:
                        data.servers[url] = ServerEntry(
                            added_at=now,
                            last_used=now,
                        )
            except (json.JSONDecodeError, KeyError, ValueError):
                log.warning("Skipping corrupt tokens.json")
            tokens_file.unlink()

        if data.servers or data.tokens:
            self._write(data)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_state_file.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/state_file.py tests/test_state_file.py
git commit -m "feat: add StateFile migration from old PID files and tokens.json"
```

---

## Task 3: StateFileSource — URL and Token Resolution

**Files:**
- Create: `src/zndraw/settings_sources.py`
- Test: `tests/test_state_file_source.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_state_file_source.py
"""Tests for StateFileSource (pydantic-settings custom source)."""
from __future__ import annotations

from datetime import UTC, datetime
from unittest.mock import patch

import pytest

from zndraw.state_file import ServerEntry, StateFile, TokenEntry


@pytest.fixture
def state_file(tmp_path):
    return StateFile(directory=tmp_path)


def _local_entry(
    pid: int = 12345,
    last_used: datetime | None = None,
    local_token: str = "local-tok",
) -> ServerEntry:
    return ServerEntry(
        added_at=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
        last_used=last_used or datetime(2026, 3, 25, 14, 0, tzinfo=UTC),
        pid=pid,
        version="0.5.0",
        local_token=local_token,
    )


def _remote_entry(last_used: datetime | None = None) -> ServerEntry:
    return ServerEntry(
        added_at=datetime(2026, 3, 24, 9, 0, tzinfo=UTC),
        last_used=last_used or datetime(2026, 3, 25, 12, 0, tzinfo=UTC),
    )


def _make_source(state_file: StateFile, current_state: dict | None = None):
    """Create a StateFileSource with a test StateFile."""
    from zndraw.client.settings import ClientSettings
    from zndraw.settings_sources import StateFileSource

    source = StateFileSource(ClientSettings, state_file=state_file)
    if current_state:
        source._current_state = current_state
    return source


# --- URL resolution ---


def test_url_resolves_healthy_local_server(state_file):
    """Healthy local server is discovered as URL."""
    state_file.add_server("http://localhost:8000", _local_entry())

    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()

    assert result["url"] == "http://localhost:8000"


def test_url_skips_dead_pid(state_file):
    """Local server with dead PID is skipped."""
    state_file.add_server("http://localhost:8000", _local_entry())

    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=False),
    ):
        result = source()

    assert result.get("url") is None


def test_url_skips_unresponsive_local(state_file):
    """Local server with alive PID but failed health check is skipped."""
    state_file.add_server("http://localhost:8000", _local_entry())

    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=False),
    ):
        result = source()

    assert result.get("url") is None


def test_url_prefers_localhost_over_remote(state_file):
    """Localhost servers are preferred over remote servers."""
    state_file.add_server("https://remote.example.com", _remote_entry(
        last_used=datetime(2026, 3, 25, 15, 0, tzinfo=UTC),  # more recent
    ))
    state_file.add_server("http://localhost:8000", _local_entry(
        last_used=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),  # less recent
    ))

    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()

    assert result["url"] == "http://localhost:8000"


def test_url_falls_back_to_remote(state_file):
    """Falls back to remote when no local server is available."""
    state_file.add_server("http://localhost:8000", _local_entry())
    state_file.add_server("https://remote.example.com", _remote_entry())

    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=False),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()

    assert result["url"] == "https://remote.example.com"


def test_url_most_recent_local_wins(state_file):
    """Among multiple local servers, most recently used wins."""
    state_file.add_server("http://localhost:8000", _local_entry(
        pid=100, last_used=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
    ))
    state_file.add_server("http://localhost:9000", _local_entry(
        pid=200, last_used=datetime(2026, 3, 25, 14, 0, tzinfo=UTC),
    ))

    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()

    assert result["url"] == "http://localhost:9000"


def test_url_dead_pid_entry_removed(state_file):
    """Local entries with dead PIDs are removed from state.json."""
    state_file.add_server("http://localhost:8000", _local_entry())

    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_pid_alive", return_value=False):
        source()

    assert state_file.get_server("http://localhost:8000") is None


def test_url_empty_state_returns_none(state_file):
    """Empty state file returns no URL."""
    source = _make_source(state_file)
    result = source()
    assert result.get("url") is None


# --- Token resolution ---


def test_token_local_server_uses_local_token(state_file):
    """Local server URL resolves to its local_token."""
    state_file.add_server("http://localhost:8000", _local_entry(local_token="my-local"))

    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()

    assert result["token"] == "my-local"


def test_token_remote_uses_stored_token(state_file):
    """Remote server URL resolves to stored access_token."""
    state_file.add_server("https://remote.example.com", _remote_entry())
    state_file.add_token(
        "https://remote.example.com",
        TokenEntry(
            access_token="stored.jwt",
            email="user@example.com",
            stored_at=datetime(2026, 3, 25, tzinfo=UTC),
        ),
    )

    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_url_healthy", return_value=True):
        result = source()

    assert result["token"] == "stored.jwt"


def test_token_no_stored_returns_none(state_file):
    """Remote server with no stored token returns None for token."""
    state_file.add_server("https://remote.example.com", _remote_entry())

    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_url_healthy", return_value=True):
        result = source()

    assert result.get("token") is None


def test_token_uses_url_from_higher_source(state_file):
    """When a higher source provides URL, token is resolved for that URL."""
    state_file.add_token(
        "https://override.example.com",
        TokenEntry(
            access_token="override.jwt",
            email="user@example.com",
            stored_at=datetime(2026, 3, 25, tzinfo=UTC),
        ),
    )
    # Also register the server so it's known
    state_file.add_server("https://override.example.com", _remote_entry())

    source = _make_source(
        state_file,
        current_state={"url": "https://override.example.com"},
    )
    result = source()

    # Should use stored token for the URL provided by a higher source
    assert result.get("token") == "override.jwt"
    # Should NOT override url (higher source already set it)
    assert "url" not in result


def test_remote_unhealthy_kept_in_state(state_file):
    """Remote entries with failed health check are kept (not removed)."""
    state_file.add_server("https://remote.example.com", _remote_entry())

    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_url_healthy", return_value=False):
        source()

    # Remote entry is NOT deleted (unlike dead local PIDs)
    assert state_file.get_server("https://remote.example.com") is not None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_state_file_source.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'zndraw.settings_sources'`

- [ ] **Step 3: Implement StateFileSource**

```python
# src/zndraw/settings_sources.py
"""Custom pydantic-settings source backed by StateFile.

Resolves URL via health-check-based server discovery and token
via local_token (localhost) or stored access_token (remote).
"""
from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING, Any

import httpx
from pydantic.fields import FieldInfo
from pydantic_settings import BaseSettings, PydanticBaseSettingsSource

from zndraw.state_file import StateFile

if TYPE_CHECKING:
    from zndraw.state_file import ServerEntry

log = logging.getLogger(__name__)


def _is_pid_alive(pid: int) -> bool:
    """Check if a process is alive via os.kill(pid, 0)."""
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    return True


def _is_url_healthy(url: str, timeout: float = 2.0) -> bool:
    """Check if a server is healthy via GET /v1/health."""
    try:
        with httpx.Client(timeout=timeout) as client:
            resp = client.get(f"{url}/v1/health")
            return resp.status_code == 200
    except httpx.RequestError:
        return False


def _is_localhost(url: str) -> bool:
    """Check if a URL points to localhost."""
    return "localhost" in url or "127.0.0.1" in url


class StateFileSource(PydanticBaseSettingsSource):
    """Pydantic-settings source backed by ~/.zndraw/state.json.

    Resolves ``url`` via health-check-based server discovery
    (localhost preferred, sorted by last_used descending) and
    ``token`` via local_token (local servers) or stored access_token
    (remote servers).

    Parameters
    ----------
    settings_cls
        The settings class being configured.
    state_file
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
        self, field: FieldInfo, field_name: str
    ) -> tuple[Any, str, bool]:
        # Not used — __call__ handles everything
        return None, field_name, False

    def __call__(self) -> dict[str, Any]:
        self._state_file.migrate_if_needed()
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

        Partitions servers by localhost vs remote, sorts by last_used
        descending, tries localhost first. Dead local PIDs are cleaned up.
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

        # Sort by last_used descending (most recent first)
        local.sort(key=lambda x: x[1].last_used, reverse=True)
        remote.sort(key=lambda x: x[1].last_used, reverse=True)

        # Try localhost first
        for url, entry in local:
            if entry.pid is not None and not _is_pid_alive(entry.pid):
                log.debug("Removing dead local server: %s (PID %d)", url, entry.pid)
                self._state_file.remove_server(url)
                continue
            if _is_url_healthy(url):
                return url

        # Then try remote
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

        Local servers use local_token. Remote servers use stored access_token.
        """
        data = self._state_file.read()

        # Check if this URL is a local server with a local_token
        server = data.servers.get(url)
        if server is not None and _is_localhost(url) and server.local_token:
            return server.local_token

        # Check for stored access token
        token_entry = data.tokens.get(url)
        if token_entry is not None:
            return token_entry.access_token

        return None
```

- [ ] **Step 4: Create ClientSettings stub (needed by tests)**

```python
# src/zndraw/client/settings.py
"""Client-side settings resolved via pydantic-settings source chain.

Sources (highest to lowest priority):
    init args > env vars (ZNDRAW_*) > pyproject.toml [tool.zndraw] > StateFileSource
"""
from __future__ import annotations

from pydantic import SecretStr
from pydantic_settings import (
    BaseSettings,
    PydanticBaseSettingsSource,
    PyprojectTomlConfigSettingsSource,
    SettingsConfigDict,
)

from zndraw.settings_sources import StateFileSource


class ClientSettings(BaseSettings):
    """Connection settings for the ZnDraw client.

    Attributes
    ----------
    url
        Server URL (e.g. http://localhost:8000).
    room
        Room ID to connect to.
    user
        User email for authentication.
    password
        Password for login.
    token
        Authentication token (JWT or local_token).
    """

    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_",
        pyproject_toml_table_header=("tool", "zndraw"),
    )

    url: str | None = None
    room: str | None = None
    user: str | None = None
    password: SecretStr | None = None
    token: str | None = None

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,  # noqa: ARG003
        file_secret_settings: PydanticBaseSettingsSource,  # noqa: ARG003
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        """Return sources: init > env > pyproject.toml > state file."""
        return (
            init_settings,
            env_settings,
            PyprojectTomlConfigSettingsSource(settings_cls),
            StateFileSource(settings_cls),
        )
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_state_file_source.py -v`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/settings_sources.py src/zndraw/client/settings.py tests/test_state_file_source.py
git commit -m "feat: add StateFileSource for URL discovery and token resolution"
```

---

## Task 4: ClientSettings — Full Source Chain Tests

**Files:**
- Modify: `src/zndraw/client/settings.py` (already created in Task 3)
- Test: `tests/test_client_settings.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_client_settings.py
"""Tests for ClientSettings source chain: init > env > pyproject.toml > state file."""
from __future__ import annotations

from unittest.mock import patch

import pytest


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Remove ZNDRAW_* client env vars to isolate tests."""
    for key in ("ZNDRAW_URL", "ZNDRAW_ROOM", "ZNDRAW_USER", "ZNDRAW_PASSWORD", "ZNDRAW_TOKEN"):
        monkeypatch.delenv(key, raising=False)


@pytest.fixture
def _no_state_file():
    """Patch StateFileSource to return empty dict (no state file)."""
    with patch("zndraw.settings_sources.StateFileSource.__call__", return_value={}):
        yield


def test_init_args_highest_priority(_no_state_file, monkeypatch):
    """Init args override everything."""
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")

    from zndraw.client.settings import ClientSettings

    settings = ClientSettings(url="http://init-server:8000")
    assert settings.url == "http://init-server:8000"


def test_env_overrides_defaults(_no_state_file, monkeypatch):
    """Env vars provide values when no init args given."""
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")
    monkeypatch.setenv("ZNDRAW_ROOM", "env-room")

    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url == "http://env-server:8000"
    assert settings.room == "env-room"


def test_all_fields_default_to_none(_no_state_file):
    """All fields default to None when no source provides values."""
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url is None
    assert settings.room is None
    assert settings.user is None
    assert settings.password is None
    assert settings.token is None


def test_pyproject_toml_provides_values(_no_state_file, tmp_path, monkeypatch):
    """Values from [tool.zndraw] in pyproject.toml are used."""
    toml_content = """\
[tool.zndraw]
url = "http://toml-server:8000"
room = "toml-room"
"""
    (tmp_path / "pyproject.toml").write_text(toml_content)
    monkeypatch.chdir(tmp_path)

    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url == "http://toml-server:8000"
    assert settings.room == "toml-room"


def test_env_overrides_pyproject_toml(_no_state_file, tmp_path, monkeypatch):
    """Env vars override pyproject.toml values."""
    toml_content = """\
[tool.zndraw]
url = "http://toml-server:8000"
"""
    (tmp_path / "pyproject.toml").write_text(toml_content)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")

    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url == "http://env-server:8000"


def test_password_coerced_to_secretstr(_no_state_file, monkeypatch):
    """String password is auto-wrapped to SecretStr."""
    monkeypatch.setenv("ZNDRAW_PASSWORD", "my-secret")

    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.password is not None
    assert settings.password.get_secret_value() == "my-secret"


def test_no_namespace_overlap_with_server(_no_state_file, monkeypatch):
    """ZNDRAW_SERVER_* env vars do NOT affect ClientSettings."""
    monkeypatch.setenv("ZNDRAW_SERVER_PORT", "9999")

    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    # ClientSettings has no 'port' or 'server_port' field
    assert not hasattr(settings, "port")
    assert not hasattr(settings, "server_port")


def test_missing_pyproject_toml_silent(_no_state_file):
    """Missing pyproject.toml does not cause an error."""
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url is None
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `uv run pytest tests/test_client_settings.py -v`
Expected: All PASS (ClientSettings was implemented in Task 3)

- [ ] **Step 3: Commit**

```bash
git add tests/test_client_settings.py
git commit -m "test: add ClientSettings source chain tests"
```

---

## Task 5: Local Admin Token — Server-Side Auth Dependency

**Files:**
- Modify: `src/zndraw/dependencies.py`
- Modify: `src/zndraw/routes/admin.py`
- Test: `tests/test_local_token_auth.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_local_token_auth.py
"""Tests for local admin token authentication."""
from __future__ import annotations

import pytest
import pytest_asyncio
from httpx import ASGITransport, AsyncClient


@pytest_asyncio.fixture
async def client_with_local_token(app_with_db):
    """AsyncClient where app has a local_token set."""
    app_with_db.state.local_token = "test-local-token-secret"
    transport = ASGITransport(app=app_with_db)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        yield client


@pytest.mark.asyncio
async def test_local_token_grants_admin_access(client_with_local_token):
    """Bearer local_token should grant admin access to admin endpoints."""
    resp = await client_with_local_token.get(
        "/v1/admin/users",
        headers={"Authorization": "Bearer test-local-token-secret"},
    )
    assert resp.status_code == 200


@pytest.mark.asyncio
async def test_wrong_local_token_rejected(client_with_local_token):
    """Wrong local_token should fall through to normal auth (and fail)."""
    resp = await client_with_local_token.get(
        "/v1/admin/users",
        headers={"Authorization": "Bearer wrong-token"},
    )
    assert resp.status_code == 401


@pytest.mark.asyncio
async def test_no_local_token_configured(app_with_db):
    """When no local_token is set, only normal auth works."""
    app_with_db.state.local_token = None
    transport = ASGITransport(app=app_with_db)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        resp = await client.get(
            "/v1/admin/users",
            headers={"Authorization": "Bearer some-random-token"},
        )
        assert resp.status_code == 401
```

> **Note:** This test depends on the `app_with_db` fixture from `conftest.py`. If this fixture doesn't exist, create a minimal one or add `local_token` to the existing test app setup. Check `tests/conftest.py` for the existing fixture pattern.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_local_token_auth.py -v`
Expected: FAIL — local token is not recognized as admin auth

- [ ] **Step 3: Implement local token auth dependency**

Add to `src/zndraw/dependencies.py`:

> **Important:** `current_superuser` is a FastAPI dependency (from fastapi-users), not a plain async function — you cannot call it directly. Use `OptionalUserDep` (which returns `User | None` without raising on invalid tokens) to compose the check.

```python
async def get_local_token_or_admin(
    request: Request,
    user: OptionalUserDep,
) -> User:
    """Accept local_token as superuser OR require normal admin auth.

    When the request carries ``Authorization: Bearer <local_token>`` and
    it matches ``app.state.local_token``, a synthetic superuser identity
    is returned without touching the database.

    Otherwise, falls through to normal JWT auth via ``OptionalUserDep``.
    """
    # Path 1: Normal JWT auth succeeded → check superuser
    if user is not None:
        if not user.is_superuser:
            raise Forbidden.exception("Not a superuser")
        return user

    # Path 2: Check local admin token
    local_token: str | None = getattr(request.app.state, "local_token", None)
    if local_token is not None:
        auth_header = request.headers.get("Authorization", "")
        if auth_header == f"Bearer {local_token}":
            return User(
                email="local-admin@localhost",
                hashed_password="",
                is_superuser=True,
            )

    raise HTTPException(status_code=401, detail="Not authenticated")


LocalTokenOrAdminDep = Annotated[User, Depends(get_local_token_or_admin)]
```

Also add the missing import at the top of `dependencies.py`:

```python
from fastapi import Depends, HTTPException, Path, Request
```

- [ ] **Step 4: Update shutdown endpoint to use LocalTokenOrAdminDep**

In `src/zndraw/routes/admin.py`, change the shutdown endpoint:

```python
# Change import
from zndraw.dependencies import AdminUserDep, LocalTokenOrAdminDep, SessionDep

# Change the shutdown endpoint signature
@router.post(
    "/shutdown",
    responses=problem_responses(Forbidden),
)
async def shutdown_server(
    _admin: LocalTokenOrAdminDep,  # was AdminUserDep
) -> ShutdownResponse:
```

- [ ] **Step 5: Initialize app.state.local_token to None**

In `src/zndraw/app.py`, add after `app.state.settings_overrides = {}`:

```python
app.state.local_token = None  # CLI populates on server start
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `uv run pytest tests/test_local_token_auth.py -v`
Expected: All PASS

- [ ] **Step 7: Commit**

```bash
git add src/zndraw/dependencies.py src/zndraw/routes/admin.py src/zndraw/app.py tests/test_local_token_auth.py
git commit -m "feat: accept local_token as superuser auth for admin endpoints"
```

---

## Task 6: CLI Server Lifecycle — StateFile + Local Token

**Files:**
- Modify: `src/zndraw/cli.py`
- Modify: `src/zndraw/server_manager.py` (update `shutdown_server`)

- [ ] **Step 1: Update `resolve_server` in `src/zndraw/cli.py` to write state.json**

In `resolve_server()`, replace the PID file write block (lines ~358-370) with StateFile write:

```python
# Before:
#   shutdown_token = secrets.token_urlsafe(32)
#   ...
#   write_server_info(ServerInfo(...))
#   fastapi_app.state.shutdown_token = shutdown_token

# After:
from zndraw.state_file import ServerEntry, StateFile

local_token = secrets.token_urlsafe(32)

if detached:
    daemonize()

state = StateFile()
now = datetime.now(UTC)
state.add_server(
    f"http://{client_host}:{effective_port}",
    ServerEntry(
        added_at=now,
        last_used=now,
        pid=os.getpid(),
        version=__version__,
        local_token=local_token,
    ),
)

fastapi_app.state.local_token = local_token
fastapi_app.state.settings_overrides = settings_kwargs
```

Also add `from datetime import UTC, datetime` at the top of the file if not already imported.

Remove the old imports: `from zndraw.server_manager import ServerInfo, write_server_info, find_running_server`

Replace `find_running_server(port)` call with StateFile-based discovery. For the "reuse existing server" check at the top of `resolve_server`, use the StateFile approach.

> **Important:** Call `migrate_if_needed()` BEFORE `read()` — otherwise old PID files won't be in state.json yet. Also, the return value's third element (`effective_port`) is used for cleanup at lines 611 and 631 (`remove_server_info(effective_port)`). Change the return type to include the server URL for StateFile cleanup, and update those cleanup lines too.

Change `resolve_server` signature to return `tuple[str, uvicorn.Server | None, str]` where the third element is the server URL (for StateFile cleanup):

```python
from zndraw.state_file import StateFile
from zndraw.settings_sources import _is_pid_alive, _is_url_healthy

state = StateFile()
state.migrate_if_needed()
data = state.read()

# Check if a matching server is already running
for url, entry in data.servers.items():
    if port is not None:
        if f":{port}" not in url:
            continue
    if entry.pid is not None and not _is_pid_alive(entry.pid):
        state.remove_server(url)
        continue
    if _is_url_healthy(url):
        typer.echo(f"Found existing server (PID: {entry.pid}, URL: {url}, Version: {entry.version})")
        if entry.version and entry.version != __version__:
            typer.echo(f"Warning: Server version ({entry.version}) differs from CLI version ({__version__})")
        return url, None, url  # reuse existing — return URL for cleanup
```

At the end of `resolve_server`, where the new server entry is written, return the URL:

```python
server_url = f"http://{client_host}:{effective_port}"
# ... (StateFile write, uvicorn config, etc.)
return server_url, uvicorn.Server(config), server_url
```

Then update the cleanup code in `main()` (lines ~611 and ~631):

```python
# Before:
#   remove_server_info(effective_port)
# After:
url, server, server_url = resolve_server(...)
# ... later, on shutdown:
StateFile().remove_server(server_url)
```

The callers at lines 611 and 631 become:

```python
# Line ~611 (blocking server, no files/browser):
try:
    server.run()
except KeyboardInterrupt:
    typer.echo("\nShutting down...")
StateFile().remove_server(server_url)
typer.echo("Server stopped.")

# Line ~631 (threaded server):
try:
    thread.join()
except KeyboardInterrupt:
    typer.echo("\nShutting down...")
StateFile().remove_server(server_url)
typer.echo("Server stopped.")
```

- [ ] **Step 2: Update `shutdown_server` in `src/zndraw/server_manager.py`**

Replace the `X-Shutdown-Token` mechanism with `local_token` Bearer auth:

```python
def shutdown_server(url: str, local_token: str | None = None) -> bool:
    """Shutdown the server gracefully via API endpoint.

    Parameters
    ----------
    url
        Full server URL (e.g. http://localhost:8000).
    local_token
        Local admin token for authorization.
    """
    headers = {}
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
```

> **Note:** This changes the `shutdown_server` signature. Update all callers in `cli.py` accordingly. The old `ServerInfo`-based signature is removed.

- [ ] **Step 3: Update server stop/shutdown in `cli.py`**

Find the `--shutdown` / `handle_shutdown` function and update it to use StateFile:

```python
def handle_shutdown(port: int | None) -> None:
    """Shut down running server(s)."""
    from zndraw.settings_sources import _is_pid_alive
    from zndraw.state_file import StateFile

    state = StateFile()
    state.migrate_if_needed()
    data = state.read()
    found = False

    for url, entry in list(data.servers.items()):
        # If port specified, only target that port
        if port is not None and f":{port}" not in url:
            continue
        if entry.pid is not None and not _is_pid_alive(entry.pid):
            state.remove_server(url)
            continue
        found = True
        if shutdown_server(url, local_token=entry.local_token):
            state.remove_server(url)
            typer.echo(f"Server at {url} shut down successfully")
        else:
            typer.echo(f"Failed to shut down server at {url}", err=True)

    if not found:
        typer.echo("No running servers found.")
    raise typer.Exit()
```

- [ ] **Step 4: Run existing CLI tests**

Run: `uv run pytest tests/test_cli.py -v`
Expected: All PASS (pure function tests should be unaffected)

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/cli.py src/zndraw/server_manager.py
git commit -m "feat: CLI server lifecycle uses StateFile + local_token for auth"
```

---

## Task 7: ZnDraw.__post_init__ Refactor + auth_utils Simplification

**Files:**
- Modify: `src/zndraw/client/core.py`
- Modify: `src/zndraw/auth_utils.py`
- Modify: `tests/test_resolve_token.py`

- [ ] **Step 1: Simplify auth_utils.py**

Replace `src/zndraw/auth_utils.py` with:

```python
"""Authentication utilities for ZnDraw client.

Provides credential validation and auth helper functions.
Token resolution chain (stored token, local_token) is handled by
the pydantic-settings source chain (StateFileSource).
"""
from __future__ import annotations

import httpx
from pydantic import SecretStr


def validate_credentials(
    token: str | None,
    user: str | None,
    password: SecretStr | str | None,
) -> None:
    """Validate credential combinations (fail fast, before any network call).

    Raises
    ------
    ValueError
        On invalid credential combinations.
    """
    if token is not None and (user is not None or password is not None):
        msg = "Cannot combine --token with --user/--password"
        raise ValueError(msg)
    if user is not None and password is None:
        msg = "Missing --password (required when --user is provided)"
        raise ValueError(msg)
    if password is not None and user is None:
        msg = "Missing --user (required when --password is provided)"
        raise ValueError(msg)


def login_with_credentials(
    base_url: str,
    user: str,
    password: SecretStr | str,
) -> str:
    """Login with email+password, return access token.

    Parameters
    ----------
    base_url
        Server URL.
    user
        User email.
    password
        User password.

    Returns
    -------
    str
        JWT access token.
    """
    raw = password.get_secret_value() if isinstance(password, SecretStr) else password
    with httpx.Client(base_url=base_url, timeout=10.0) as client:
        resp = client.post(
            "/v1/auth/jwt/login",
            data={"username": user, "password": raw},
        )
        resp.raise_for_status()
        return resp.json()["access_token"]


def guest_login(base_url: str) -> str:
    """Create a guest session, return access token.

    Parameters
    ----------
    base_url
        Server URL.

    Returns
    -------
    str
        Guest JWT access token.
    """
    with httpx.Client(base_url=base_url, timeout=10.0) as client:
        resp = client.post("/v1/auth/guest")
        resp.raise_for_status()
        return resp.json()["access_token"]
```

- [ ] **Step 2: Refactor ZnDraw.__post_init__**

In `src/zndraw/client/core.py`, replace `__post_init__` (lines 129-190):

```python
    def __post_init__(self) -> None:
        """Initialize the client (REST-only, socket connects lazily)."""
        import atexit

        from zndraw.auth_utils import guest_login, login_with_credentials
        from zndraw.client.settings import ClientSettings

        # Normalize password to SecretStr
        if isinstance(self.password, str):
            self.password = SecretStr(self.password)

        # Resolve via pydantic-settings chain
        overrides = {
            k: v
            for k, v in {
                "url": self.url,
                "room": self.room,
                "user": self.user,
                "password": self.password,
                "token": self.token,
            }.items()
            if v is not None
        }
        resolved = ClientSettings(**overrides)

        # url is required — if still None after full chain, error out
        if resolved.url is None:
            raise ConnectionError(
                "No ZnDraw server found. Pass url=, set ZNDRAW_URL, "
                "add [tool.zndraw] url to pyproject.toml, or start a local server."
            )
        self.url = resolved.url
        self.room = resolved.room or str(uuid.uuid4())

        # Token resolution: settings chain > user/password login > guest
        if resolved.token is not None:
            self.token = resolved.token
        elif resolved.user and resolved.password:
            self.token = login_with_credentials(
                self.url, resolved.user, resolved.password
            )
        else:
            self.token = guest_login(self.url)

        # Create API manager
        self.api = APIManager(url=self.url, room_id=self.room, token=self.token)

        # Populate self.user for guest/stored-token sessions
        if self.user is None and self.token is not None:
            resp = self.api.http.get(
                "/v1/auth/users/me",
                headers={"Authorization": f"Bearer {self.token}"},
            )
            if resp.status_code == 200:
                self.user = resp.json().get("email")

        # Create socket manager (no connection yet — connects lazily)
        self.socket = SocketManager(zndraw=self)

        # Create job manager (zero-cost until first register())
        self._jobs = JobManager(
            api=self.api,
            tsio=self.socket.tsio,
            execute=self._execute_task if self.auto_pickup else None,
            heartbeat_interval=self.heartbeat_interval,
            polling_interval=self.polling_interval,
        )

        # Verify/create room via REST and seed frame count cache
        try:
            info = self.api.get_room_info()
            self.cached_length = info.get("frame_count", 0)
        except KeyError:
            if not self.create_if_missing:
                raise
            self.api.create_room(copy_from=self.copy_from)
            self.cached_length = 0

        atexit.register(self.disconnect)
```

Also **delete** the `_resolve_url` static method (lines ~493-517).

- [ ] **Step 2b: Update `list_rooms()` and `login()` class methods**

These class methods (lines ~519-574) also use `_resolve_url` and the old `resolve_token`. Update both:

```python
    @classmethod
    def list_rooms(
        cls,
        url: str | None = None,
        *,
        token: str | None = None,
        search: str | None = None,
    ) -> list[dict[str, Any]]:
        """List all rooms on the server.

        Parameters
        ----------
        url
            Server URL. If None, auto-discovers via state file.
        token
            JWT token. If None, uses stored token or creates guest session.
        search
            Optional search filter.
        """
        from zndraw.auth_utils import guest_login
        from zndraw.client.settings import ClientSettings

        overrides = {k: v for k, v in {"url": url, "token": token}.items() if v is not None}
        resolved = ClientSettings(**overrides)
        if resolved.url is None:
            raise ConnectionError(
                "No ZnDraw server found. Pass url= or start a local server."
            )
        resolved_token = resolved.token or guest_login(resolved.url)
        api = APIManager(url=resolved.url, room_id="", token=resolved_token)
        try:
            return api.list_rooms(search=search)
        finally:
            api.close()

    @classmethod
    def login(
        cls,
        url: str | None = None,
        username: str = "",
        password: str = "",
    ) -> str:
        """Authenticate and return a JWT token.

        Parameters
        ----------
        url
            Server URL. If None, auto-discovers via state file.
        username
            User email.
        password
            User password.
        """
        from zndraw.auth_utils import login_with_credentials
        from zndraw.client.settings import ClientSettings

        overrides = {k: v for k, v in {"url": url}.items() if v is not None}
        resolved = ClientSettings(**overrides)
        if resolved.url is None:
            raise ConnectionError(
                "No ZnDraw server found. Pass url= or start a local server."
            )
        return login_with_credentials(resolved.url, username, password)
```

- [ ] **Step 3: Update tests**

Update `tests/test_resolve_token.py` to test the simplified auth_utils:

```python
# tests/test_resolve_token.py
"""Tests for auth_utils credential validation and login helpers."""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import httpx
import pytest
from pydantic import SecretStr

from zndraw.auth_utils import guest_login, login_with_credentials, validate_credentials


# --- Validation tests (unchanged) ---


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"token": "t", "user": "u", "password": "p"}, "Cannot combine"),
        ({"token": "t", "password": "p"}, "Cannot combine"),
        ({"user": "u"}, "Missing --password"),
        ({"password": "p"}, "Missing --user"),
    ],
    ids=["token+user+password", "token+password", "user-only", "password-only"],
)
def test_invalid_credential_combinations_raise(kwargs: dict, match: str):
    with pytest.raises(ValueError, match=match):
        validate_credentials(**kwargs)


def test_valid_combinations_pass():
    validate_credentials(token="t", user=None, password=None)
    validate_credentials(token=None, user="u", password="p")
    validate_credentials(token=None, user=None, password=None)


# --- login_with_credentials ---


def test_login_with_credentials_success():
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.jwt"}
    mock_resp.raise_for_status = MagicMock()
    mock_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        result = login_with_credentials("http://localhost:8000", "user@test.com", "pass")

    assert result == "login.jwt"


def test_login_with_secretstr_password():
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.jwt"}
    mock_resp.raise_for_status = MagicMock()
    mock_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        result = login_with_credentials(
            "http://localhost:8000", "user@test.com", SecretStr("pass")
        )

    assert result == "login.jwt"
    mock_client.post.assert_called_once_with(
        "/v1/auth/jwt/login",
        data={"username": "user@test.com", "password": "pass"},
    )


# --- guest_login ---


def test_guest_login_success():
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"access_token": "guest.jwt"}
    mock_resp.raise_for_status = MagicMock()
    mock_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        result = guest_login("http://localhost:8000")

    assert result == "guest.jwt"
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_resolve_token.py tests/test_state_file.py tests/test_client_settings.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/auth_utils.py src/zndraw/client/core.py tests/test_resolve_token.py
git commit -m "refactor: ZnDraw uses ClientSettings, auth_utils simplified to login helpers"
```

---

## Task 8: CLI Connection + Auth Commands Refactor

**Files:**
- Modify: `src/zndraw/cli_agent/connection.py`
- Modify: `src/zndraw/cli_agent/auth.py`

- [ ] **Step 1: Refactor connection.py type aliases**

Remove `envvar=` from all Typer option type aliases:

```python
# Before:
UrlOpt = Annotated[str | None, typer.Option("--url", envvar="ZNDRAW_URL", help="...")]
# After:
UrlOpt = Annotated[str | None, typer.Option("--url", help="ZnDraw server URL [env: ZNDRAW_URL].")]
TokenOpt = Annotated[str | None, typer.Option("--token", help="Auth token [env: ZNDRAW_TOKEN].")]
RoomOpt = Annotated[str | None, typer.Option("--room", help="Room ID [env: ZNDRAW_ROOM].")]
UserOpt = Annotated[str | None, typer.Option("--user", help="User email [env: ZNDRAW_USER].")]
PasswordOpt = Annotated[str | None, typer.Option("--password", help="Password [env: ZNDRAW_PASSWORD].")]
```

- [ ] **Step 2: Replace resolve_url and resolve_token with ClientSettings**

Replace `resolve_url()`, `resolve_token()`, `get_connection()`, `get_zndraw()`:

```python
def get_connection(
    url: str | None,
    token: str | None,
    user: str | None = None,
    password: str | None = None,
) -> Connection:
    """Create a Connection using ClientSettings for resolution."""
    from zndraw.auth_utils import guest_login, login_with_credentials
    from zndraw.client.settings import ClientSettings

    overrides = {
        k: v
        for k, v in {"url": url, "token": token, "user": user, "password": password}.items()
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

    # Resolve token
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
    """Create a ZnDraw instance from CLI context."""
    from zndraw import ZnDraw

    return ZnDraw(
        url=url, room=room, token=token,
        user=user, password=password,
        create_if_missing=False,
    )
```

Remove the old `resolve_url()`, `resolve_token()`, `resolve_room()` functions. Also remove `get_token_store()` and the import of `find_running_server`, `TokenStore`.

Keep `resolve_room()` if it's still used elsewhere, but convert it to use ClientSettings if appropriate. Actually, `resolve_room` just validates that room is not None — keep it as-is since it's a simple check.

- [ ] **Step 3: Refactor auth.py to use StateFile**

In `src/zndraw/cli_agent/auth.py`:

```python
# Replace imports
from zndraw.state_file import StateFile, ServerEntry, TokenEntry as StateTokenEntry

# In login():
#   Replace: store = get_token_store()
#   With:    state = StateFile()
#
#   Replace: store.set(resolved_url, TokenEntry(...))
#   With:    state.add_token(resolved_url, StateTokenEntry(...))
#            state.add_server(resolved_url, ServerEntry(added_at=now, last_used=now))

# In logout():
#   Replace: store = get_token_store()
#            store.delete(resolved_url)
#   With:    state = StateFile()
#            state.remove_token(resolved_url)
#            state.remove_server(resolved_url)
```

Update `login()`:
```python
@auth_app.command("login")
def login(
    url: UrlOpt = None,
    code: bool = typer.Option(False, "--code", help="Print URL only, don't open browser"),
) -> None:
    """Login via browser approval (device-code flow)."""
    with cli_error_handler():
        from zndraw.auth_utils import guest_login
        from zndraw.client.settings import ClientSettings
        from zndraw.state_file import ServerEntry, StateFile
        from zndraw.state_file import TokenEntry as StateTokenEntry

        # Resolve URL via settings chain (but don't require auth)
        overrides = {"url": url} if url is not None else {}
        settings = ClientSettings(**overrides)
        if settings.url is None:
            from zndraw.cli_agent.connection import EXIT_CONNECTION_ERROR, die

            die(
                "No Server Found",
                "No running zndraw server found. Start one or pass --url.",
                503,
                EXIT_CONNECTION_ERROR,
            )
        resolved_url = settings.url

        state = StateFile()

        with httpx.Client(base_url=resolved_url, timeout=30.0) as client:
            # 1. Create challenge
            resp = client.post("/v1/auth/cli-login")
            resp.raise_for_status()
            challenge = resp.json()

            code_str = challenge["code"]
            secret = challenge["secret"]
            approve_url = f"{resolved_url}/auth/cli?code={code_str}"

            typer.echo(f"\n  Your code: {code_str}\n")
            if code:
                typer.echo(f"  Visit: {approve_url}")
            else:
                typer.echo(f"  Opening browser... (or visit: {approve_url})")
                webbrowser.open(approve_url)

            typer.echo("  Waiting for approval...\n")

            # 2. Poll for approval
            for _ in range(300):
                time.sleep(1)
                poll = client.get(
                    f"/v1/auth/cli-login/{code_str}",
                    params={"secret": secret},
                )

                if poll.status_code == 404:
                    typer.echo("Login rejected.", err=True)
                    raise typer.Exit(code=1)

                if poll.status_code == 410:
                    typer.echo("Login challenge expired.", err=True)
                    raise typer.Exit(code=1)

                data = poll.json()
                if data["status"] == "approved" and data["token"]:
                    token = data["token"]

                    me_resp = client.get(
                        "/v1/auth/users/me",
                        headers={"Authorization": f"Bearer {token}"},
                    )
                    email = me_resp.json().get("email", "unknown")

                    now = datetime.now(UTC)
                    state.add_token(
                        resolved_url,
                        StateTokenEntry(
                            access_token=token,
                            email=email,
                            stored_at=now,
                        ),
                    )
                    # Also register server
                    if state.get_server(resolved_url) is None:
                        state.add_server(
                            resolved_url,
                            ServerEntry(added_at=now, last_used=now),
                        )

                    typer.echo(f"Logged in as {email}")
                    typer.echo(f"Token saved to {state.path}")
                    return

            typer.echo("Login timed out.", err=True)
            raise typer.Exit(code=1)
```

Update `status()`:
```python
@auth_app.command("status")
def status(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
) -> None:
    """Show current authentication identity."""
    with cli_error_handler():
        from zndraw.auth_utils import validate_credentials
        from zndraw.client.settings import ClientSettings
        from zndraw.state_file import StateFile

        validate_credentials(token, user, password)

        overrides = {
            k: v
            for k, v in {"url": url, "token": token, "user": user, "password": password}.items()
            if v is not None
        }
        settings = ClientSettings(**overrides)
        if settings.url is None:
            from zndraw.cli_agent.connection import EXIT_CONNECTION_ERROR, die

            die("No Server Found", "Pass --url.", 503, EXIT_CONNECTION_ERROR)
        resolved_url = settings.url

        state = StateFile()

        # Determine token and source — no guest fallback
        if token is not None:
            active_token = token
            token_source = "flag"
        elif user is not None and password is not None:
            from zndraw.auth_utils import login_with_credentials

            active_token = login_with_credentials(resolved_url, user, password)
            token_source = "login"
        elif settings.token is not None:
            active_token = settings.token
            # Determine source: local_token or stored
            server = state.get_server(resolved_url)
            if server and server.local_token == active_token:
                token_source = "local_token"
            else:
                token_source = "stored"
        else:
            json_print(
                {
                    "server": resolved_url,
                    "user_id": None,
                    "email": None,
                    "is_superuser": False,
                    "token_source": "none",
                }
            )
            return

        with httpx.Client(
            base_url=resolved_url,
            headers={"Authorization": f"Bearer {active_token}"},
            timeout=10.0,
        ) as client:
            resp = client.get("/v1/auth/users/me")
            if resp.status_code != 200:
                if resp.status_code in (401, 403) and token_source == "stored":
                    state.remove_token(resolved_url)
                json_print(
                    {
                        "server": resolved_url,
                        "user_id": None,
                        "email": None,
                        "is_superuser": False,
                        "token_source": "expired"
                        if resp.status_code in (401, 403)
                        else "error",
                    }
                )
                return
            user_data = resp.json()

        json_print(
            {
                "server": resolved_url,
                "user_id": user_data.get("id"),
                "email": user_data.get("email"),
                "is_superuser": user_data.get("is_superuser", False),
                "token_source": token_source,
            }
        )
```

Update `logout()`:
```python
@auth_app.command("logout")
def logout(url: UrlOpt = None) -> None:
    """Remove stored token and server entry for the given server."""
    with cli_error_handler():
        from zndraw.client.settings import ClientSettings
        from zndraw.state_file import StateFile

        overrides = {"url": url} if url is not None else {}
        settings = ClientSettings(**overrides)
        if settings.url is None:
            from zndraw.cli_agent.connection import EXIT_CONNECTION_ERROR, die

            die("No Server Found", "Pass --url.", 503, EXIT_CONNECTION_ERROR)
        resolved_url = settings.url

        state = StateFile()
        state.remove_token(resolved_url)
        state.remove_server(resolved_url)
        typer.echo(f"Logged out from {resolved_url}")
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_cli.py tests/test_cli_auth.py -v`
Expected: PASS (update test assertions if needed for the new flow)

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/cli_agent/connection.py src/zndraw/cli_agent/auth.py
git commit -m "refactor: CLI uses ClientSettings + StateFile, remove envvar= from Typer"
```

---

## Task 9: Cleanup — Remove Old Code from server_manager.py

**Files:**
- Modify: `src/zndraw/server_manager.py`
- Delete: `tests/test_token_store.py`
- Modify: all files importing old symbols

- [ ] **Step 1: Strip server_manager.py to retained utilities**

Remove from `src/zndraw/server_manager.py`:
- `ServerInfo` dataclass
- `get_zndraw_dir()`
- `get_pid_file_path()`
- `read_server_info()`
- `write_server_info()`
- `remove_server_info()`
- `list_all_pid_files()`
- `find_running_server()`
- `get_server_status()`
- `shutdown_server()` (moved/rewritten, or keep updated version)
- `TokenEntry`
- `TokenStore`

**Keep:**
- `is_process_running(pid)` — used by StateFileSource
- `is_server_responsive(port)` — may still be useful for wait_for_server_ready
- `wait_for_server_ready(url)` — used by CLI server start
- `DEFAULT_PORT` — used by tests and CLI

The retained file should look like:

```python
"""Server management utilities for ZnDraw CLI.

Retained utilities: process checking and server readiness polling.
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
    """Check if a process with the given PID is running."""
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    return True


def is_server_responsive(port: int, timeout: float = 2.0) -> bool:
    """Check if a server is responsive on the given port."""
    try:
        with httpx.Client(timeout=timeout) as client:
            response = client.get(f"http://localhost:{port}/v1/health")
            return response.status_code == 200
    except httpx.RequestError:
        return False


def wait_for_server_ready(
    url: str, timeout: float = 30.0, poll_interval: float = 0.2
) -> bool:
    """Wait for a server to become ready by polling the health endpoint."""
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
                        "Server ready after %.2fs (%s attempts)", elapsed, attempt
                    )
                    return True
            except httpx.RequestError:
                pass

            time.sleep(current_interval)
            current_interval = min(current_interval * 1.5, max_interval)

    log.warning("Server not ready after %ss (%s attempts)", timeout, attempt)
    return False
```

- [ ] **Step 2: Fix all broken imports**

Search for all files importing from `server_manager` and update:

```bash
uv run ruff check --select I --fix .
```

Then grep for broken references:
- `from zndraw.server_manager import ServerInfo` → remove
- `from zndraw.server_manager import TokenStore, TokenEntry` → use `from zndraw.state_file import ...`
- `from zndraw.server_manager import write_server_info, read_server_info, remove_server_info` → remove
- `from zndraw.server_manager import find_running_server` → remove

- [ ] **Step 3: Delete old test file**

```bash
git rm tests/test_token_store.py
```

- [ ] **Step 4: Run full test suite**

Run: `uv run pytest tests/ -v --timeout=120`
Expected: All PASS

- [ ] **Step 5: Run type checker**

Run: `uv run pyright .`
Expected: No new errors (existing Redis/socketio false positives are acceptable)

- [ ] **Step 6: Format and lint**

```bash
uv run ruff format .
uv run ruff check --select I --fix .
```

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor: remove ServerInfo, PID files, TokenStore from server_manager"
```

---

## Verification Checklist

After all tasks complete, verify these end-to-end scenarios:

- [ ] `ZnDraw()` with zero args discovers local server and connects as admin
- [ ] `ZnDraw(url="http://localhost:8000")` connects to explicit URL
- [ ] `ZNDRAW_URL=http://... uv run python -c "from zndraw import ZnDraw; ZnDraw()"` reads env var
- [ ] `pyproject.toml` `[tool.zndraw] url = "..."` is picked up by `ZnDraw()`
- [ ] `zndraw-cli auth login --url https://remote` stores token in state.json
- [ ] `zndraw-cli auth logout --url https://remote` removes token and server entry
- [ ] `zndraw --shutdown` sends local_token Bearer auth to shutdown endpoint
- [ ] `state.json` is mode 0600
- [ ] Old `server-*.pid` and `tokens.json` files are auto-migrated
- [ ] `uv run pytest tests/ -v` — all tests pass
- [ ] `uv run pyright .` — no new type errors
