"""Tests for unified StateFile (replaces PID files + tokens.json)."""

from __future__ import annotations

import json
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


# --- Migration from old format ---


def test_migrate_pid_files(tmp_path):
    """Old server-{port}.pid files are migrated to state.json."""
    pid_data = {"pid": 99999, "port": 8000, "version": "0.4.0", "shutdown_token": "old"}
    (tmp_path / "server-8000.pid").write_text(json.dumps(pid_data))

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()

    entry = sf.get_server("http://localhost:8000")
    assert entry is not None
    assert entry.pid == 99999
    assert entry.version == "0.4.0"
    assert entry.local_token is None
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
    assert server.pid is None


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
    assert not sf.path.exists()


def test_migrate_idempotent(tmp_path):
    """Running migration twice is safe."""
    pid_data = {"pid": 99999, "port": 8000, "version": "0.4.0"}
    (tmp_path / "server-8000.pid").write_text(json.dumps(pid_data))

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()
    sf.migrate_if_needed()

    assert sf.get_server("http://localhost:8000") is not None


def test_migrate_corrupt_pid_file_skipped(tmp_path):
    """Corrupt PID files are skipped during migration."""
    (tmp_path / "server-8000.pid").write_text("not valid json")

    sf = StateFile(directory=tmp_path)
    sf.migrate_if_needed()

    assert sf.get_server("http://localhost:8000") is None
    assert not (tmp_path / "server-8000.pid").exists()
