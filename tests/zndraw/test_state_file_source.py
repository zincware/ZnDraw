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
    local_token: str = "local-tok",  # noqa: S107
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
    state_file.add_server("http://localhost:8000", _local_entry())
    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result["url"] == "http://localhost:8000"


def test_url_skips_dead_pid(state_file):
    state_file.add_server("http://localhost:8000", _local_entry())
    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_pid_alive", return_value=False):
        result = source()
    assert result.get("url") is None


def test_url_skips_unresponsive_local(state_file):
    state_file.add_server("http://localhost:8000", _local_entry())
    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=False),
    ):
        result = source()
    assert result.get("url") is None


def test_url_prefers_localhost_over_remote(state_file):
    state_file.add_server(
        "https://remote.example.com",
        _remote_entry(
            last_used=datetime(2026, 3, 25, 15, 0, tzinfo=UTC),
        ),
    )
    state_file.add_server(
        "http://localhost:8000",
        _local_entry(
            last_used=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
        ),
    )
    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result["url"] == "http://localhost:8000"


def test_url_falls_back_to_remote(state_file):
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
    state_file.add_server(
        "http://localhost:8000",
        _local_entry(
            pid=100,
            last_used=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
        ),
    )
    state_file.add_server(
        "http://localhost:9000",
        _local_entry(
            pid=200,
            last_used=datetime(2026, 3, 25, 14, 0, tzinfo=UTC),
        ),
    )
    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result["url"] == "http://localhost:9000"


def test_url_dead_pid_entry_removed(state_file):
    state_file.add_server("http://localhost:8000", _local_entry())
    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_pid_alive", return_value=False):
        source()
    assert state_file.get_server("http://localhost:8000") is None


def test_url_empty_state_returns_none(state_file):
    source = _make_source(state_file)
    result = source()
    assert result.get("url") is None


# --- Token resolution ---


def test_token_local_server_uses_access_token(state_file):
    entry = _local_entry(local_token="raw-admin")
    entry.access_token = "real.jwt"
    state_file.add_server("http://localhost:8000", entry)
    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result["token"] == "real.jwt"


def test_token_local_server_no_access_token_returns_none(state_file):
    state_file.add_server(
        "http://localhost:8000", _local_entry(local_token="raw-only")
    )
    source = _make_source(state_file)
    with (
        patch("zndraw.settings_sources._is_pid_alive", return_value=True),
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result.get("token") is None


def test_token_remote_uses_stored_token(state_file):
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
    state_file.add_server("https://remote.example.com", _remote_entry())
    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_url_healthy", return_value=True):
        result = source()
    assert result.get("token") is None


def test_token_uses_url_from_higher_source(state_file):
    state_file.add_token(
        "https://override.example.com",
        TokenEntry(
            access_token="override.jwt",
            email="user@example.com",
            stored_at=datetime(2026, 3, 25, tzinfo=UTC),
        ),
    )
    state_file.add_server("https://override.example.com", _remote_entry())
    source = _make_source(
        state_file, current_state={"url": "https://override.example.com"}
    )
    result = source()
    assert result.get("token") == "override.jwt"
    assert "url" not in result


def test_remote_unhealthy_kept_in_state(state_file):
    state_file.add_server("https://remote.example.com", _remote_entry())
    source = _make_source(state_file)
    with patch("zndraw.settings_sources._is_url_healthy", return_value=False):
        source()
    assert state_file.get_server("https://remote.example.com") is not None
