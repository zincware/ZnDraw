"""Tests for StateFileSource (pydantic-settings custom source)."""

from __future__ import annotations

import os
from datetime import UTC, datetime
from unittest.mock import patch

import pytest

from zndraw.state_file import ServerEntry, StateFile, TokenEntry

# PID that is guaranteed to not exist on any reasonable system
_DEAD_PID = 999999999

# URL guaranteed to not respond
_DEAD_URL = "http://127.0.0.1:1"


@pytest.fixture
def state_file(tmp_path):
    return StateFile(directory=tmp_path)


def _local_entry(
    pid: int = _DEAD_PID,
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


def test_url_resolves_healthy_local_server(state_file, server_factory):
    """A healthy local server is discovered and returned."""
    instance = server_factory({})
    url = instance.url
    pid = os.getpid()
    state_file.add_server(url, _local_entry(pid=pid))
    source = _make_source(state_file)
    result = source()
    assert result["url"] == url


def test_url_skips_dead_pid(state_file):
    """A local entry with a dead PID is skipped (and removed)."""
    state_file.add_server("http://localhost:8000", _local_entry(pid=_DEAD_PID))
    source = _make_source(state_file)
    result = source()
    assert result.get("url") is None


def test_url_skips_unresponsive_local(state_file):
    """A local server with alive PID but unresponsive URL is skipped."""
    pid = os.getpid()
    state_file.add_server(_DEAD_URL, _local_entry(pid=pid))
    source = _make_source(state_file)
    result = source()
    assert result.get("url") is None


def test_url_prefers_localhost_over_remote(state_file, server_factory):
    """A healthy local server is preferred over a more-recently-used remote."""
    instance = server_factory({})
    local_url = instance.url
    pid = os.getpid()

    state_file.add_server(
        _DEAD_URL,  # remote (unreachable — won't be selected)
        _remote_entry(
            last_used=datetime(2026, 3, 25, 15, 0, tzinfo=UTC),
        ),
    )
    state_file.add_server(
        local_url,
        _local_entry(
            pid=pid,
            last_used=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
        ),
    )
    source = _make_source(state_file)
    result = source()
    assert result["url"] == local_url


def test_url_falls_back_to_remote(state_file):
    """When local server has dead PID, falls back to healthy remote.

    Uses a fake non-localhost URL for the 'remote' entry since server_factory
    always binds to 127.0.0.1 (which StateFileSource classifies as local).
    The _is_url_healthy patch simulates the remote being reachable so we can
    test the local→remote fallback logic without a real remote server.
    """
    state_file.add_server("http://localhost:8000", _local_entry(pid=_DEAD_PID))
    state_file.add_server("https://remote.example.com", _remote_entry())
    source = _make_source(state_file)
    with (
        # why: pure-logic test of URL-healthy decision path without requiring remote server
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result["url"] == "https://remote.example.com"


def test_url_most_recent_local_wins(state_file, server_factory):
    """When two local servers are healthy, the most recently used one is returned."""
    instance1 = server_factory({})
    instance2 = server_factory({})
    url1 = instance1.url
    url2 = instance2.url
    pid = os.getpid()

    state_file.add_server(
        url1,
        _local_entry(
            pid=pid,
            last_used=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
        ),
    )
    state_file.add_server(
        url2,
        _local_entry(
            pid=pid,
            last_used=datetime(2026, 3, 25, 14, 0, tzinfo=UTC),
        ),
    )
    source = _make_source(state_file)
    result = source()
    # url2 has later last_used — should be preferred
    assert result["url"] == url2


def test_url_dead_pid_entry_removed(state_file):
    """A local entry with a dead PID is removed from state after discovery."""
    state_file.add_server("http://localhost:8000", _local_entry(pid=_DEAD_PID))
    source = _make_source(state_file)
    source()
    assert state_file.get_server("http://localhost:8000") is None


def test_url_empty_state_returns_none(state_file):
    """Empty state file returns no URL."""
    source = _make_source(state_file)
    result = source()
    assert result.get("url") is None


# --- Token resolution ---


def test_token_local_server_uses_access_token(state_file, server_factory):
    """Local server with access_token returns it as token."""
    instance = server_factory({})
    url = instance.url
    pid = os.getpid()
    entry = _local_entry(pid=pid, local_token="raw-admin")
    entry.access_token = "real.jwt"
    state_file.add_server(url, entry)
    source = _make_source(state_file)
    result = source()
    assert result["token"] == "real.jwt"


def test_token_local_server_no_access_token_returns_none(state_file, server_factory):
    """Local server without access_token returns no token."""
    instance = server_factory({})
    url = instance.url
    pid = os.getpid()
    state_file.add_server(url, _local_entry(pid=pid, local_token="raw-only"))
    source = _make_source(state_file)
    result = source()
    assert result.get("token") is None


def test_token_remote_uses_stored_token(state_file):
    """Remote server with stored token entry returns it."""
    state_file.add_server(_DEAD_URL, _remote_entry())
    state_file.add_token(
        _DEAD_URL,
        TokenEntry(
            access_token="stored.jwt",
            email="user@example.com",
            stored_at=datetime(2026, 3, 25, tzinfo=UTC),
        ),
    )
    source = _make_source(state_file)
    with (
        # why: pure-logic test of URL-healthy decision path without requiring remote server
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result["token"] == "stored.jwt"


def test_token_no_stored_returns_none(state_file):
    """Remote server without stored token returns no token."""
    state_file.add_server(_DEAD_URL, _remote_entry())
    source = _make_source(state_file)
    with (
        # why: pure-logic test of URL-healthy decision path without requiring remote server
        patch("zndraw.settings_sources._is_url_healthy", return_value=True),
    ):
        result = source()
    assert result.get("token") is None


def test_token_uses_url_from_higher_source(state_file, server_factory):
    """Token is resolved for URL provided by higher-priority source."""
    instance = server_factory({})
    url = instance.url
    state_file.add_token(
        url,
        TokenEntry(
            access_token="override.jwt",
            email="user@example.com",
            stored_at=datetime(2026, 3, 25, tzinfo=UTC),
        ),
    )
    state_file.add_server(url, _remote_entry())
    source = _make_source(
        state_file, current_state={"url": url}
    )
    result = source()
    assert result.get("token") == "override.jwt"
    assert "url" not in result


def test_remote_unhealthy_kept_in_state(state_file):
    """An unhealthy remote server entry is NOT removed from state."""
    state_file.add_server(_DEAD_URL, _remote_entry())
    source = _make_source(state_file)
    source()
    assert state_file.get_server(_DEAD_URL) is not None
