"""E2E test: StateFileSource discovers a running server via state.json.

Uses a real uvicorn server via server_factory to verify that StateFileSource
can locate the server from a state.json file backed by a real running process.
"""

from __future__ import annotations

import os
from datetime import UTC, datetime
from pathlib import Path

from zndraw.state_file import ServerEntry, StateFile


def _make_source(state_file: StateFile):
    """Create a StateFileSource pointing at the given StateFile."""
    from zndraw.client.settings import ClientSettings
    from zndraw.settings_sources import StateFileSource

    return StateFileSource(ClientSettings, state_file=state_file)


# =============================================================================
# Server discovery via state.json
# =============================================================================


def test_statefile_discovers_running_server(server_factory, tmp_path: Path):
    """StateFileSource returns the URL of a registered running server.

    Steps:
    1. Start a real server via server_factory
    2. Create a StateFile pointing at tmp_path
    3. Register the server with the current PID (alive process)
    4. Call StateFileSource — it should return the server URL
    """
    instance = server_factory({})
    url = instance.url

    state_file = StateFile(directory=tmp_path)
    entry = ServerEntry(
        added_at=datetime.now(UTC),
        last_used=datetime.now(UTC),
        pid=os.getpid(),
    )
    state_file.add_server(url, entry)

    source = _make_source(state_file)
    result = source()

    assert result.get("url") == url


def test_statefile_discovers_server_with_access_token(server_factory, tmp_path: Path):
    """StateFileSource returns token stored in server entry for local server.

    A local server entry with access_token set should have that token
    returned as the resolved token.
    """
    instance = server_factory({})
    url = instance.url

    state_file = StateFile(directory=tmp_path)
    entry = ServerEntry(
        added_at=datetime.now(UTC),
        last_used=datetime.now(UTC),
        pid=os.getpid(),
        access_token="real.jwt.token",
    )
    state_file.add_server(url, entry)

    source = _make_source(state_file)
    result = source()

    assert result.get("url") == url
    assert result.get("token") == "real.jwt.token"


def test_statefile_prefers_most_recent_server(server_factory, tmp_path: Path):
    """When multiple servers are registered, the most recently used one wins."""
    instance1 = server_factory({})
    instance2 = server_factory({})
    url1 = instance1.url
    url2 = instance2.url
    pid = os.getpid()

    state_file = StateFile(directory=tmp_path)
    # Register url1 with an older last_used timestamp
    state_file.add_server(
        url1,
        ServerEntry(
            added_at=datetime.now(UTC),
            last_used=datetime(2026, 3, 25, 10, 0, tzinfo=UTC),
            pid=pid,
        ),
    )
    # Register url2 with a more recent last_used timestamp
    state_file.add_server(
        url2,
        ServerEntry(
            added_at=datetime.now(UTC),
            last_used=datetime(2026, 3, 25, 14, 0, tzinfo=UTC),
            pid=pid,
        ),
    )

    source = _make_source(state_file)
    result = source()

    # url2 has later last_used — should be preferred
    assert result.get("url") == url2


def test_statefile_no_url_when_dead_process(tmp_path: Path):
    """StateFileSource returns no URL if the registered process is dead."""
    state_file = StateFile(directory=tmp_path)
    # Use a PID that is guaranteed to not exist
    dead_pid = 999_999_999
    state_file.add_server(
        "http://localhost:9999",
        ServerEntry(
            added_at=datetime.now(UTC),
            last_used=datetime.now(UTC),
            pid=dead_pid,
        ),
    )

    source = _make_source(state_file)
    result = source()

    assert result.get("url") is None


def test_statefile_server_removed_after_dead_pid(tmp_path: Path):
    """StateFileSource removes dead-PID entries from state.json after discovery."""
    state_file = StateFile(directory=tmp_path)
    dead_pid = 999_999_999
    state_file.add_server(
        "http://localhost:9999",
        ServerEntry(
            added_at=datetime.now(UTC),
            last_used=datetime.now(UTC),
            pid=dead_pid,
        ),
    )

    source = _make_source(state_file)
    source()

    # Dead entry should have been cleaned up
    assert state_file.get_server("http://localhost:9999") is None


def test_statefile_url_healthy_check(server_factory, tmp_path: Path):
    """StateFileSource uses health check to verify server is reachable."""
    import httpx

    instance = server_factory({})
    url = instance.url

    # Verify the server is actually reachable before registering
    resp = httpx.get(f"{url}/v1/health", timeout=5.0)
    assert resp.status_code == 200, "Server must be healthy for this E2E test"

    state_file = StateFile(directory=tmp_path)
    state_file.add_server(
        url,
        ServerEntry(
            added_at=datetime.now(UTC),
            last_used=datetime.now(UTC),
            pid=os.getpid(),
        ),
    )

    source = _make_source(state_file)
    result = source()

    assert result.get("url") == url
