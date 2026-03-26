"""Regression tests: local_token must not be used as client auth token.

When the CLI starts a local server, StateFileSource resolves the auth token
from state.json. Previously it returned the raw local_token (a random string)
which the server rejected as an invalid JWT (401 Unauthorized).

The fix stores a real admin JWT in ServerEntry.access_token and
_resolve_token() returns that instead of the raw local_token.
"""

from __future__ import annotations

import os
from datetime import UTC, datetime

import pytest

from zndraw.state_file import ServerEntry, StateFile


def _local_entry(
    *,
    local_token: str = "raw-local-admin-token",  # noqa: S107
    access_token: str | None = None,
) -> ServerEntry:
    """Create a local server entry for testing."""
    return ServerEntry(
        added_at=datetime(2026, 3, 26, 10, 0, tzinfo=UTC),
        last_used=datetime(2026, 3, 26, 14, 0, tzinfo=UTC),
        pid=os.getpid(),
        version="test",
        local_token=local_token,
        access_token=access_token,
    )


def _make_source(state_file: StateFile, current_state: dict | None = None):
    """Create a StateFileSource backed by a test StateFile."""
    from zndraw.client.settings import ClientSettings
    from zndraw.settings_sources import StateFileSource

    source = StateFileSource(ClientSettings, state_file=state_file)
    if current_state:
        source._current_state = current_state
    return source


# ---------------------------------------------------------------------------
# Unit: _resolve_token must not return raw local_token for client auth
# ---------------------------------------------------------------------------


def test_resolve_token_never_returns_raw_local_token(tmp_path):
    """local_token is for admin operations only — never for client auth."""
    state = StateFile(directory=tmp_path)
    state.add_server(
        "http://localhost:8000",
        _local_entry(local_token="not-a-jwt"),
    )
    source = _make_source(state, current_state={"url": "http://localhost:8000"})
    result = source()
    token = result.get("token")
    assert token != "not-a-jwt", (
        "_resolve_token returned raw local_token — must return access_token or None"
    )


def test_resolve_token_returns_access_token_when_present(tmp_path):
    """access_token (real JWT) is preferred over local_token."""
    state = StateFile(directory=tmp_path)
    state.add_server(
        "http://localhost:8000",
        _local_entry(
            local_token="raw-admin-token",
            access_token="real.jwt.token",
        ),
    )
    source = _make_source(state, current_state={"url": "http://localhost:8000"})
    result = source()
    assert result.get("token") == "real.jwt.token"


# ---------------------------------------------------------------------------
# E2E: full CLI path — _acquire_admin_jwt → state.json → ZnDraw client
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_e2e_dev_mode_zndraw_client_connects(server, tmp_path, monkeypatch):
    """E2E dev mode: CLI acquires JWT → stores in state → client connects.

    Exercises the REAL path: _acquire_admin_jwt() → _store_jwt_in_state()
    → StateFileSource resolves access_token → ZnDraw client authenticates.
    """
    from zndraw.cli import _acquire_admin_jwt
    from zndraw.client import ZnDraw

    # Step 1: CLI acquires admin JWT (the actual function the CLI calls)
    jwt = _acquire_admin_jwt(server)
    assert jwt is not None, "_acquire_admin_jwt must return a valid JWT"

    # Step 2: Store in state.json (same as CLI does)
    state = StateFile(directory=tmp_path)
    state.add_server(
        server,
        _local_entry(local_token="raw-admin-token", access_token=jwt),
    )

    # Step 3: ZnDraw client resolves token via StateFileSource
    # Monkeypatch so StateFileSource reads our test state file
    monkeypatch.setattr(
        "zndraw.settings_sources.StateFile",
        lambda: StateFile(directory=tmp_path),
    )  # why: redirects StateFile to tmp_path for filesystem isolation; test uses real server via server_factory

    client = ZnDraw(url=server, room="test-e2e-dev")
    try:
        # If auth failed, this would raise PermissionError (the original bug)
        assert len(client) == 0
    finally:
        client.disconnect()


@pytest.mark.integration
def test_e2e_production_mode_zndraw_client_connects(server_auth, tmp_path, monkeypatch):
    """E2E production mode: CLI logs in as admin → stores JWT → client connects.

    Uses server_auth fixture which sets DEFAULT_ADMIN_EMAIL/PASSWORD.
    Exercises the login_with_credentials path of _acquire_admin_jwt().
    """
    from zndraw.cli import _acquire_admin_jwt
    from zndraw.client import ZnDraw

    # Set production mode auth env vars (same as server_auth fixture)
    monkeypatch.setenv("ZNDRAW_AUTH_DEFAULT_ADMIN_EMAIL", "admin@local.test")
    monkeypatch.setenv("ZNDRAW_AUTH_DEFAULT_ADMIN_PASSWORD", "adminpassword")

    # Step 1: CLI acquires admin JWT via login_with_credentials
    jwt = _acquire_admin_jwt(server_auth)
    assert jwt is not None, "_acquire_admin_jwt must return admin JWT"

    # Step 2: Store in state.json
    state = StateFile(directory=tmp_path)
    state.add_server(
        server_auth,
        _local_entry(local_token="raw-admin-token", access_token=jwt),
    )

    # Step 3: ZnDraw client resolves token via StateFileSource
    monkeypatch.setattr(
        "zndraw.settings_sources.StateFile",
        lambda: StateFile(directory=tmp_path),
    )  # why: redirects StateFile to tmp_path for filesystem isolation; test uses real server via server_factory

    client = ZnDraw(url=server_auth, room="test-e2e-prod")
    try:
        assert len(client) == 0
    finally:
        client.disconnect()
