"""Tests for extension persistence across server restarts.

These tests verify that the ZnDraw client auto-reconnects and re-registers
extensions after a server restart.
"""

import time

import requests

from zndraw import ZnDraw
from zndraw.extensions import Category, Extension


class PersistenceTestExtension(Extension):
    """A test extension for verifying re-registration after server restart."""

    category = Category.MODIFIER

    def run(self, vis: ZnDraw, **kwargs):
        """Dummy run method."""
        pass


def test_extension_reregisters_after_restart(server_provider):
    """Test that client auto-reconnects and re-registers extensions after restart.

    1. Register extension (client stays connected)
    2. Restart server (Redis persists - session remains valid)
    3. Client should auto-reconnect and re-register extension
    """
    vis = ZnDraw(url=server_provider.url, room="test-room", user="worker")
    vis.register_extension(PersistenceTestExtension, public=True)
    time.sleep(0.5)

    # Verify registered
    response = requests.get(f"{server_provider.url}/api/rooms/test-room/schema/modifiers")
    assert "PersistenceTestExtension" in [ext["name"] for ext in response.json()]

    # Normal restart (Redis persists, session valid) - client should auto-reconnect
    server_provider.restart()

    # Wait for client to auto-reconnect
    for _ in range(20):  # Up to 10 seconds
        time.sleep(0.5)
        if vis.socket.connected:
            break
    else:
        raise AssertionError("Client did not reconnect within 10 seconds")

    time.sleep(1.0)

    # Extension should be available (persisted in Redis + re-registered by client)
    response = requests.get(f"{server_provider.url}/api/rooms/test-room/schema/modifiers")
    assert response.status_code == 200
    extension_names = [ext["name"] for ext in response.json()]
    assert "PersistenceTestExtension" in extension_names, (
        f"Extension not available after restart. Available: {extension_names}"
    )


def test_fresh_restart_invalidates_session(server_provider):
    """Test that fresh_restart invalidates client session (Redis is flushed).

    When Redis is flushed, the client's session is invalidated and it cannot
    reconnect with its old credentials. This is expected behavior.
    """
    vis = ZnDraw(url=server_provider.url, room="test-room", user="worker")
    vis.register_extension(PersistenceTestExtension, public=True)
    time.sleep(0.5)

    # Verify registered
    response = requests.get(f"{server_provider.url}/api/rooms/test-room/schema/modifiers")
    assert "PersistenceTestExtension" in [ext["name"] for ext in response.json()]

    # Fresh restart flushes Redis - session becomes invalid
    server_provider.fresh_restart()

    # Client cannot reconnect (session invalidated)
    time.sleep(2.0)
    assert vis.socket.connected is False, "Client should NOT reconnect after fresh restart"

    # Extension is gone (Redis was flushed)
    response = requests.get(f"{server_provider.url}/api/rooms/test-room/schema/modifiers")
    assert response.status_code == 200
    extension_names = [ext["name"] for ext in response.json()]
    assert "PersistenceTestExtension" not in extension_names, (
        "Extension should be cleared after fresh restart"
    )
