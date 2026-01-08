"""Unit tests for ServerProvider class."""

import requests

from zndraw import ZnDraw


def test_server_provider_url(server_provider):
    """ServerProvider.url returns correct URL format."""
    assert server_provider.url.startswith("http://127.0.0.1:")
    assert str(server_provider.port) in server_provider.url


def test_server_provider_is_running(server_provider):
    """ServerProvider.is_running returns True when server is running."""
    assert server_provider.is_running is True


def test_server_provider_health_endpoint(server_provider):
    """Server responds to health endpoint."""
    response = requests.get(f"{server_provider.url}/health", timeout=5)
    assert response.status_code == 200


def test_server_provider_stop(server_provider):
    """ServerProvider.stop() gracefully stops the server."""
    # Verify server is running
    assert server_provider.is_running is True

    # Stop the server
    server_provider.stop()

    # Verify server is not running
    assert server_provider.is_running is False
    assert server_provider._process is None


def test_server_provider_stop_when_not_running(server_provider):
    """ServerProvider.stop() is a no-op when server is not running."""
    server_provider.stop()  # First stop
    server_provider.stop()  # Second stop should not raise


def test_server_provider_kill(server_provider):
    """ServerProvider.kill() force kills the server."""
    assert server_provider.is_running is True

    server_provider.kill()

    assert server_provider.is_running is False
    assert server_provider._process is None


def test_server_provider_kill_when_not_running(server_provider):
    """ServerProvider.kill() is a no-op when server is not running."""
    server_provider.kill()  # First kill
    server_provider.kill()  # Second kill should not raise


def test_server_provider_restart(server_provider):
    """ServerProvider.restart() stops and starts the server."""
    original_process = server_provider._process

    url = server_provider.restart()

    assert url == server_provider.url
    assert server_provider.is_running is True
    # Process should be different after restart
    assert server_provider._process is not original_process


def test_server_provider_zndraw_connection(server_provider):
    """ZnDraw can connect to ServerProvider server."""
    vis = ZnDraw(url=server_provider.url, room="test-provider", user="TestUser")

    assert vis.url == server_provider.url
    assert vis.user == "TestUser"

    vis.disconnect()


def test_server_provider_restart_preserves_storage(server_provider):
    """Data in storage persists across restarts."""
    import ase

    # Create a room and add data
    vis = ZnDraw(url=server_provider.url, room="persistence-test", user="TestUser")
    vis.append(ase.Atoms("H2O"))
    vis.disconnect()

    # Restart the server
    server_provider.restart()

    # Reconnect and verify data persists
    vis2 = ZnDraw(url=server_provider.url, room="persistence-test", user="TestUser")
    assert len(vis2) == 1
    assert vis2[0].get_chemical_formula() == "H2O"
    vis2.disconnect()


def test_server_provider_flush_redis(server_provider):
    """ServerProvider.flush_redis() clears Redis data."""
    import redis

    # Add some data to Redis
    client = redis.Redis.from_url(server_provider.redis_url, decode_responses=True)
    client.set("test-key", "test-value")
    assert client.get("test-key") == "test-value"

    # Flush Redis
    server_provider.flush_redis()

    # Verify data is cleared
    assert client.get("test-key") is None


def test_server_provider_fresh_restart(server_provider):
    """ServerProvider.fresh_restart() clears all data and restarts."""
    import ase

    # Create a room and add data
    vis = ZnDraw(url=server_provider.url, room="fresh-test", user="TestUser")
    vis.append(ase.Atoms("H2O"))
    vis.disconnect()

    # Fresh restart (clears Redis)
    server_provider.fresh_restart()

    # Server should be running
    assert server_provider.is_running is True

    # Room data should be cleared (room won't exist or will be empty)
    vis2 = ZnDraw(url=server_provider.url, room="fresh-test", user="TestUser")
    # After fresh restart, room is recreated empty
    assert len(vis2) == 0
    vis2.disconnect()
