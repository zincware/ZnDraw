"""Test auto-discovery feature for ZnDraw client."""

from unittest.mock import patch

import pytest

from zndraw import ZnDraw


def test_auto_discovery_with_running_server(server):
    """Test that ZnDraw auto-discovers a running local server."""
    # The server fixture is running, so auto-discovery should find it
    # Note: This only works if the server writes its info to the expected location
    vis = ZnDraw(url=None, room="test")

    # Should have discovered the server URL
    assert vis.url == server


def test_explicit_url_bypasses_discovery(server):
    """Test that providing explicit URL bypasses auto-discovery."""
    with patch("zndraw.zndraw.find_running_server") as mock_find:
        # Create ZnDraw instance with explicit URL (using real server)
        vis = ZnDraw(url=server, room="test")

        # Should use the provided URL
        assert vis.url == server
        # Should NOT call find_running_server
        mock_find.assert_not_called()


def test_auto_discovery_no_server():
    """Test that ZnDraw raises error when no server is running."""
    with patch("zndraw.zndraw.find_running_server") as mock_find:
        mock_find.return_value = None

        with pytest.raises(RuntimeError, match="No local ZnDraw server found"):
            ZnDraw(url=None, room="test")


def test_auto_discovery_with_stale_server():
    """Test that ZnDraw raises error when server info exists but server is not running."""
    # When find_running_server returns None, it means no running server was found
    # (stale server info is automatically cleaned up by find_running_server)
    with patch("zndraw.zndraw.find_running_server") as mock_find:
        mock_find.return_value = None

        with pytest.raises(RuntimeError, match="No local ZnDraw server found"):
            ZnDraw(url=None, room="test")
