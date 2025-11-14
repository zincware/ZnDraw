"""Test auto-discovery feature for ZnDraw client."""

from unittest.mock import patch

import pytest

from zndraw import ZnDraw
from zndraw.server_manager import ServerInfo


def test_auto_discovery_with_running_server(server):
    """Test that ZnDraw auto-discovers a running local server."""
    # The server fixture is running, so auto-discovery should find it
    # Note: This only works if the server writes its info to the expected location
    vis = ZnDraw(url=None, room="test")

    # Should have discovered the server URL
    assert vis.url == server


def test_explicit_url_bypasses_discovery(server):
    """Test that providing explicit URL bypasses auto-discovery."""
    with patch("zndraw.zndraw.get_server_status") as mock_get_status:
        # Create ZnDraw instance with explicit URL (using real server)
        vis = ZnDraw(url=server, room="test")

        # Should use the provided URL
        assert vis.url == server
        # Should NOT call get_server_status
        mock_get_status.assert_not_called()


def test_auto_discovery_no_server():
    """Test that ZnDraw raises error when no server is running."""
    with patch("zndraw.zndraw.get_server_status") as mock_get_status:
        mock_get_status.return_value = (False, None, "No server found")

        with pytest.raises(RuntimeError, match="No local ZnDraw server found"):
            ZnDraw(url=None, room="test")


def test_auto_discovery_with_stale_server():
    """Test that ZnDraw raises error when server info exists but server is not running."""
    mock_server_info = ServerInfo(pid=12345, port=5000, version="1.0.0")

    with patch("zndraw.zndraw.get_server_status") as mock_get_status:
        mock_get_status.return_value = (
            False,
            mock_server_info,
            "Process 12345 is not running",
        )

        with pytest.raises(RuntimeError, match="No local ZnDraw server found"):
            ZnDraw(url=None, room="test")
