"""Test auto-discovery feature for ZnDraw client."""

from unittest.mock import patch

import pytest

from zndraw import ZnDraw
from zndraw.server_manager import ServerInfo


def test_auto_discovery_with_running_server():
    """Test that ZnDraw auto-discovers a running local server."""
    # Mock a running server
    mock_server_info = ServerInfo(pid=12345, port=5000, version="0.6.0a4")

    with patch("zndraw.zndraw.get_server_status") as mock_get_status:
        mock_get_status.return_value = (True, mock_server_info, "Server is running")

        with patch("zndraw.api_manager.APIManager.get_version") as mock_version:
            mock_version.return_value = "0.6.0a4"

            with patch("zndraw.api_manager.APIManager.login") as mock_login:
                mock_login.return_value = {
                    "status": "success",
                    "token": "test-jwt-token",
                    "userName": "test-user",
                    "role": "guest",
                }

                with patch("zndraw.api_manager.APIManager.join_room") as mock_join:
                    mock_join.return_value = {
                        "userName": "test-user",
                        "selection": None,
                        "frame_selection": None,
                        "bookmarks": {},
                        "step": 0,
                        "geometries": None,
                        "frameCount": 0,
                    }

                    with patch("zndraw.socket_manager.SocketManager.connect"):
                        # Create ZnDraw instance with url=None
                        vis = ZnDraw(url=None, room="test")

                        # Should auto-discover and set the URL
                        assert vis.url == "http://localhost:5000"
                        mock_get_status.assert_called_once()


def test_auto_discovery_no_server():
    """Test that ZnDraw raises error when no server is running."""
    with patch("zndraw.zndraw.get_server_status") as mock_get_status:
        mock_get_status.return_value = (False, None, "No server found")

        with pytest.raises(RuntimeError, match="No local ZnDraw server found"):
            ZnDraw(url=None, room="test")


def test_explicit_url_bypasses_discovery():
    """Test that providing explicit URL bypasses auto-discovery."""
    with patch("zndraw.zndraw.get_server_status") as mock_get_status:
        with patch("zndraw.api_manager.APIManager.get_version") as mock_version:
            mock_version.return_value = "0.6.0a4"

            with patch("zndraw.api_manager.APIManager.login") as mock_login:
                mock_login.return_value = {
                    "status": "success",
                    "token": "test-jwt-token",
                    "userName": "test-user",
                    "role": "guest",
                }

                with patch("zndraw.api_manager.APIManager.join_room") as mock_join:
                    mock_join.return_value = {
                        "userName": "test-user",
                        "selection": None,
                        "frame_selection": None,
                        "bookmarks": {},
                        "step": 0,
                        "geometries": None,
                        "frameCount": 0,
                    }

                    with patch("zndraw.socket_manager.SocketManager.connect"):
                        # Create ZnDraw instance with explicit URL
                        vis = ZnDraw(url="http://example.com:8000", room="test")

                        # Should use the provided URL
                        assert vis.url == "http://example.com:8000"
                        # Should NOT call get_server_status
                        mock_get_status.assert_not_called()


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
