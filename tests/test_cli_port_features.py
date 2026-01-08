"""Tests for CLI port-related features (--status, --shutdown, --port)."""

from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from zndraw.cli import app
from zndraw.server_manager import DEFAULT_PORT, ServerInfo

runner = CliRunner()


# =============================================================================
# --status flag tests
# =============================================================================


def test_status_no_server_running():
    """--status with no server returns exit code 1."""
    with patch("zndraw.cli.find_running_server", return_value=None):
        result = runner.invoke(app, ["--status"])

    assert result.exit_code == 1
    assert "No local ZnDraw server is running" in result.output


def test_status_server_running():
    """--status with running server returns exit code 0 and displays info."""
    mock_info = ServerInfo(pid=1234, port=5000, version="1.0.0")

    with patch("zndraw.cli.find_running_server", return_value=mock_info):
        result = runner.invoke(app, ["--status"])

    assert result.exit_code == 0
    assert "Server running" in result.output
    assert "PID: 1234" in result.output
    assert "Port: 5000" in result.output
    assert "Version: 1.0.0" in result.output
    assert "http://localhost:5000" in result.output


def test_status_with_port_no_server():
    """--status --port 5001 with no server on that port returns exit code 1."""
    with patch("zndraw.cli.find_running_server", return_value=None) as mock_find:
        result = runner.invoke(app, ["--status", "--port", "5001"])

    assert result.exit_code == 1
    assert "No ZnDraw server running on port 5001" in result.output
    # Verify find_running_server was called with the specific port
    mock_find.assert_called_once_with(5001)


def test_status_with_port_server_running():
    """--status --port 5001 with server on that port returns info."""
    mock_info = ServerInfo(pid=5678, port=5001, version="2.0.0")

    with patch("zndraw.cli.find_running_server", return_value=mock_info) as mock_find:
        result = runner.invoke(app, ["--status", "--port", "5001"])

    assert result.exit_code == 0
    assert "Port: 5001" in result.output
    mock_find.assert_called_once_with(5001)


def test_status_auto_discovers_without_port():
    """--status without --port calls find_running_server with None (auto-discover)."""
    mock_info = ServerInfo(pid=1234, port=DEFAULT_PORT, version="1.0.0")

    with patch("zndraw.cli.find_running_server", return_value=mock_info) as mock_find:
        result = runner.invoke(app, ["--status"])

    assert result.exit_code == 0
    # Should call with None to trigger auto-discovery
    mock_find.assert_called_once_with(None)


# =============================================================================
# --shutdown flag tests
# =============================================================================


def test_shutdown_no_server_running():
    """--shutdown with no server exits gracefully."""
    with patch("zndraw.cli.find_running_server", return_value=None):
        result = runner.invoke(app, ["--shutdown"])

    assert result.exit_code == 0
    assert "No running server found" in result.output


def test_shutdown_server_running_success():
    """--shutdown with running server shuts it down."""
    mock_info = ServerInfo(pid=1234, port=5000, version="1.0.0")

    with (
        patch("zndraw.cli.find_running_server", return_value=mock_info),
        patch("zndraw.cli.shutdown_server", return_value=True) as mock_shutdown,
    ):
        result = runner.invoke(app, ["--shutdown"])

    assert result.exit_code == 0
    assert "Shutting down server" in result.output
    assert "PID: 1234" in result.output
    assert "shut down successfully" in result.output
    mock_shutdown.assert_called_once_with(mock_info)


def test_shutdown_server_running_failure():
    """--shutdown that fails returns exit code 1."""
    mock_info = ServerInfo(pid=1234, port=5000, version="1.0.0")

    with (
        patch("zndraw.cli.find_running_server", return_value=mock_info),
        patch("zndraw.cli.shutdown_server", return_value=False),
    ):
        result = runner.invoke(app, ["--shutdown"])

    assert result.exit_code == 1
    assert "Failed to shut down server" in result.output


def test_shutdown_with_port_no_server():
    """--shutdown --port 5001 with no server on that port exits gracefully."""
    with patch("zndraw.cli.find_running_server", return_value=None) as mock_find:
        result = runner.invoke(app, ["--shutdown", "--port", "5001"])

    assert result.exit_code == 0
    assert "No server running on port 5001" in result.output
    mock_find.assert_called_once_with(5001)


def test_shutdown_with_port_server_running():
    """--shutdown --port 5001 shuts down server on that specific port."""
    mock_info = ServerInfo(pid=5678, port=5001, version="2.0.0")

    with (
        patch("zndraw.cli.find_running_server", return_value=mock_info) as mock_find,
        patch("zndraw.cli.shutdown_server", return_value=True),
    ):
        result = runner.invoke(app, ["--shutdown", "--port", "5001"])

    assert result.exit_code == 0
    assert "Port: 5001" in result.output
    mock_find.assert_called_once_with(5001)


def test_shutdown_auto_discovers_without_port():
    """--shutdown without --port calls find_running_server with None."""
    mock_info = ServerInfo(pid=1234, port=DEFAULT_PORT, version="1.0.0")

    with (
        patch("zndraw.cli.find_running_server", return_value=mock_info) as mock_find,
        patch("zndraw.cli.shutdown_server", return_value=True),
    ):
        result = runner.invoke(app, ["--shutdown"])

    assert result.exit_code == 0
    mock_find.assert_called_once_with(None)


# =============================================================================
# --port help text test
# =============================================================================


def test_port_help_text():
    """--port help text explains auto-discovery behavior."""
    # Set TERM=dumb to disable Rich's fancy rendering (boxes, colors, cursor control)
    # and COLUMNS for consistent width. This ensures consistent output in CI
    # environments where Rich may render incorrectly due to no TTY.
    result = runner.invoke(app, ["--help"], env={"TERM": "dumb", "COLUMNS": "200"})

    assert result.exit_code == 0
    # Verify help text mentions --port option
    assert "--port" in result.output


# =============================================================================
# Integration tests with real server (subprocess-based)
# =============================================================================


def test_status_with_real_server(server):
    """--status finds a real running server."""
    # Extract port from server URL
    port = int(server.split(":")[-1])

    result = runner.invoke(app, ["--status", "--port", str(port)])

    assert result.exit_code == 0
    assert "Server running" in result.output
    assert f"Port: {port}" in result.output


def test_status_wrong_port_with_real_server(server, get_free_port):
    """--status --port with wrong port shows no server."""
    # Get a different port that's not running
    other_port = get_free_port()

    result = runner.invoke(app, ["--status", "--port", str(other_port)])

    assert result.exit_code == 1
    assert f"No ZnDraw server running on port {other_port}" in result.output
