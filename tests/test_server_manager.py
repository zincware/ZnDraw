"""Unit tests for server_manager.py port-based features."""

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from zndraw.server_manager import (
    DEFAULT_PORT,
    ServerInfo,
    find_running_server,
    get_pid_file_path,
    get_server_status,
    get_zndraw_dir,
    list_all_pid_files,
    read_server_info,
    remove_server_info,
    write_server_info,
)


@pytest.fixture
def mock_zndraw_dir(tmp_path, monkeypatch):
    """Mock get_zndraw_dir to use tmp_path."""
    monkeypatch.setattr("zndraw.server_manager.get_zndraw_dir", lambda: tmp_path)
    return tmp_path


# --- PID file path tests ---


def test_get_pid_file_path_returns_port_specific_path(mock_zndraw_dir):
    """PID file path includes port number."""
    path = get_pid_file_path(5000)
    assert path == mock_zndraw_dir / "server-5000.pid"


def test_get_pid_file_path_different_ports_different_files(mock_zndraw_dir):
    """Different ports produce different PID file paths."""
    path1 = get_pid_file_path(5000)
    path2 = get_pid_file_path(5001)
    assert path1 != path2
    assert "5000" in str(path1)
    assert "5001" in str(path2)


# --- Read/Write/Remove tests ---


def test_write_and_read_server_info(mock_zndraw_dir):
    """Write server info and read it back."""
    info = ServerInfo(pid=1234, port=5000, version="1.0.0", shutdown_token="abc123")
    write_server_info(info)

    read_info = read_server_info(5000)
    assert read_info is not None
    assert read_info.pid == 1234
    assert read_info.port == 5000
    assert read_info.version == "1.0.0"
    assert read_info.shutdown_token == "abc123"


def test_read_server_info_nonexistent_returns_none(mock_zndraw_dir):
    """Reading nonexistent PID file returns None."""
    result = read_server_info(9999)
    assert result is None


def test_read_server_info_invalid_json_returns_none(mock_zndraw_dir):
    """Reading invalid JSON returns None."""
    pid_file = mock_zndraw_dir / "server-5000.pid"
    pid_file.write_text("not valid json")

    result = read_server_info(5000)
    assert result is None


def test_read_server_info_missing_fields_returns_none(mock_zndraw_dir):
    """Reading JSON with missing required fields returns None."""
    pid_file = mock_zndraw_dir / "server-5000.pid"
    pid_file.write_text('{"pid": 1234}')  # Missing port and version

    result = read_server_info(5000)
    assert result is None


def test_remove_server_info_cleans_up_file(mock_zndraw_dir):
    """Remove deletes the PID file."""
    info = ServerInfo(pid=1234, port=5000, version="1.0.0")
    write_server_info(info)

    pid_file = mock_zndraw_dir / "server-5000.pid"
    assert pid_file.exists()

    remove_server_info(5000)
    assert not pid_file.exists()


def test_remove_server_info_nonexistent_is_noop(mock_zndraw_dir):
    """Remove on nonexistent file doesn't raise."""
    # Should not raise
    remove_server_info(9999)


# --- list_all_pid_files tests ---


def test_list_all_pid_files_empty_directory(mock_zndraw_dir):
    """Empty directory returns empty list."""
    result = list_all_pid_files()
    assert result == []


def test_list_all_pid_files_multiple_servers(mock_zndraw_dir):
    """Multiple PID files are all returned."""
    # Create multiple PID files
    for port in [5000, 5001, 5002]:
        info = ServerInfo(pid=port, port=port, version="1.0.0")
        write_server_info(info)

    result = list_all_pid_files()
    assert len(result) == 3
    # Should be sorted
    assert result[0].name == "server-5000.pid"
    assert result[1].name == "server-5001.pid"
    assert result[2].name == "server-5002.pid"


def test_list_all_pid_files_ignores_other_files(mock_zndraw_dir):
    """Non-PID files are ignored."""
    # Create a PID file
    info = ServerInfo(pid=1234, port=5000, version="1.0.0")
    write_server_info(info)

    # Create other files that shouldn't be matched
    (mock_zndraw_dir / "config.json").write_text("{}")
    (mock_zndraw_dir / "server.log").write_text("log data")
    (mock_zndraw_dir / "not-server-5000.pid").write_text("{}")

    result = list_all_pid_files()
    assert len(result) == 1
    assert result[0].name == "server-5000.pid"


# --- find_running_server tests ---


def test_find_running_server_specific_port_not_running(mock_zndraw_dir):
    """Specific port with no server returns None."""
    result = find_running_server(port=5001)
    assert result is None


def test_find_running_server_specific_port_running(mock_zndraw_dir):
    """Specific port with running server returns info."""
    info = ServerInfo(pid=1234, port=5001, version="1.0.0")
    write_server_info(info)

    with (
        patch("zndraw.server_manager.is_process_running", return_value=True),
        patch("zndraw.server_manager.is_server_responsive", return_value=True),
    ):
        result = find_running_server(port=5001)

    assert result is not None
    assert result.port == 5001


def test_find_running_server_auto_discovery_default_port_first(mock_zndraw_dir):
    """Auto-discovery checks default port (5000) first."""
    # Create servers on default and non-default ports
    for port in [DEFAULT_PORT, 5001]:
        info = ServerInfo(pid=port, port=port, version="1.0.0")
        write_server_info(info)

    with (
        patch("zndraw.server_manager.is_process_running", return_value=True),
        patch("zndraw.server_manager.is_server_responsive", return_value=True),
    ):
        result = find_running_server()

    assert result is not None
    assert result.port == DEFAULT_PORT


def test_find_running_server_auto_discovery_smallest_port_fallback(mock_zndraw_dir):
    """Auto-discovery returns smallest port when default not running."""
    # Create servers on non-default ports only
    for port in [5003, 5001, 5002]:  # Out of order intentionally
        info = ServerInfo(pid=port, port=port, version="1.0.0")
        write_server_info(info)

    def is_process_running(pid):
        # All processes "running"
        return True

    def is_server_responsive(port):
        # All servers "responsive"
        return True

    with (
        patch(
            "zndraw.server_manager.is_process_running", side_effect=is_process_running
        ),
        patch(
            "zndraw.server_manager.is_server_responsive",
            side_effect=is_server_responsive,
        ),
    ):
        result = find_running_server()

    assert result is not None
    assert result.port == 5001  # Smallest non-default port


def test_find_running_server_no_servers_returns_none(mock_zndraw_dir):
    """No running servers returns None."""
    result = find_running_server()
    assert result is None


def test_find_running_server_cleans_stale_pid_files(mock_zndraw_dir):
    """Stale PID files are cleaned up during discovery."""
    # Create a PID file for a non-running process
    info = ServerInfo(pid=99999, port=5001, version="1.0.0")
    write_server_info(info)

    pid_file = mock_zndraw_dir / "server-5001.pid"
    assert pid_file.exists()

    with patch("zndraw.server_manager.is_process_running", return_value=False):
        result = find_running_server(port=5001)

    assert result is None
    # Stale file should be cleaned up
    assert not pid_file.exists()


# --- get_server_status tests ---


def test_get_server_status_no_pid_file(mock_zndraw_dir):
    """No PID file returns not running status."""
    is_running, info, message = get_server_status(5000)

    assert is_running is False
    assert info is None
    assert "No PID file" in message


def test_get_server_status_stale_pid_cleans_up(mock_zndraw_dir):
    """Stale PID file (process not running) is cleaned up."""
    info = ServerInfo(pid=99999, port=5000, version="1.0.0")
    write_server_info(info)

    pid_file = mock_zndraw_dir / "server-5000.pid"
    assert pid_file.exists()

    with patch("zndraw.server_manager.is_process_running", return_value=False):
        is_running, returned_info, message = get_server_status(5000)

    assert is_running is False
    assert returned_info is None
    assert "cleaned up" in message
    assert not pid_file.exists()


def test_get_server_status_running_server(mock_zndraw_dir):
    """Running server returns correct status."""
    info = ServerInfo(pid=1234, port=5000, version="1.0.0")
    write_server_info(info)

    with (
        patch("zndraw.server_manager.is_process_running", return_value=True),
        patch("zndraw.server_manager.is_server_responsive", return_value=True),
    ):
        is_running, returned_info, message = get_server_status(5000)

    assert is_running is True
    assert returned_info is not None
    assert returned_info.pid == 1234
    assert "running" in message.lower()


def test_get_server_status_process_running_but_not_responsive(mock_zndraw_dir):
    """Process running but server not responsive returns not running."""
    info = ServerInfo(pid=1234, port=5000, version="1.0.0")
    write_server_info(info)

    with (
        patch("zndraw.server_manager.is_process_running", return_value=True),
        patch("zndraw.server_manager.is_server_responsive", return_value=False),
    ):
        is_running, returned_info, message = get_server_status(5000)

    assert is_running is False
    assert returned_info is not None  # Info still returned for diagnostics
    assert "not responsive" in message


# --- DEFAULT_PORT constant test ---


def test_default_port_is_5000():
    """Verify default port constant."""
    assert DEFAULT_PORT == 5000
