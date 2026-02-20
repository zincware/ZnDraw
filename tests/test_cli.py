"""Tests for ZnDraw CLI."""

import os

import pytest
from typer.testing import CliRunner

from zndraw import __version__
from zndraw.cli import (
    app,
    get_room_names,
    open_browser_to,
    path_to_room,
    sanitize_room_name,
)
from zndraw.server_manager import ServerInfo

runner = CliRunner()


# ── 1. Pure function unit tests ─────────────────────────────────────


@pytest.mark.parametrize(
    ("input_", "expected"),
    [
        ("hello world", "hello_world"),
        ("file.xyz", "file_xyz"),
        ("path/to/file", "path_to_file"),
        ("special@chars!", "special_chars_"),
        ("keep-hyphens_underscores", "keep-hyphens_underscores"),
        ("Abc123", "Abc123"),
    ],
    ids=["spaces", "dots", "slashes", "special", "hyphens", "alphanumeric"],
)
def test_sanitize_room_name(input_: str, expected: str):
    assert sanitize_room_name(input_) == expected


def test_path_to_room_unique_produces_different_names():
    a = path_to_room("file.xyz", unique=True)
    b = path_to_room("file.xyz", unique=True)
    assert a != b
    assert a.startswith("file_xyz_")
    assert b.startswith("file_xyz_")


def test_path_to_room_deterministic():
    a = path_to_room("file.xyz", unique=False)
    b = path_to_room("file.xyz", unique=False)
    assert a == b == "file_xyz"


def test_get_room_names_default_unique():
    names = get_room_names(["a.xyz", "b.xyz"], room=None, append=False)
    assert len(names) == 2
    assert names[0] != names[1]


def test_get_room_names_append_deterministic():
    names = get_room_names(["a.xyz", "b.xyz"], room=None, append=True)
    assert names == ["a_xyz", "b_xyz"]


def test_get_room_names_explicit_room():
    names = get_room_names(["a.xyz", "b.xyz", "c.xyz"], room="my-room", append=False)
    assert names == ["my-room", "my-room", "my-room"]


def test_get_room_names_empty_list():
    assert get_room_names([], room=None, append=False) == []


# ── 2. open_browser_to URL construction ─────────────────────────────


@pytest.mark.parametrize(
    ("room", "copy_from", "expected_url"),
    [
        (None, None, "http://localhost:8000"),
        ("my-room", None, "http://localhost:8000/rooms/my-room"),
        ("my-room", "@none", "http://localhost:8000/rooms/my-room?copy_from=@none"),
    ],
    ids=["root", "room", "room-copy-from-none"],
)
def test_open_browser_to_url(monkeypatch, room, copy_from, expected_url):
    opened_urls: list[str] = []
    monkeypatch.setattr("zndraw.cli.webbrowser.open", opened_urls.append)

    open_browser_to("http://localhost:8000", room, browser=True, copy_from=copy_from)
    assert opened_urls == [expected_url]


def test_open_browser_to_noop_when_disabled(monkeypatch):
    opened_urls: list[str] = []
    monkeypatch.setattr("zndraw.cli.webbrowser.open", opened_urls.append)

    open_browser_to("http://localhost:8000", "room", browser=False)
    assert opened_urls == []


# ── 3. CLI flag tests via CliRunner ──────────────────────────────────


def test_version_flag():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


@pytest.mark.parametrize(
    ("args", "error_fragment"),
    [
        (["--detached", "--status"], "--detached cannot be used with"),
        (["--detached", "--shutdown"], "--detached cannot be used with"),
        (
            ["--detached", "--connect", "http://x"],
            "--detached cannot be used with",
        ),
        (["--append", "--room", "r", "f.xyz"], "--append and --room cannot be used"),
        (["--append"], "--append and --room require file path"),
        (["--room", "r"], "--append and --room require file path"),
    ],
    ids=[
        "detached-status",
        "detached-shutdown",
        "detached-connect",
        "append-room",
        "append-no-file",
        "room-no-file",
    ],
)
def test_invalid_flag_combos(args: list[str], error_fragment: str):
    result = runner.invoke(app, args)
    assert result.exit_code == 1
    assert error_fragment in result.output


def test_file_not_found():
    result = runner.invoke(app, ["nonexistent_file.xyz", "--no-browser"])
    assert result.exit_code == 1
    assert "File not found" in result.output


# ── 4. --status / --shutdown tests ──────────────────────────────────


def test_status_no_server(monkeypatch):
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: None)
    result = runner.invoke(app, ["--status"])
    assert result.exit_code == 1
    assert "No local ZnDraw server is running" in result.output


def test_status_no_server_specific_port(monkeypatch):
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: None)
    result = runner.invoke(app, ["--status", "--port", "9999"])
    assert result.exit_code == 1
    assert "No ZnDraw server running on port 9999" in result.output


def test_status_server_running(monkeypatch):
    info = ServerInfo(pid=1234, port=8000, version="1.0.0")
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: info)
    result = runner.invoke(app, ["--status"])
    assert result.exit_code == 0
    assert "PID: 1234" in result.output
    assert "Port: 8000" in result.output


def test_shutdown_no_server(monkeypatch):
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: None)
    result = runner.invoke(app, ["--shutdown"])
    assert result.exit_code == 0
    assert "Nothing to shut down" in result.output


# ── 5. Browser-before-upload ordering ───────────────────────────────


def test_browser_before_upload_new_server(monkeypatch, tmp_path):
    """New server path: browser opens before upload_file is called."""
    dummy = tmp_path / "test.xyz"
    dummy.write_text("dummy")

    call_order: list[str] = []

    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: None)
    monkeypatch.setattr("zndraw.cli.wait_for_server_ready", lambda *_a, **_kw: True)
    monkeypatch.setattr("uvicorn.Server.run", lambda _self: None)
    monkeypatch.setattr(
        "zndraw.cli.webbrowser.open", lambda _url: call_order.append("browser")
    )
    monkeypatch.setattr(
        "zndraw.cli.upload_file", lambda *_a, **_kw: call_order.append("upload")
    )

    result = runner.invoke(app, [str(dummy)])
    assert result.exit_code == 0
    assert call_order == ["browser", "upload"]


def test_browser_before_upload_existing_server(monkeypatch, tmp_path):
    """Existing server path: browser opens before upload_file is called."""
    dummy = tmp_path / "test.xyz"
    dummy.write_text("dummy")

    call_order: list[str] = []
    info = ServerInfo(pid=1234, port=8000, version=__version__)
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: info)
    monkeypatch.setattr(
        "zndraw.cli.webbrowser.open", lambda _url: call_order.append("browser")
    )
    monkeypatch.setattr(
        "zndraw.cli.upload_file", lambda *_a, **_kw: call_order.append("upload")
    )

    result = runner.invoke(app, [str(dummy)])
    assert result.exit_code == 0
    assert call_order == ["browser", "upload"]


def test_browser_before_upload_remote(monkeypatch, tmp_path):
    """Remote server path: browser opens before upload_file is called."""
    dummy = tmp_path / "test.xyz"
    dummy.write_text("dummy")

    call_order: list[str] = []
    monkeypatch.setattr(
        "zndraw.cli.webbrowser.open", lambda _url: call_order.append("browser")
    )
    monkeypatch.setattr(
        "zndraw.cli.upload_file", lambda *_a, **_kw: call_order.append("upload")
    )

    result = runner.invoke(app, ["--connect", "http://example.com", str(dummy)])
    assert result.exit_code == 0
    assert call_order == ["browser", "upload"]


# ── 6. Existing env var propagation tests ────────────────────────────


def test_cli_sets_host_and_port_env_vars(monkeypatch):
    """CLI writes ZNDRAW_HOST and ZNDRAW_PORT env vars before starting server."""
    monkeypatch.delenv("ZNDRAW_HOST", raising=False)
    monkeypatch.delenv("ZNDRAW_PORT", raising=False)

    captured_host = None
    captured_port = None

    def fake_server_run(_self):
        nonlocal captured_host, captured_port
        captured_host = os.environ.get("ZNDRAW_HOST")
        captured_port = os.environ.get("ZNDRAW_PORT")

    monkeypatch.setattr("uvicorn.Server.run", fake_server_run)
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: None)

    result = runner.invoke(
        app, ["--port", "9999", "--host", "127.0.0.1", "--no-browser"]
    )
    assert result.exit_code == 0

    assert captured_host == "127.0.0.1"
    assert captured_port == "9999"


def test_cli_default_port_env_var(monkeypatch):
    """CLI writes default port (8000) to ZNDRAW_PORT when --port is not specified."""
    monkeypatch.delenv("ZNDRAW_HOST", raising=False)
    monkeypatch.delenv("ZNDRAW_PORT", raising=False)

    captured_port = None

    def fake_server_run(_self):
        nonlocal captured_port
        captured_port = os.environ.get("ZNDRAW_PORT")

    monkeypatch.setattr("uvicorn.Server.run", fake_server_run)
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: None)

    result = runner.invoke(app, ["--no-browser"])
    assert result.exit_code == 0

    assert captured_port == "8000"


def test_cli_reads_host_from_env(monkeypatch):
    """CLI reads ZNDRAW_HOST from env when --host is not specified."""
    monkeypatch.setenv("ZNDRAW_HOST", "192.168.1.1")
    monkeypatch.delenv("ZNDRAW_PORT", raising=False)

    captured_host = None

    def fake_server_run(_self):
        nonlocal captured_host
        captured_host = os.environ.get("ZNDRAW_HOST")

    monkeypatch.setattr("uvicorn.Server.run", fake_server_run)
    monkeypatch.setattr("zndraw.cli.find_running_server", lambda _port: None)

    result = runner.invoke(app, ["--no-browser"])
    assert result.exit_code == 0

    assert captured_host == "192.168.1.1"
