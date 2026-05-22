"""Tests for ZnDraw CLI."""

from datetime import UTC, datetime
from types import SimpleNamespace
from uuid import UUID

import httpx
import pytest
import typer
from httpx import MockTransport, Response
from typer.testing import CliRunner

from zndraw import __version__
from zndraw.cli import (
    _acquire_token,
    _resolve_owner_id,
    _validate_room_arg,
    app,
    get_room_names,
    open_browser_to,
    path_to_room,
    sanitize_room_name,
)
from zndraw.config import Settings
from zndraw.state_file import ServerEntry, StateFile

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
    monkeypatch.setattr(
        "zndraw.cli.webbrowser.open", opened_urls.append
    )  # why: prevents browser window during test

    open_browser_to("http://localhost:8000", room, browser=True, copy_from=copy_from)
    assert opened_urls == [expected_url]


def test_open_browser_to_noop_when_disabled(monkeypatch):
    opened_urls: list[str] = []
    monkeypatch.setattr(
        "zndraw.cli.webbrowser.open", opened_urls.append
    )  # why: prevents browser window during test

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


def _empty_state(monkeypatch, tmp_path):
    """Point StateFile at an empty tmp dir, disable health checks, stub auth."""
    monkeypatch.setattr(
        "zndraw.cli.StateFile", lambda: StateFile(directory=tmp_path)
    )  # why: isolates state to tmp_path for filesystem isolation
    monkeypatch.setattr(
        "zndraw.cli._is_url_healthy", lambda _url: False
    )  # why: simulates no existing server for StateFile logic
    monkeypatch.setattr(
        "zndraw.cli._acquire_token", lambda _url: "tok"
    )  # why: avoid real HTTP auth calls in unit tests
    monkeypatch.setattr(
        "zndraw.cli._resolve_owner_id",
        lambda _url, _tok: UUID("12345678-1234-5678-1234-567812345678"),
    )  # why: avoid real HTTP auth calls in unit tests


def test_status_no_server(monkeypatch, tmp_path):
    _empty_state(monkeypatch, tmp_path)
    result = runner.invoke(app, ["--status"])
    assert result.exit_code == 1
    assert "No local ZnDraw server is running" in result.output


def test_status_no_server_specific_port(monkeypatch, tmp_path):
    _empty_state(monkeypatch, tmp_path)
    result = runner.invoke(app, ["--status", "--port", "9999"])
    assert result.exit_code == 1
    assert "No ZnDraw server running on port 9999" in result.output


def test_status_server_running(monkeypatch, tmp_path):
    state = StateFile(directory=tmp_path)
    now = datetime.now(UTC)
    state.add_server(
        "http://localhost:8000",
        ServerEntry(added_at=now, last_used=now, pid=1234, version="1.0.0"),
    )
    monkeypatch.setattr(
        "zndraw.cli.StateFile", lambda: state
    )  # why: isolates state to tmp_path for filesystem isolation
    monkeypatch.setattr(
        "zndraw.cli._is_url_healthy", lambda _url: True
    )  # why: simulates no existing server for StateFile logic

    result = runner.invoke(app, ["--status"])
    assert result.exit_code == 0
    assert "PID: 1234" in result.output


def test_shutdown_no_server(monkeypatch, tmp_path):
    _empty_state(monkeypatch, tmp_path)
    result = runner.invoke(app, ["--shutdown"])
    assert result.exit_code == 0
    assert "Nothing to shut down" in result.output


# ── 5. Browser-before-upload ordering ───────────────────────────────


def test_browser_before_upload_new_server(monkeypatch, tmp_path):
    """New server path: browser opens before upload_file is called."""
    dummy = tmp_path / "test.xyz"
    dummy.write_text("dummy")

    call_order: list[str] = []
    state_dir = tmp_path / "state"
    state_dir.mkdir()

    monkeypatch.setattr(
        "zndraw.cli.StateFile", lambda: StateFile(directory=state_dir)
    )  # why: isolates state to tmp_path for filesystem isolation
    monkeypatch.setattr(
        "zndraw.cli._is_url_healthy", lambda _url: False
    )  # why: simulates no existing server for StateFile logic
    monkeypatch.setattr(
        "zndraw.cli.wait_for_server_ready", lambda *_a, **_kw: True
    )  # why: skips server polling in unit test
    monkeypatch.setattr(
        "uvicorn.Server.run", lambda _self: None
    )  # why: prevents real server startup in unit test
    monkeypatch.setattr(
        "zndraw.cli.webbrowser.open", lambda _url: call_order.append("browser")
    )  # why: prevents browser window during test
    monkeypatch.setattr(
        "zndraw.cli.upload_file", lambda *_a, **_kw: call_order.append("upload")
    )  # why: tracks call order (browser-before-upload orchestration test)
    monkeypatch.setattr(
        "zndraw.cli._acquire_token", lambda _url: "tok"
    )  # why: avoid real /v1/auth/* HTTP from the new-server post-startup hook
    monkeypatch.setattr(
        "zndraw.cli._resolve_owner_id",
        lambda _url, _tok: UUID("12345678-1234-5678-1234-567812345678"),
    )  # why: same reason — Task 4 will need this too

    result = runner.invoke(app, [str(dummy)])
    assert result.exit_code == 0
    assert call_order == ["browser", "upload"]


def test_browser_before_upload_existing_server(monkeypatch, tmp_path):
    """Existing server path: browser opens before upload_file is called."""
    dummy = tmp_path / "test.xyz"
    dummy.write_text("dummy")

    call_order: list[str] = []
    state_dir = tmp_path / "state"
    state = StateFile(directory=state_dir)
    now = datetime.now(UTC)
    state.add_server(
        "http://localhost:8000",
        ServerEntry(added_at=now, last_used=now, pid=1234, version=__version__),
    )
    monkeypatch.setattr(
        "zndraw.cli.StateFile", lambda: state
    )  # why: isolates state to tmp_path for filesystem isolation
    monkeypatch.setattr(
        "zndraw.cli._is_url_healthy", lambda _url: True
    )  # why: simulates no existing server for StateFile logic
    monkeypatch.setattr(
        "zndraw.cli.webbrowser.open", lambda _url: call_order.append("browser")
    )  # why: prevents browser window during test
    monkeypatch.setattr(
        "zndraw.cli.upload_file", lambda *_a, **_kw: call_order.append("upload")
    )  # why: tracks call order (browser-before-upload orchestration test)
    monkeypatch.setattr(
        "zndraw.cli._acquire_token", lambda _url: "tok"
    )  # why: avoid real HTTP auth calls in unit tests
    monkeypatch.setattr(
        "zndraw.cli._resolve_owner_id",
        lambda _url, _tok: UUID("12345678-1234-5678-1234-567812345678"),
    )  # why: avoid real HTTP auth calls in unit tests

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
    )  # why: prevents browser window during test
    monkeypatch.setattr(
        "zndraw.cli.upload_file", lambda *_a, **_kw: call_order.append("upload")
    )  # why: tracks call order (browser-before-upload orchestration test)
    monkeypatch.setattr(
        "zndraw.cli._acquire_token", lambda _url: "tok"
    )  # why: avoid real HTTP auth calls in unit tests
    monkeypatch.setattr(
        "zndraw.cli._resolve_owner_id",
        lambda _url, _tok: UUID("12345678-1234-5678-1234-567812345678"),
    )  # why: avoid real HTTP auth calls in unit tests

    result = runner.invoke(app, ["--connect", "http://example.com", str(dummy)])
    assert result.exit_code == 0
    assert call_order == ["browser", "upload"]


# ── 6. Settings propagation tests ────────────────────────────────────


@pytest.fixture
def capture_settings(monkeypatch, tmp_path):
    """Spy on Settings instantiation and stub out server startup."""
    captured: list[Settings] = []
    original_init = Settings.__init__

    def spy_init(self, **kwargs):
        original_init(self, **kwargs)
        captured.append(self)

    monkeypatch.setattr(
        Settings, "__init__", spy_init
    )  # why: spy on Settings instantiation to verify config propagation
    monkeypatch.setattr(
        "uvicorn.Server.run", lambda _self: None
    )  # why: prevents real server startup in unit test
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    monkeypatch.setattr(
        "zndraw.cli.StateFile", lambda: StateFile(directory=state_dir)
    )  # why: isolates state to tmp_path for filesystem isolation
    monkeypatch.setattr(
        "zndraw.cli._is_url_healthy", lambda _url: False
    )  # why: simulates no existing server for StateFile logic
    monkeypatch.setattr(
        "zndraw.cli.wait_for_server_ready",
        lambda _url, timeout=30.0: True,  # why: skips server polling in unit test
    )
    monkeypatch.setattr(
        "zndraw.cli._acquire_token", lambda _url: "test-token"
    )  # why: unit test of Settings propagation, not auth flow
    monkeypatch.setattr(
        "zndraw.cli._resolve_owner_id",
        lambda _url, _tok: UUID("12345678-1234-5678-1234-567812345678"),
    )  # why: avoid live HTTP from Settings-propagation tests
    return captured


def test_cli_passes_host_and_port_to_settings(capture_settings):
    """CLI --host and --port are forwarded to Settings via app.state."""
    result = runner.invoke(
        app, ["--port", "9999", "--host", "127.0.0.1", "--no-browser"]
    )
    assert result.exit_code == 0
    assert any(s.host == "127.0.0.1" and s.port == 9999 for s in capture_settings)


def test_cli_default_port_from_settings(capture_settings):
    """When --port is not specified, Settings default (8000) is used."""
    result = runner.invoke(app, ["--no-browser"])
    assert result.exit_code == 0
    assert any(s.port == 8000 for s in capture_settings)


def test_cli_reads_host_from_env(monkeypatch, capture_settings):
    """Settings reads ZNDRAW_SERVER_HOST from env when --host is not specified."""
    monkeypatch.setenv("ZNDRAW_SERVER_HOST", "192.168.1.1")
    monkeypatch.delenv("ZNDRAW_SERVER_PORT", raising=False)

    result = runner.invoke(app, ["--no-browser"])
    assert result.exit_code == 0
    assert any(s.host == "192.168.1.1" for s in capture_settings)


# ── 7. _resolve_owner_id ────────────────────────────────────────────


def _mock_httpx_with(monkeypatch, handler):
    """Force every httpx.Client built inside cli.py to use a MockTransport."""
    real_client = httpx.Client

    def fake_client(*args, **kwargs):
        kwargs["transport"] = MockTransport(handler)
        return real_client(*args, **kwargs)

    monkeypatch.setattr("zndraw.cli.httpx.Client", fake_client)


def test_resolve_owner_id_returns_uuid(monkeypatch):
    captured: dict[str, str] = {}

    def handler(request):
        captured["url"] = str(request.url)
        captured["auth"] = request.headers.get("Authorization", "")
        return Response(200, json={"id": "12345678-1234-5678-1234-567812345678"})

    _mock_httpx_with(monkeypatch, handler)

    result = _resolve_owner_id("http://test", "tok")
    assert result == UUID("12345678-1234-5678-1234-567812345678")
    assert captured["url"].endswith("/v1/auth/users/me")
    assert captured["auth"] == "Bearer tok"


def test_resolve_owner_id_raises_on_http_error(monkeypatch):
    _mock_httpx_with(monkeypatch, lambda _req: Response(403))
    with pytest.raises(httpx.HTTPStatusError):
        _resolve_owner_id("http://test", "tok")


# ── 8. _validate_room_arg ───────────────────────────────────────────


def test_validate_room_arg_accepts_composed():
    owner = UUID("12345678-1234-5678-1234-567812345678")
    # composed-form: <uuid>/<name>; should not raise
    _validate_room_arg(f"{owner}/proj", owner)


def test_validate_room_arg_rejects_bare():
    owner = UUID("12345678-1234-5678-1234-567812345678")
    with pytest.raises(typer.BadParameter) as exc:
        _validate_room_arg("foo", owner)
    msg = str(exc.value)
    assert "must be '<owner_uuid>/<name>'" in msg
    assert "foo" in msg
    assert str(owner) in msg  # surfaces caller's UUID for copy-paste


def test_validate_room_arg_rejects_bad_uuid():
    owner = UUID("12345678-1234-5678-1234-567812345678")
    with pytest.raises(typer.BadParameter) as exc:
        _validate_room_arg("not-a-uuid/foo", owner)
    assert "not a valid UUID" in str(exc.value)


def test_validate_room_arg_rejects_bad_name_chars():
    owner = UUID("12345678-1234-5678-1234-567812345678")
    with pytest.raises(typer.BadParameter) as exc:
        _validate_room_arg(f"{owner}/bad name", owner)
    assert "invalid characters" in str(exc.value)


# ── 9. _acquire_token chain ─────────────────────────────────────────


def test_acquire_token_uses_stored_token(monkeypatch):
    monkeypatch.setattr(
        "zndraw.cli.ClientSettings", lambda **_: SimpleNamespace(token="stored")
    )
    monkeypatch.setattr(
        "zndraw.cli.resolve_or_refresh_token",
        lambda url, tok: f"resolved:{tok}",
    )

    def _fail(*_a, **_kw):
        pytest.fail("should not be called when a stored token is present")

    monkeypatch.setattr("zndraw.cli.guest_login", _fail)
    monkeypatch.setattr("zndraw.cli.login_with_credentials", _fail)

    assert _acquire_token("http://x") == "resolved:stored"


def test_acquire_token_dev_mode_uses_guest_login(monkeypatch):
    monkeypatch.setattr(
        "zndraw.cli.ClientSettings", lambda **_: SimpleNamespace(token=None)
    )
    monkeypatch.setattr(
        "zndraw.cli.AuthSettings",
        lambda: SimpleNamespace(
            is_dev_mode=True,
            default_admin_email=None,
            default_admin_password=None,
        ),
    )
    monkeypatch.setattr("zndraw.cli.guest_login", lambda _url: "guest-tok")

    def _fail(*_a, **_kw):
        pytest.fail("login_with_credentials must not be called in dev mode")

    monkeypatch.setattr("zndraw.cli.login_with_credentials", _fail)

    assert _acquire_token("http://x") == "guest-tok"


def test_acquire_token_admin_creds(monkeypatch):
    monkeypatch.setattr(
        "zndraw.cli.ClientSettings", lambda **_: SimpleNamespace(token=None)
    )
    monkeypatch.setattr(
        "zndraw.cli.AuthSettings",
        lambda: SimpleNamespace(
            is_dev_mode=False,
            default_admin_email="admin@example.com",
            default_admin_password="pw",
        ),
    )

    monkeypatch.setattr(
        "zndraw.cli.login_with_credentials",
        lambda url, user, pw: f"admin:{user}:{pw}",
    )

    def _fail(_url):
        pytest.fail("guest_login must not be called when admin creds exist")

    monkeypatch.setattr("zndraw.cli.guest_login", _fail)

    assert _acquire_token("http://x") == "admin:admin@example.com:pw"


def test_acquire_token_falls_back_to_guest(monkeypatch):
    """No stored token, not dev mode, missing admin password → guest_login."""
    monkeypatch.setattr(
        "zndraw.cli.ClientSettings", lambda **_: SimpleNamespace(token=None)
    )
    monkeypatch.setattr(
        "zndraw.cli.AuthSettings",
        lambda: SimpleNamespace(
            is_dev_mode=False,
            default_admin_email="admin@example.com",
            default_admin_password=None,
        ),
    )
    monkeypatch.setattr("zndraw.cli.guest_login", lambda _url: "guest-fallback")

    assert _acquire_token("http://x") == "guest-fallback"


# ── 10. main() composes rooms via owner_id ──────────────────────────


FIXED_OWNER = UUID("12345678-1234-5678-1234-567812345678")


def _stub_token_chain(monkeypatch):
    """Stub server discovery + auth helpers used by main()."""
    monkeypatch.setattr("zndraw.cli._acquire_token", lambda _url: "tok")
    monkeypatch.setattr("zndraw.cli._resolve_owner_id", lambda _url, _tok: FIXED_OWNER)


def test_main_composes_room_from_path(monkeypatch, tmp_path):
    """Without --room, paths produce '<owner>/<sanitized_path>_<rand>' names."""
    dummy = tmp_path / "trj.xyz"
    dummy.write_text("dummy")

    _empty_state(monkeypatch, tmp_path)
    monkeypatch.setattr("zndraw.cli.wait_for_server_ready", lambda *_a, **_kw: True)
    monkeypatch.setattr("uvicorn.Server.run", lambda _self: None)
    _stub_token_chain(monkeypatch)

    uploads: list[tuple[str, str, str]] = []
    monkeypatch.setattr(
        "zndraw.cli.upload_file",
        lambda path, url, room_id, *_a, **_kw: uploads.append((path, url, room_id)),
    )
    monkeypatch.setattr("zndraw.cli.webbrowser.open", lambda _url: None)

    result = runner.invoke(app, [str(dummy)])
    assert result.exit_code == 0, result.output
    assert len(uploads) == 1
    _, _, room_id = uploads[0]
    assert room_id.startswith(f"{FIXED_OWNER}/trj_xyz_")


def test_main_uses_explicit_composed_room_without_reprefixing(monkeypatch, tmp_path):
    """--room <other_uuid>/proj must NOT be re-prefixed with the caller's UUID."""
    dummy = tmp_path / "trj.xyz"
    dummy.write_text("dummy")

    _empty_state(monkeypatch, tmp_path)
    monkeypatch.setattr("zndraw.cli.wait_for_server_ready", lambda *_a, **_kw: True)
    monkeypatch.setattr("uvicorn.Server.run", lambda _self: None)
    _stub_token_chain(monkeypatch)

    uploads: list[str] = []
    monkeypatch.setattr(
        "zndraw.cli.upload_file",
        lambda path, url, room_id, *_a, **_kw: uploads.append(room_id),
    )
    monkeypatch.setattr("zndraw.cli.webbrowser.open", lambda _url: None)

    explicit = "aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa/proj"
    result = runner.invoke(app, ["--room", explicit, str(dummy)])
    assert result.exit_code == 0, result.output
    assert uploads == [explicit]


def test_main_errors_on_bare_room(monkeypatch, tmp_path):
    """zndraw --room foo <file> → typer.BadParameter; message shows caller UUID."""
    dummy = tmp_path / "trj.xyz"
    dummy.write_text("dummy")

    _empty_state(monkeypatch, tmp_path)
    monkeypatch.setattr("zndraw.cli.wait_for_server_ready", lambda *_a, **_kw: True)
    monkeypatch.setattr("uvicorn.Server.run", lambda _self: None)
    _stub_token_chain(monkeypatch)
    monkeypatch.setattr("zndraw.cli.upload_file", lambda *_a, **_kw: None)
    monkeypatch.setattr("zndraw.cli.webbrowser.open", lambda _url: None)

    result = runner.invoke(app, ["--room", "foo", str(dummy)], env={"NO_COLOR": "1"})
    assert result.exit_code != 0
    assert "must be '<owner_uuid>/<name>'" in result.output
    assert str(FIXED_OWNER) in result.output


def test_main_synthesizes_workspace_when_no_paths_no_room(monkeypatch, tmp_path):
    """zndraw (no args, --no-browser) → server starts; no upload."""
    _empty_state(monkeypatch, tmp_path)
    monkeypatch.setattr("zndraw.cli.wait_for_server_ready", lambda *_a, **_kw: True)
    monkeypatch.setattr("uvicorn.Server.run", lambda _self: None)
    _stub_token_chain(monkeypatch)
    upload_calls: list[object] = []
    monkeypatch.setattr(
        "zndraw.cli.upload_file",
        lambda *_a, **_kw: upload_calls.append(_a),
    )
    monkeypatch.setattr("zndraw.cli.webbrowser.open", lambda _url: None)

    result = runner.invoke(app, ["--no-browser"])
    assert result.exit_code == 0, result.output
    assert upload_calls == []
