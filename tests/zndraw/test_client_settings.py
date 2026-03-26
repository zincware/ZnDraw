"""Tests for ClientSettings source chain: init > env > pyproject.toml > state file."""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Remove ZNDRAW_* client env vars to isolate tests."""
    for key in (
        "ZNDRAW_URL",
        "ZNDRAW_ROOM",
        "ZNDRAW_USER",
        "ZNDRAW_PASSWORD",
        "ZNDRAW_TOKEN",
    ):
        monkeypatch.delenv(key, raising=False)


@pytest.fixture(autouse=True)
def _no_state_file(tmp_path, monkeypatch):
    """Redirect StateFile to an empty tmp_path directory.

    why: redirects StateFile to empty tmp_path — StateFileSource naturally returns {}
    because no state.json exists, simulating 'no server found' without patching
    """
    from zndraw.settings_sources import StateFile

    monkeypatch.setattr(
        "zndraw.settings_sources.StateFile",
        lambda: StateFile(directory=tmp_path),
    )


def test_init_args_highest_priority(monkeypatch):
    """Init args override everything."""
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings(url="http://init-server:8000")
    assert settings.url == "http://init-server:8000"


def test_env_overrides_defaults(monkeypatch):
    """Env vars provide values when no init args given."""
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")
    monkeypatch.setenv("ZNDRAW_ROOM", "env-room")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url == "http://env-server:8000"
    assert settings.room == "env-room"


def test_all_fields_default_to_none(tmp_path, monkeypatch):
    """All fields default to None when no source provides values."""
    monkeypatch.chdir(tmp_path)
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url is None
    assert settings.room is None
    assert settings.user is None
    assert settings.password is None
    assert settings.token is None


def test_pyproject_toml_provides_values(tmp_path, monkeypatch):
    """Values from [tool.zndraw] in pyproject.toml are used."""
    toml_content = """\
[tool.zndraw]
url = "http://toml-server:8000"
room = "toml-room"
"""
    (tmp_path / "pyproject.toml").write_text(toml_content)
    monkeypatch.chdir(tmp_path)
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url == "http://toml-server:8000"
    assert settings.room == "toml-room"


def test_env_overrides_pyproject_toml(tmp_path, monkeypatch):
    """Env vars override pyproject.toml values."""
    toml_content = """\
[tool.zndraw]
url = "http://toml-server:8000"
"""
    (tmp_path / "pyproject.toml").write_text(toml_content)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url == "http://env-server:8000"


def test_password_coerced_to_secretstr(monkeypatch):
    """String password is auto-wrapped to SecretStr."""
    monkeypatch.setenv("ZNDRAW_PASSWORD", "my-secret")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.password is not None
    assert settings.password.get_secret_value() == "my-secret"


def test_no_namespace_overlap_with_server(monkeypatch):
    """ZNDRAW_SERVER_* env vars do NOT affect ClientSettings."""
    monkeypatch.setenv("ZNDRAW_SERVER_PORT", "9999")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert not hasattr(settings, "port")
    assert not hasattr(settings, "server_port")


def test_missing_pyproject_toml_silent(tmp_path, monkeypatch):
    """Missing pyproject.toml does not cause an error."""
    monkeypatch.chdir(tmp_path)
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url is None
