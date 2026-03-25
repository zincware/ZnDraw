"""Tests for ClientSettings source chain: init > env > pyproject.toml > state file."""

from __future__ import annotations

from unittest.mock import patch

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


@pytest.fixture
def _no_state_file():
    """Patch StateFileSource to return empty dict (no state file)."""
    with patch("zndraw.settings_sources.StateFileSource.__call__", return_value={}):
        yield


def test_init_args_highest_priority(_no_state_file, monkeypatch):
    """Init args override everything."""
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings(url="http://init-server:8000")
    assert settings.url == "http://init-server:8000"


def test_env_overrides_defaults(_no_state_file, monkeypatch):
    """Env vars provide values when no init args given."""
    monkeypatch.setenv("ZNDRAW_URL", "http://env-server:8000")
    monkeypatch.setenv("ZNDRAW_ROOM", "env-room")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url == "http://env-server:8000"
    assert settings.room == "env-room"


def test_all_fields_default_to_none(_no_state_file):
    """All fields default to None when no source provides values."""
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url is None
    assert settings.room is None
    assert settings.user is None
    assert settings.password is None
    assert settings.token is None


def test_pyproject_toml_provides_values(_no_state_file, tmp_path, monkeypatch):
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


def test_env_overrides_pyproject_toml(_no_state_file, tmp_path, monkeypatch):
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


def test_password_coerced_to_secretstr(_no_state_file, monkeypatch):
    """String password is auto-wrapped to SecretStr."""
    monkeypatch.setenv("ZNDRAW_PASSWORD", "my-secret")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.password is not None
    assert settings.password.get_secret_value() == "my-secret"


def test_no_namespace_overlap_with_server(_no_state_file, monkeypatch):
    """ZNDRAW_SERVER_* env vars do NOT affect ClientSettings."""
    monkeypatch.setenv("ZNDRAW_SERVER_PORT", "9999")
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert not hasattr(settings, "port")
    assert not hasattr(settings, "server_port")


def test_missing_pyproject_toml_silent(_no_state_file):
    """Missing pyproject.toml does not cause an error."""
    from zndraw.client.settings import ClientSettings

    settings = ClientSettings()
    assert settings.url is None
