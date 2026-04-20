"""Tests for application configuration."""

import os
from pathlib import Path

import pytest
from httpx import AsyncClient

from zndraw.config import Settings


class TestStorageUri:
    """Test storage URI configuration."""

    def test_default_storage_is_memory(self) -> None:
        """Default storage should be memory://."""
        settings = Settings()
        assert settings.storage == "memory://"


class TestGuestPassword:
    """Test guest password configuration."""

    def test_default_guest_password(self) -> None:
        """Default guest_password should be 'zndraw'."""
        settings = Settings()
        assert settings.guest_password.get_secret_value() == "zndraw"


class TestMediaAndServerSettings:
    """Test media path and server configuration."""

    def test_default_media_path(self) -> None:
        """Default media_path should be 'zndraw-media'."""
        settings = Settings()
        assert settings.media_path == Path("zndraw-media")

    def test_default_host_and_port(self) -> None:
        """Default host should be '0.0.0.0' and port should be 8000."""
        settings = Settings()
        assert settings.host == "0.0.0.0"
        assert settings.port == 8000


class TestSimgenEnabled:
    """Test SiMGen feature flag configuration."""

    def test_default_simgen_disabled(self) -> None:
        """SiMGen should be disabled by default."""
        settings = Settings()
        assert settings.simgen_enabled is False


class TestWorkerEnabled:
    """Test worker_enabled configuration."""

    def test_default_worker_enabled(self) -> None:
        """Worker should be enabled by default."""
        settings = Settings()
        assert settings.worker_enabled is True


class TestInternalUrl:
    """Test internal_url configuration."""

    def test_default_internal_url_is_none(self) -> None:
        """internal_url should be None by default."""
        settings = Settings()
        assert settings.internal_url is None


class TestGlobalSettingsEndpoint:
    """Test /v1/config/global-settings endpoint."""

    @pytest.mark.asyncio
    async def test_simgen_disabled_by_default(self, client: AsyncClient) -> None:
        """Endpoint returns simgen.enabled=false with default settings."""
        response = await client.get("/v1/config/global-settings")
        assert response.status_code == 200
        assert response.json() == {"simgen": {"enabled": False}}

    @pytest.mark.asyncio
    async def test_simgen_enabled_from_settings(self, client: AsyncClient) -> None:
        """Endpoint returns simgen.enabled=true when configured."""
        from zndraw.app import app

        original = app.state.settings
        app.state.settings = Settings(simgen_enabled=True)
        try:
            response = await client.get("/v1/config/global-settings")
            assert response.status_code == 200
            assert response.json() == {"simgen": {"enabled": True}}
        finally:
            app.state.settings = original


def test_settings_from_pyproject_toml(tmp_path, monkeypatch):
    """Settings should load from [tool.zndraw.server] in pyproject.toml."""
    toml_content = """\
[tool.zndraw.server]
port = 9999
host = "192.168.1.1"
storage = "/data/frames.lmdb"
"""
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(toml_content)
    monkeypatch.chdir(tmp_path)

    from zndraw.config import Settings

    settings = Settings()
    assert settings.port == 9999
    assert settings.host == "192.168.1.1"
    assert settings.storage == "/data/frames.lmdb"


def test_env_overrides_pyproject_toml(tmp_path, monkeypatch):
    """Env vars should take priority over pyproject.toml."""
    toml_content = """\
[tool.zndraw.server]
port = 9999
"""
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(toml_content)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ZNDRAW_SERVER_PORT", "7777")

    from zndraw.config import Settings

    settings = Settings()
    assert settings.port == 7777


def test_init_overrides_env_and_pyproject(tmp_path, monkeypatch):
    """Init args should take priority over env and pyproject.toml."""
    toml_content = """\
[tool.zndraw.server]
port = 9999
"""
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(toml_content)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ZNDRAW_SERVER_PORT", "7777")

    from zndraw.config import Settings

    settings = Settings(port=5555)
    assert settings.port == 5555


def test_missing_pyproject_toml_is_silent(tmp_path, monkeypatch):
    """Settings should work fine without a pyproject.toml."""
    monkeypatch.chdir(tmp_path)

    from zndraw.config import Settings

    settings = Settings()
    assert settings.port == 8000


class TestFilebrowserConfig:
    """Test filebrowser_enabled and filebrowser_path configuration."""

    def test_default_filebrowser_enabled(self) -> None:
        """Filebrowser should be enabled by default."""
        settings = Settings()
        assert settings.filebrowser_enabled is True

    def test_default_filebrowser_path_is_cwd(self) -> None:
        """Default filebrowser_path should be '.'."""
        settings = Settings()
        assert settings.filebrowser_path == "."

    def test_filebrowser_path_none_str_is_literal(self) -> None:
        """Without env_parse_none_str, lowercase 'none' is a literal path string."""
        os.environ["ZNDRAW_SERVER_FILEBROWSER_PATH"] = "none"
        try:
            settings = Settings()
            assert settings.filebrowser_path == "none"
        finally:
            os.environ.pop("ZNDRAW_SERVER_FILEBROWSER_PATH", None)


def test_guest_password_literal_none_not_coerced(monkeypatch):
    """Dropping env_parse_none_str means 'none' is a literal string
    everywhere, not an implicit None sentinel.
    """
    from pydantic import SecretStr

    from zndraw.config import Settings

    monkeypatch.setenv("ZNDRAW_SERVER_GUEST_PASSWORD", "none")
    s = Settings()
    assert isinstance(s.guest_password, SecretStr)
    assert s.guest_password.get_secret_value() == "none"


def test_task_queue_name_default():
    s = Settings()
    assert s.task_queue_name == "zndraw:tasks"


def test_result_backend_key_prefix_default():
    s = Settings()
    assert s.result_backend_key_prefix == "zndraw"


def test_provider_executor_timeout_default():
    s = Settings()
    assert s.provider_executor_timeout == 30.0
