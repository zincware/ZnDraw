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

    def test_storage_from_env(self) -> None:
        """Storage URI should be configurable via env var."""
        os.environ["ZNDRAW_SERVER_STORAGE"] = "/tmp/test.lmdb"

        try:
            settings = Settings()
            assert settings.storage == "/tmp/test.lmdb"
        finally:
            os.environ.pop("ZNDRAW_SERVER_STORAGE", None)

    def test_mongodb_storage_from_env(self) -> None:
        """MongoDB URI should be configurable via env var."""
        os.environ["ZNDRAW_SERVER_STORAGE"] = "mongodb://mongo:27017/zndraw"

        try:
            settings = Settings()
            assert settings.storage == "mongodb://mongo:27017/zndraw"
        finally:
            os.environ.pop("ZNDRAW_SERVER_STORAGE", None)


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

    def test_media_path_from_env(self) -> None:
        """media_path should be configurable via env var."""
        os.environ["ZNDRAW_SERVER_MEDIA_PATH"] = "/var/data/media"

        try:
            settings = Settings()
            assert settings.media_path == Path("/var/data/media")
        finally:
            os.environ.pop("ZNDRAW_SERVER_MEDIA_PATH", None)

    def test_default_host_and_port(self) -> None:
        """Default host should be '0.0.0.0' and port should be 8000."""
        settings = Settings()
        assert settings.host == "0.0.0.0"
        assert settings.port == 8000

    def test_host_and_port_from_env(self) -> None:
        """host and port should be configurable via env vars."""
        os.environ["ZNDRAW_SERVER_HOST"] = "127.0.0.1"
        os.environ["ZNDRAW_SERVER_PORT"] = "8080"

        try:
            settings = Settings()
            assert settings.host == "127.0.0.1"
            assert settings.port == 8080
        finally:
            os.environ.pop("ZNDRAW_SERVER_HOST", None)
            os.environ.pop("ZNDRAW_SERVER_PORT", None)


class TestSimgenEnabled:
    """Test SiMGen feature flag configuration."""

    def test_default_simgen_disabled(self) -> None:
        """SiMGen should be disabled by default."""
        settings = Settings()
        assert settings.simgen_enabled is False

    def test_simgen_enabled_from_env(self) -> None:
        """SiMGen should be configurable via ZNDRAW_SERVER_SIMGEN_ENABLED."""
        os.environ["ZNDRAW_SERVER_SIMGEN_ENABLED"] = "true"

        try:
            settings = Settings()
            assert settings.simgen_enabled is True
        finally:
            os.environ.pop("ZNDRAW_SERVER_SIMGEN_ENABLED", None)


class TestWorkerEnabled:
    """Test worker_enabled configuration."""

    def test_default_worker_enabled(self) -> None:
        """Worker should be enabled by default."""
        settings = Settings()
        assert settings.worker_enabled is True

    def test_worker_disabled_from_env(self) -> None:
        """Worker can be disabled via ZNDRAW_SERVER_WORKER_ENABLED."""
        os.environ["ZNDRAW_SERVER_WORKER_ENABLED"] = "false"
        try:
            settings = Settings()
            assert settings.worker_enabled is False
        finally:
            os.environ.pop("ZNDRAW_SERVER_WORKER_ENABLED", None)


class TestInternalUrl:
    """Test internal_url configuration."""

    def test_default_internal_url_is_none(self) -> None:
        """internal_url should be None by default."""
        settings = Settings()
        assert settings.internal_url is None

    def test_internal_url_from_env(self) -> None:
        """internal_url should be configurable via ZNDRAW_SERVER_INTERNAL_URL."""
        os.environ["ZNDRAW_SERVER_INTERNAL_URL"] = "http://nginx"
        try:
            settings = Settings()
            assert settings.internal_url == "http://nginx"
        finally:
            os.environ.pop("ZNDRAW_SERVER_INTERNAL_URL", None)


class TestSettingsFromEnv:
    """Test that Settings reads from environment variables."""

    def test_settings_returns_settings(self) -> None:
        """Settings() should return a Settings instance."""
        settings = Settings()
        assert isinstance(settings, Settings)

    def test_settings_picks_up_env_changes(self) -> None:
        """New Settings() instance should pick up env var changes."""
        os.environ["ZNDRAW_SERVER_PORT"] = "9999"

        try:
            settings = Settings()
            assert settings.port == 9999
        finally:
            os.environ.pop("ZNDRAW_SERVER_PORT", None)


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


class TestFilebrowserPath:
    """Test filebrowser_path configuration."""

    def test_default_filebrowser_path_is_cwd(self) -> None:
        """Default filebrowser_path should be '.'."""
        settings = Settings()
        assert settings.filebrowser_path == "."

    def test_filebrowser_path_from_env(self) -> None:
        """filebrowser_path should be configurable via ZNDRAW_SERVER_FILEBROWSER_PATH."""
        os.environ["ZNDRAW_SERVER_FILEBROWSER_PATH"] = "/data"
        try:
            settings = Settings()
            assert settings.filebrowser_path == "/data"
        finally:
            os.environ.pop("ZNDRAW_SERVER_FILEBROWSER_PATH", None)

    def test_filebrowser_path_none_sentinel(self) -> None:
        """Sentinel 'none' (case-insensitive) disables the default provider."""
        os.environ["ZNDRAW_SERVER_FILEBROWSER_PATH"] = "NONE"
        try:
            settings = Settings()
            assert settings.filebrowser_path == "NONE"
            assert settings.filebrowser_path.lower() == "none"
        finally:
            os.environ.pop("ZNDRAW_SERVER_FILEBROWSER_PATH", None)
