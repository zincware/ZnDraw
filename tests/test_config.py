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
        os.environ["ZNDRAW_STORAGE"] = "/tmp/test.lmdb"

        try:
            settings = Settings()
            assert settings.storage == "/tmp/test.lmdb"
        finally:
            os.environ.pop("ZNDRAW_STORAGE", None)

    def test_mongodb_storage_from_env(self) -> None:
        """MongoDB URI should be configurable via env var."""
        os.environ["ZNDRAW_STORAGE"] = "mongodb://mongo:27017/zndraw"

        try:
            settings = Settings()
            assert settings.storage == "mongodb://mongo:27017/zndraw"
        finally:
            os.environ.pop("ZNDRAW_STORAGE", None)


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
        os.environ["ZNDRAW_MEDIA_PATH"] = "/var/data/media"

        try:
            settings = Settings()
            assert settings.media_path == Path("/var/data/media")
        finally:
            os.environ.pop("ZNDRAW_MEDIA_PATH", None)

    def test_default_host_and_port(self) -> None:
        """Default host should be '0.0.0.0' and port should be 5000."""
        settings = Settings()
        assert settings.host == "0.0.0.0"
        assert settings.port == 5000

    def test_host_and_port_from_env(self) -> None:
        """host and port should be configurable via env vars."""
        os.environ["ZNDRAW_HOST"] = "127.0.0.1"
        os.environ["ZNDRAW_PORT"] = "8080"

        try:
            settings = Settings()
            assert settings.host == "127.0.0.1"
            assert settings.port == 8080
        finally:
            os.environ.pop("ZNDRAW_HOST", None)
            os.environ.pop("ZNDRAW_PORT", None)


class TestSimgenEnabled:
    """Test SiMGen feature flag configuration."""

    def test_default_simgen_disabled(self) -> None:
        """SiMGen should be disabled by default."""
        settings = Settings()
        assert settings.simgen_enabled is False

    def test_simgen_enabled_from_env(self) -> None:
        """SiMGen should be configurable via ZNDRAW_SIMGEN_ENABLED."""
        os.environ["ZNDRAW_SIMGEN_ENABLED"] = "true"

        try:
            settings = Settings()
            assert settings.simgen_enabled is True
        finally:
            os.environ.pop("ZNDRAW_SIMGEN_ENABLED", None)


class TestWorkerEnabled:
    """Test worker_enabled configuration."""

    def test_default_worker_enabled(self) -> None:
        """Worker should be enabled by default."""
        settings = Settings()
        assert settings.worker_enabled is True

    def test_worker_disabled_from_env(self) -> None:
        """Worker can be disabled via ZNDRAW_WORKER_ENABLED."""
        os.environ["ZNDRAW_WORKER_ENABLED"] = "false"
        try:
            settings = Settings()
            assert settings.worker_enabled is False
        finally:
            os.environ.pop("ZNDRAW_WORKER_ENABLED", None)


class TestServerUrl:
    """Test server_url configuration."""

    def test_default_server_url_is_none(self) -> None:
        """server_url should be None by default."""
        settings = Settings()
        assert settings.server_url is None

    def test_server_url_from_env(self) -> None:
        """server_url should be configurable via ZNDRAW_SERVER_URL."""
        os.environ["ZNDRAW_SERVER_URL"] = "http://nginx"
        try:
            settings = Settings()
            assert settings.server_url == "http://nginx"
        finally:
            os.environ.pop("ZNDRAW_SERVER_URL", None)


class TestSettingsFromEnv:
    """Test that Settings reads from environment variables."""

    def test_settings_returns_settings(self) -> None:
        """Settings() should return a Settings instance."""
        settings = Settings()
        assert isinstance(settings, Settings)

    def test_settings_picks_up_env_changes(self) -> None:
        """New Settings() instance should pick up env var changes."""
        os.environ["ZNDRAW_PORT"] = "9999"

        try:
            settings = Settings()
            assert settings.port == 9999
        finally:
            os.environ.pop("ZNDRAW_PORT", None)


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
