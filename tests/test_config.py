"""Tests for storage discriminated union configuration."""

import os
from pathlib import Path

import pytest
from httpx import AsyncClient

from zndraw.config import (
    LMDBStorage,
    MemoryStorage,
    MongoDBStorage,
    Settings,
)


class TestStorageDiscriminatedUnion:
    """Test storage configuration discriminated union."""

    def test_default_storage_is_memory(self) -> None:
        """Default storage should be MemoryStorage."""
        settings = Settings()
        assert isinstance(settings.storage, MemoryStorage)
        assert settings.storage.type == "memory"

    def test_memory_storage_explicit(self) -> None:
        """MemoryStorage can be explicitly created."""
        storage = MemoryStorage()
        assert storage.type == "memory"

    def test_lmdb_storage_default_map_size(self) -> None:
        """LMDBStorage should have default map_size of 1GB."""
        storage = LMDBStorage(path=Path("/tmp/lmdb"))
        assert storage.type == "lmdb"
        assert storage.path == Path("/tmp/lmdb")
        assert storage.map_size == 1024 * 1024 * 1024  # 1GB

    def test_lmdb_storage_custom_map_size(self) -> None:
        """LMDBStorage should accept custom map_size."""
        storage = LMDBStorage(path=Path("/tmp/lmdb"), map_size=512 * 1024 * 1024)
        assert storage.map_size == 512 * 1024 * 1024

    def test_mongodb_storage_default_database(self) -> None:
        """MongoDBStorage should have default database 'zndraw'."""
        storage = MongoDBStorage(url="mongodb://localhost:27017")
        assert storage.type == "mongodb"
        assert storage.url == "mongodb://localhost:27017"
        assert storage.database == "zndraw"

    def test_mongodb_storage_custom_database(self) -> None:
        """MongoDBStorage should accept custom database name."""
        storage = MongoDBStorage(url="mongodb://localhost:27017", database="custom_db")
        assert storage.database == "custom_db"


class TestStorageFromEnvVars:
    """Test storage configuration from environment variables."""

    def test_lmdb_storage_from_env(self) -> None:
        """LMDB storage should be configurable via env vars."""
        os.environ["ZNDRAW_STORAGE__TYPE"] = "lmdb"
        os.environ["ZNDRAW_STORAGE__PATH"] = "/tmp/test-lmdb"
        os.environ["ZNDRAW_STORAGE__MAP_SIZE"] = "536870912"  # 512MB

        try:
            settings = Settings()
            assert isinstance(settings.storage, LMDBStorage)
            assert settings.storage.type == "lmdb"
            assert settings.storage.path == Path("/tmp/test-lmdb")
            assert settings.storage.map_size == 536870912
        finally:
            os.environ.pop("ZNDRAW_STORAGE__TYPE", None)
            os.environ.pop("ZNDRAW_STORAGE__PATH", None)
            os.environ.pop("ZNDRAW_STORAGE__MAP_SIZE", None)

    def test_mongodb_storage_from_env(self) -> None:
        """MongoDB storage should be configurable via env vars."""
        os.environ["ZNDRAW_STORAGE__TYPE"] = "mongodb"
        os.environ["ZNDRAW_STORAGE__URL"] = "mongodb://mongo.example.com:27017"
        os.environ["ZNDRAW_STORAGE__DATABASE"] = "test_db"

        try:
            settings = Settings()
            assert isinstance(settings.storage, MongoDBStorage)
            assert settings.storage.type == "mongodb"
            assert settings.storage.url == "mongodb://mongo.example.com:27017"
            assert settings.storage.database == "test_db"
        finally:
            os.environ.pop("ZNDRAW_STORAGE__TYPE", None)
            os.environ.pop("ZNDRAW_STORAGE__URL", None)
            os.environ.pop("ZNDRAW_STORAGE__DATABASE", None)


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
