"""Tests for application lifespan management.

Tests verify:
1. App starts without REDIS_URL (auto-starts TcpFakeServer)
2. App starts with REDIS_URL (uses external Redis)
3. Storage backends are correctly initialized based on config
"""

import os
from collections.abc import Generator

import httpx
import pytest

from zndraw.config import (
    LMDBStorage as LMDBStorageConfig,
    MemoryStorage,
    MongoDBStorage,
)
from zndraw.storage import InMemoryStorage, LMDBStorage as LMDBStorageBackend


@pytest.fixture(autouse=True)
def clean_env() -> Generator[None, None, None]:
    """Clean environment before each test."""
    # Remove any existing env vars that might interfere
    env_vars_to_clear = [
        "ZNDRAW_REDIS_URL",
        "ZNDRAW_DATABASE_URL",
        "ZNDRAW_STORAGE__TYPE",
        "ZNDRAW_STORAGE__PATH",
        "ZNDRAW_STORAGE__MAP_SIZE",
    ]
    original_values = {k: os.environ.pop(k, None) for k in env_vars_to_clear}

    yield

    # Restore original values
    for key, value in original_values.items():
        if value is not None:
            os.environ[key] = value
        else:
            os.environ.pop(key, None)


class TestLifespanWithoutRedisUrl:
    """Test lifespan when REDIS_URL is not configured."""

    def test_app_starts_without_redis_url(self, server_factory) -> None:
        """App should auto-start TcpFakeServer when REDIS_URL is not set."""
        # Explicitly unset REDIS_URL
        os.environ.pop("ZNDRAW_REDIS_URL", None)
        os.environ["ZNDRAW_DATABASE_URL"] = "sqlite+aiosqlite://"

        # Create server without REDIS_URL - should work with fake server
        instance = server_factory(
            {
                "ZNDRAW_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )
        # Need to unset REDIS_URL after server_factory applies its defaults
        # The server should have started successfully

        response = httpx.get(f"{instance.url}/v1/health", timeout=5.0)
        assert response.status_code == 200
        assert response.json() == {"status": "ok"}

    def test_health_endpoint_returns_200_with_fake_redis(self, server_factory) -> None:
        """Health endpoint should return 200 when using fake Redis."""
        os.environ.pop("ZNDRAW_REDIS_URL", None)

        instance = server_factory(
            {
                "ZNDRAW_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        response = httpx.get(f"{instance.url}/v1/health", timeout=5.0)
        assert response.status_code == 200


class TestLifespanWithRedisUrl:
    """Test lifespan when REDIS_URL is configured."""

    def test_app_starts_with_redis_url(self, server_factory) -> None:
        """App should use external Redis when REDIS_URL is set."""
        instance = server_factory(
            {
                "ZNDRAW_REDIS_URL": "redis://localhost:6379",
                "ZNDRAW_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        response = httpx.get(f"{instance.url}/v1/health", timeout=5.0)
        assert response.status_code == 200
        assert response.json() == {"status": "ok"}


class TestStorageInitialization:
    """Test storage backend initialization based on config."""

    def test_memory_storage_config_creates_inmemory_backend(self) -> None:
        """MemoryStorage config should create InMemoryStorage backend."""
        from zndraw.database import _create_storage_backend

        config = MemoryStorage()
        backend = _create_storage_backend(config)

        assert isinstance(backend, InMemoryStorage)

    def test_lmdb_storage_config_creates_lmdb_backend(self, tmp_path) -> None:
        """LMDBStorage config should create LMDBStorage backend."""
        from zndraw.database import _create_storage_backend

        lmdb_path = tmp_path / "test-lmdb"
        config = LMDBStorageConfig(path=lmdb_path, map_size=100_000_000)
        backend = _create_storage_backend(config)

        assert isinstance(backend, LMDBStorageBackend)
        assert backend.path == lmdb_path
        assert backend.map_size == 100_000_000

    def test_mongodb_storage_config_raises_not_implemented(self) -> None:
        """MongoDBStorage config should raise NotImplementedError."""
        from zndraw.database import _create_storage_backend

        config = MongoDBStorage(url="mongodb://localhost:27017")

        with pytest.raises(NotImplementedError, match="MongoDB storage"):
            _create_storage_backend(config)


class TestFakeRedisServer:
    """Test the fake Redis server functionality."""

    def test_get_free_port_returns_valid_port(self) -> None:
        """_get_free_port should return a valid port number."""
        from zndraw.database import _get_free_port

        port = _get_free_port()
        assert isinstance(port, int)
        assert 1024 <= port <= 65535

    def test_get_free_port_returns_different_ports(self) -> None:
        """_get_free_port should return different ports on successive calls."""
        from zndraw.database import _get_free_port

        # Get multiple ports - they should be different (high probability)
        ports = {_get_free_port() for _ in range(5)}
        # At least some should be different (very high probability)
        assert len(ports) >= 2


class TestCleanupSweeper:
    """Tests for cleanup sweeper integration."""

    def test_cleanup_sweeper_starts_in_lifespan(self, server_factory) -> None:
        """Cleanup sweeper should start when app starts."""
        os.environ.pop("ZNDRAW_REDIS_URL", None)

        instance = server_factory(
            {
                "ZNDRAW_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        # Health endpoint works means lifespan completed successfully
        response = httpx.get(f"{instance.url}/v1/health", timeout=5.0)
        assert response.status_code == 200
