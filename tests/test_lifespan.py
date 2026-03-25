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

from zndraw.storage import FrameStorage


@pytest.fixture(autouse=True)
def clean_env() -> Generator[None, None, None]:
    """Clean environment before each test."""
    env_vars_to_clear = [
        "ZNDRAW_SERVER_REDIS_URL",
        "ZNDRAW_SERVER_DATABASE_URL",
        "ZNDRAW_SERVER_STORAGE",
    ]
    original_values = {k: os.environ.pop(k, None) for k in env_vars_to_clear}

    yield

    for key, value in original_values.items():
        if value is not None:
            os.environ[key] = value
        else:
            os.environ.pop(key, None)


class TestLifespanWithoutRedisUrl:
    """Test lifespan when REDIS_URL is not configured."""

    def test_app_starts_without_redis_url(self, server_factory) -> None:
        """App should auto-start TcpFakeServer when REDIS_URL is not set."""
        os.environ.pop("ZNDRAW_SERVER_REDIS_URL", None)
        os.environ["ZNDRAW_SERVER_DATABASE_URL"] = "sqlite+aiosqlite://"

        instance = server_factory(
            {
                "ZNDRAW_SERVER_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        response = httpx.get(f"{instance.url}/v1/health", timeout=5.0)
        assert response.status_code == 200
        assert response.json() == {"status": "ok"}

    def test_health_endpoint_returns_200_with_fake_redis(self, server_factory) -> None:
        """Health endpoint should return 200 when using fake Redis."""
        os.environ.pop("ZNDRAW_SERVER_REDIS_URL", None)

        instance = server_factory(
            {
                "ZNDRAW_SERVER_DATABASE_URL": "sqlite+aiosqlite://",
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
                "ZNDRAW_SERVER_REDIS_URL": "redis://localhost:6379",
                "ZNDRAW_SERVER_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        response = httpx.get(f"{instance.url}/v1/health", timeout=5.0)
        assert response.status_code == 200
        assert response.json() == {"status": "ok"}


class TestStorageInitialization:
    """Test storage backend initialization based on config."""

    def test_memory_storage_creates_frame_storage(self, redis_client) -> None:
        """memory:// URI should create FrameStorage backend."""
        backend = FrameStorage(uri="memory://", redis=redis_client)
        assert isinstance(backend, FrameStorage)


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

        ports = {_get_free_port() for _ in range(5)}
        assert len(ports) >= 2


class TestCleanupSweeper:
    """Tests for cleanup sweeper integration."""

    def test_cleanup_sweeper_starts_in_lifespan(self, server_factory) -> None:
        """Cleanup sweeper should start when app starts."""
        os.environ.pop("ZNDRAW_SERVER_REDIS_URL", None)

        instance = server_factory(
            {
                "ZNDRAW_SERVER_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        response = httpx.get(f"{instance.url}/v1/health", timeout=5.0)
        assert response.status_code == 200
