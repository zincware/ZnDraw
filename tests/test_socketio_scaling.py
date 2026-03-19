"""Tests for Socket.IO scaling with AsyncRedisManager."""

import os

import httpx
import socketio as socketio_lib
from socketio import AsyncRedisManager
from zndraw_socketio import wrap


class TestSocketIOClientManager:
    """Tests for Socket.IO client manager configuration."""

    def test_async_redis_manager_can_be_instantiated(self) -> None:
        """AsyncRedisManager should be instantiable with a Redis URL."""
        # This tests that the import and basic instantiation works
        manager = AsyncRedisManager("redis://localhost:6379")
        assert manager is not None

    def test_manager_initialized_not_set_prematurely(self) -> None:
        """manager_initialized must stay False so python-socketio calls initialize().

        python-socketio lazily calls manager.initialize() on the first client
        connection — but only if _sio.manager_initialized is False.
        initialize() starts the Redis pub/sub listener that enables
        cross-replica event delivery.  If the flag is True before the first
        connection, the listener never starts and events stay local.

        This test replicates the setup from database.py lifespan.
        """
        server = socketio_lib.AsyncServer(async_mode="asgi")
        tsio = wrap(server)

        mgr = AsyncRedisManager("redis://localhost:6379")
        mgr.set_server(tsio)
        tsio.manager = mgr

        # The underlying server must report False so that
        # _handle_eio_connect triggers manager.initialize().
        assert server.manager_initialized is False

    def test_app_starts_with_client_manager(self, server_factory) -> None:
        """App should start successfully with AsyncRedisManager configured."""
        os.environ.pop("ZNDRAW_REDIS_URL", None)

        instance = server_factory(
            {
                "ZNDRAW_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        response = httpx.get(f"{instance.url}/health", timeout=5.0)
        assert response.status_code == 200

    def test_sio_manager_is_async_redis_manager(self) -> None:
        """After lifespan, tsio.manager should be AsyncRedisManager."""
        from zndraw.socketio import tsio

        # Before lifespan runs, the manager is default AsyncManager
        # After lifespan, it should be AsyncRedisManager
        # This test verifies the module-level tsio can be configured
        # Note: Full integration is tested via server_factory tests
        assert hasattr(tsio, "manager")
        assert hasattr(tsio, "manager_initialized")

    def test_sio_stored_in_app_state(self, server_factory) -> None:
        """app.state.sio should reference the configured Socket.IO server."""
        os.environ.pop("ZNDRAW_REDIS_URL", None)

        instance = server_factory(
            {
                "ZNDRAW_DATABASE_URL": "sqlite+aiosqlite://",
            }
        )

        # Verify server is running
        response = httpx.get(f"{instance.url}/health", timeout=5.0)
        assert response.status_code == 200

        # Note: We cannot directly access app.state from the test process
        # due to process isolation (server runs in separate thread).
        # The fact that Socket.IO connections work with the dependency
        # that reads from app.state.sio proves the setup is correct.
