"""Test fixtures for async database sessions and integration testing."""

import os
import socket
import threading
import time
from collections.abc import AsyncIterator, Callable, Generator
from contextlib import asynccontextmanager
from dataclasses import dataclass

import httpx
import pytest
import pytest_asyncio
import uvicorn
from helpers import InMemoryResultBackend, MockSioServer, create_test_user_model
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.config import Settings
from zndraw.storage import FrameStorage
from zndraw_auth import User
from zndraw_auth.settings import AuthSettings


@pytest.fixture(autouse=True)
def test_settings() -> None:
    """Configure test settings for each test.

    Sets test environment variables so that Settings() created anywhere
    (lifespan, CLI, etc.) picks up test values.
    Note: REDIS_URL is not set - the lifespan will auto-start TcpFakeServer.
    """
    # Database (in-memory for tests)
    os.environ["ZNDRAW_SERVER_DATABASE_URL"] = "sqlite+aiosqlite://"
    # Don't set REDIS_URL - lifespan will auto-start TcpFakeServer
    os.environ.pop("ZNDRAW_SERVER_REDIS_URL", None)
    # Remove host/port so Settings() reads its own defaults
    os.environ.pop("ZNDRAW_SERVER_HOST", None)
    os.environ.pop("ZNDRAW_SERVER_PORT", None)

    return


def get_free_port() -> int:
    """Find a free port on localhost."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


# =============================================================================
# REST API Test Fixtures (async sessions)
# =============================================================================


@pytest_asyncio.fixture(name="session")
async def session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session for each test."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )

    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)

        async_session_factory = async_sessionmaker(
            bind=engine,
            class_=AsyncSession,
            expire_on_commit=False,
        )

        async with async_session_factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="client")
async def client_fixture(
    session: AsyncSession,
    redis_client,
    mock_sio: MockSioServer,
    frame_storage: FrameStorage,
    result_backend: InMemoryResultBackend,
) -> AsyncIterator[AsyncClient]:
    """Async test client with real Redis, real DB, MockSioServer.

    All route tests share this fixture. No per-file client definitions.
    """
    from zndraw.app import app
    from zndraw.dependencies import (
        get_frame_storage,
        get_joblib_settings,
        get_redis,
        get_result_backend,
        get_tsio,
    )
    from zndraw_auth import get_session
    from zndraw_joblib.settings import JobLibSettings

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield session

    @asynccontextmanager
    async def test_session_maker():
        yield session

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_redis] = lambda: redis_client
    app.dependency_overrides[get_tsio] = lambda: mock_sio
    app.dependency_overrides[get_frame_storage] = lambda: frame_storage
    app.dependency_overrides[get_result_backend] = lambda: result_backend
    app.dependency_overrides[get_joblib_settings] = lambda: JobLibSettings()
    app.state.session_maker = test_session_maker
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

    # Ensure the internal worker user exists for WorkerTokenDep resolution.
    # Cache it on app.state so get_worker_token doesn't open its own session
    # (which would deadlock under SQLite serialization in routes that already
    # hold a yield-based SessionDep).
    from sqlmodel import select as _select

    from zndraw.database import ensure_internal_worker

    await ensure_internal_worker(session, app.state.settings.internal_worker_email)
    result = await session.exec(
        _select(User).where(User.email == app.state.settings.internal_worker_email)  # type: ignore[arg-type]
    )
    app.state.internal_worker_user = result.one()

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client

    app.dependency_overrides.clear()


@pytest_asyncio.fixture(name="test_user")
async def test_user_fixture(session: AsyncSession) -> User:
    """Create a test user in the database."""
    user = create_test_user_model()
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user


@pytest.fixture(name="mock_sio")
def mock_sio_fixture() -> MockSioServer:
    """MockSioServer for route tests — shared across client and test assertions."""
    return MockSioServer()


@pytest.fixture(name="settings")
def settings_fixture() -> Settings:
    """Return a Settings instance for the current test environment."""
    return Settings()


@pytest.fixture(name="result_backend")
def result_backend_fixture() -> InMemoryResultBackend:
    """In-memory result backend for route tests."""
    return InMemoryResultBackend()


# =============================================================================
# Redis Test Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="redis_client")
async def redis_client_fixture():
    """Create a Redis client for testing.

    Connects to localhost Redis and flushes the database before/after each test.
    """
    from redis.asyncio import Redis

    redis: Redis = Redis.from_url("redis://localhost", decode_responses=True)
    await redis.flushdb()
    yield redis
    await redis.flushdb()
    await redis.aclose()


@pytest_asyncio.fixture(name="frame_storage")
async def frame_storage_fixture(
    redis_client,
) -> AsyncIterator[FrameStorage]:
    """FrameStorage backed by real Redis (flushed per test)."""
    storage = FrameStorage("memory://", redis_client)
    yield storage
    await storage.close()


# =============================================================================
# Socket.IO Integration Test Fixtures (Factory Pattern)
# =============================================================================


@dataclass
class ServerInstance:
    """Holds a running server instance for cleanup."""

    url: str
    server: uvicorn.Server
    thread: threading.Thread
    env_overrides: dict[str, str]


# Type alias for server factory
ServerFactory = Callable[[dict[str, str]], ServerInstance]


@pytest.fixture(name="server_factory")
def server_factory_fixture() -> Generator[ServerFactory, None, None]:
    """Factory fixture that creates servers with custom settings.

    Usage:
        def test_something(server_factory):
            server = server_factory({"ZNDRAW_SERVER_EDIT_LOCK_TTL": "2"})
            # server.url contains the server URL

    The factory handles cleanup of all created servers automatically.
    """
    from zndraw.app import socket_app

    created_servers: list[ServerInstance] = []
    original_env: dict[str, str | None] = {}

    def _create_server(env_overrides: dict[str, str] | None = None) -> ServerInstance:
        # Always use real Redis for testing (FakeRedis is only for standalone mode)
        port = get_free_port()
        host = "127.0.0.1"

        defaults = {
            "ZNDRAW_SERVER_REDIS_URL": "redis://localhost",
            "ZNDRAW_SERVER_DATABASE_URL": "sqlite+aiosqlite://",
            "ZNDRAW_SERVER_HOST": host,
            "ZNDRAW_SERVER_PORT": str(port),
            "ZNDRAW_SERVER_TASK_QUEUE_NAME": f"zndraw:tasks:{port}",
            "ZNDRAW_SERVER_RESULT_BACKEND_KEY_PREFIX": f"zndraw:{port}",
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "false",
        }
        defaults.update(env_overrides or {})
        env_overrides = defaults

        # Store original values and apply overrides
        for key, value in env_overrides.items():
            if key not in original_env:
                original_env[key] = os.environ.get(key)
            os.environ[key] = value
        url = f"http://{host}:{port}"

        config = uvicorn.Config(
            socket_app,
            host=host,
            port=port,
            log_level="error",
        )
        server = uvicorn.Server(config)

        thread = threading.Thread(target=server.run)
        thread.daemon = True
        thread.start()

        # Wait for server to be ready
        timeout = 5.0
        start_time = time.time()
        while True:
            if time.time() - start_time > timeout:
                raise RuntimeError("Server timed out waiting for start")
            try:
                response = httpx.get(f"{url}/v1/health", timeout=1.0)
                if response.status_code == 200:
                    break
            except httpx.RequestError:
                time.sleep(0.1)

        instance = ServerInstance(
            url=url,
            server=server,
            thread=thread,
            env_overrides=env_overrides,
        )
        created_servers.append(instance)
        return instance

    yield _create_server

    # Cleanup all created servers
    for instance in created_servers:
        instance.server.should_exit = True
        instance.thread.join(timeout=1)

    # Restore original environment
    for key, original_value in original_env.items():
        if original_value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = original_value


@pytest.fixture(name="server")
def server_fixture(server_factory: ServerFactory) -> str:
    """Start a real uvicorn server for Socket.IO integration testing."""
    instance = server_factory({})
    return instance.url


@pytest.fixture(name="server_auth")
def server_auth_fixture(server_factory: ServerFactory) -> str:
    """Start a server with authentication enabled (production-like roles).

    Sets default admin email to enable production mode:
    - New users get regular user role
    - Only the admin user gets superuser
    """
    instance = server_factory(
        {
            "ZNDRAW_AUTH_DEFAULT_ADMIN_EMAIL": "admin@local.test",
            "ZNDRAW_AUTH_DEFAULT_ADMIN_PASSWORD": "adminpassword",
        }
    )
    return instance.url


@pytest_asyncio.fixture(name="http_client")
async def http_client_fixture(server: str) -> AsyncIterator[AsyncClient]:
    """Provide an async HTTP client for Socket.IO integration tests."""
    async with AsyncClient(base_url=server) as client:
        yield client


@pytest_asyncio.fixture(name="http_client_auth")
async def http_client_auth_fixture(
    server_auth: str,
) -> AsyncIterator[AsyncClient]:
    """Provide an async HTTP client for auth-enabled server."""
    async with AsyncClient(base_url=server_auth) as client:
        yield client
