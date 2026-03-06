"""Test fixtures for async database sessions and integration testing."""

import os
import socket
import threading
import time
from collections.abc import AsyncIterator, Callable, Generator
from contextlib import asynccontextmanager
from dataclasses import dataclass
from typing import Any

import httpx
import msgpack
import pytest
import pytest_asyncio
import uvicorn
from fastapi_users.password import PasswordHelper
from httpx import ASGITransport, AsyncClient
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from zndraw_auth import User
from zndraw_auth.settings import AuthSettings

from zndraw.config import Settings
from zndraw.models import Room
from zndraw.storage.base import RawFrame

_password_helper = PasswordHelper()


def make_raw_frame(data: dict) -> RawFrame:
    """Convert a simple dict to RawFrame format for test assertions.

    The storage converts input dicts to raw bytes format.
    This helper creates the expected output for simple test dicts.
    """
    result: RawFrame = {}
    for k, v in data.items():
        key: bytes = k.encode() if isinstance(k, str) else k
        packed = msgpack.packb(v)
        val: bytes = packed if packed is not None else b""
        result[key] = val
    return result


def create_test_user_model(
    email: str = "testuser@local.test",
    password: str = "testpassword",
    is_superuser: bool = False,
) -> User:
    """Create a User model instance with hashed password for tests."""
    return User(
        email=email,
        hashed_password=_password_helper.hash(password),
        is_active=True,
        is_superuser=is_superuser,
        is_verified=True,
    )


def create_test_token(user: User) -> str:
    """Create a JWT token for a test user."""
    from fastapi_users.jwt import generate_jwt

    settings = AuthSettings()
    data = {"sub": str(user.id), "aud": "fastapi-users:auth"}
    return generate_jwt(
        data,
        settings.secret_key.get_secret_value(),
        settings.token_lifetime_seconds,
    )


def decode_msgpack_response(content: bytes) -> list[dict[bytes, bytes]]:
    """Decode a MessagePack response body to a list of raw frames."""
    return msgpack.unpackb(content, raw=True)


async def create_test_user_in_db(
    session: AsyncSession, email: str = "testuser@local.test"
) -> tuple[User, str]:
    """Create a user in the DB and return (user, token)."""
    user = create_test_user_model(email=email)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user, create_test_token(user)


async def create_test_room(
    session: AsyncSession, user: User, description: str = "Test Room"
) -> Room:
    """Create a room with user as owner and return it."""
    from zndraw.models import MemberRole, RoomMembership

    room = Room(
        description=description,
        created_by_id=user.id,  # type: ignore[arg-type]
        is_public=True,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user.id,  # type: ignore[arg-type]
        role=MemberRole.OWNER,
    )
    session.add(membership)
    await session.commit()
    return room


def auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


@pytest.fixture(autouse=True)
def test_settings() -> Generator[None, None, None]:
    """Configure test settings for each test.

    Sets test environment variables so that Settings() created anywhere
    (lifespan, CLI, etc.) picks up test values.
    Note: REDIS_URL is not set - the lifespan will auto-start TcpFakeServer.
    """
    # Database (in-memory for tests)
    os.environ["ZNDRAW_DATABASE_URL"] = "sqlite+aiosqlite://"
    # Don't set REDIS_URL - lifespan will auto-start TcpFakeServer
    os.environ.pop("ZNDRAW_REDIS_URL", None)
    # Remove host/port so Settings() reads its own defaults
    os.environ.pop("ZNDRAW_HOST", None)
    os.environ.pop("ZNDRAW_PORT", None)

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


class MockSioServer:
    """Mock Socket.IO server for testing broadcasts.

    Compatible with zndraw-socketio's AsyncServerWrapper emit pattern:
    - sio.emit(PydanticModel(), room=...) - model as first arg
    - sio.emit("event", data, room=...) - classic pattern
    """

    def __init__(self) -> None:
        self.emitted: list[dict[str, Any]] = []
        self.rooms: dict[str, set[str]] = {}
        self.sessions: dict[str, dict[str, Any]] = {}

    async def emit(
        self,
        event_or_model: str | BaseModel,
        data: Any = None,
        *,
        room: str | None = None,
        skip_sid: str | None = None,
        to: str | None = None,
        **kwargs: Any,
    ) -> None:
        if isinstance(event_or_model, BaseModel):
            # zndraw-socketio pattern: model class name -> snake_case event
            cls_name = type(event_or_model).__name__
            event = "".join(
                f"_{c.lower()}" if c.isupper() else c for c in cls_name
            ).lstrip("_")
            data = event_or_model.model_dump()
        else:
            event = event_or_model
        self.emitted.append({"event": event, "data": data, "room": room, "to": to})

    async def enter_room(self, sid: str, room: str) -> None:
        if room not in self.rooms:
            self.rooms[room] = set()
        self.rooms[room].add(sid)

    async def get_session(self, sid: str) -> dict[str, Any]:
        return self.sessions.get(sid, {})

    async def save_session(self, sid: str, session: dict[str, Any]) -> None:
        self.sessions[sid] = session


@pytest_asyncio.fixture(name="client")
async def client_fixture(session: AsyncSession) -> AsyncIterator[AsyncClient]:
    """Create an async test client with the session dependency overridden."""
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    mock_sio = MockSioServer()

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield session

    def get_sio_override() -> MockSioServer:
        return mock_sio

    # Create test session_maker for Socket.IO handlers
    @asynccontextmanager
    async def test_session_maker():
        yield session

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_redis] = lambda: (
        None
    )  # tests that need redis override this
    app.dependency_overrides[get_tsio] = get_sio_override
    app.state.session_maker = test_session_maker
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

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
            server = server_factory({"ZNDRAW_PRESENCE_TTL": "2"})
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
            "ZNDRAW_REDIS_URL": "redis://localhost",
            "ZNDRAW_DATABASE_URL": "sqlite+aiosqlite://",
            "ZNDRAW_HOST": host,
            "ZNDRAW_PORT": str(port),
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


@pytest.fixture(name="server_short_ttl")
def server_short_ttl_fixture(server_factory: ServerFactory) -> str:
    """Start a server with short presence TTL (2s) for testing expiration."""
    instance = server_factory({"ZNDRAW_PRESENCE_TTL": "2"})
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


@pytest_asyncio.fixture(name="http_client_short_ttl")
async def http_client_short_ttl_fixture(
    server_short_ttl: str,
) -> AsyncIterator[AsyncClient]:
    """Provide an async HTTP client for short TTL server."""
    async with AsyncClient(base_url=server_short_ttl) as client:
        yield client
