# tests/conftest.py
"""Shared test fixtures for DRY tests."""

import asyncio
import uuid
from collections.abc import AsyncGenerator
from contextlib import asynccontextmanager
from typing import ClassVar
from unittest.mock import MagicMock

import httpx
import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw_auth import Base, User
from zndraw_joblib.dependencies import get_result_backend
from zndraw_joblib.exceptions import ProblemError, problem_exception_handler
from zndraw_joblib.provider import Provider
from zndraw_joblib.router import router
from zndraw_joblib.settings import JobLibSettings


@pytest.fixture(scope="session")
def anyio_backend():
    return "asyncio"


@pytest.fixture
def test_user_id():
    """Fixed UUID for test user."""
    return uuid.UUID("12345678-1234-5678-1234-567812345678")


@pytest.fixture
def test_user(test_user_id):
    """Create a mock User for testing."""
    user = MagicMock(spec=User)
    user.id = test_user_id
    user.email = "test@example.com"
    user.is_active = True
    user.is_superuser = True
    user.is_verified = True
    return user


@pytest.fixture
async def async_engine():
    """Create an async in-memory SQLite database engine."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    yield engine
    await engine.dispose()


@pytest.fixture
async def async_session_factory(async_engine):
    """Create tables and return an async session factory."""
    async with async_engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)

    return async_sessionmaker(async_engine, class_=AsyncSession, expire_on_commit=False)


@pytest.fixture
def db_session(async_session_factory):
    """Return an async session generator for the sweeper and other non-DI uses."""

    async def get_test_session() -> AsyncGenerator[AsyncSession, None]:
        async with async_session_factory() as session:
            yield session

    return get_test_session


@pytest.fixture
def mock_current_user(test_user):
    """Factory that returns the current user dependency."""

    async def get_current_user():
        return test_user

    return get_current_user


class InMemoryResultBackend:
    """In-memory result backend for testing."""

    def __init__(self):
        self._store: dict[str, bytes] = {}
        self._inflight: set[str] = set()
        self._waiters: dict[str, list[asyncio.Event]] = {}

    async def store(self, key: str, data: bytes, _ttl: int) -> None:
        self._store[key] = data

    async def get(self, key: str) -> bytes | None:
        return self._store.get(key)

    async def delete(self, key: str) -> None:
        self._store.pop(key, None)

    async def acquire_inflight(self, key: str, _ttl: int) -> bool:
        if key in self._inflight:
            return False
        self._inflight.add(key)
        return True

    async def release_inflight(self, key: str) -> None:
        self._inflight.discard(key)

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:  # noqa: ASYNC109
        cached = self._store.get(key)
        if cached is not None:
            return cached
        event = asyncio.Event()
        self._waiters.setdefault(key, []).append(event)
        try:
            await asyncio.wait_for(event.wait(), timeout=timeout)
            return self._store.get(key)
        except TimeoutError:
            return None
        finally:
            waiters = self._waiters.get(key, [])
            if event in waiters:
                waiters.remove(event)

    async def notify_key(self, key: str) -> None:
        for event in self._waiters.pop(key, []):
            event.set()


class _MockClientApi:
    """Adapter to make TestClient work with JobManager's ApiManager protocol."""

    def __init__(self, test_client):
        self._client = test_client

    @property
    def http(self):
        return self._client

    @property
    def base_url(self) -> str:
        return ""

    def get_headers(self) -> dict[str, str]:
        return {}

    def raise_for_status(self, response) -> None:
        """Conform to ApiManager exception contract.

        Raises
        ------
        KeyError
            404 responses.
        PermissionError
            401/403 responses.
        ValueError
            409/422 responses.
        """
        if response.status_code < 400:
            return
        try:
            detail = response.json().get("detail", response.text)
        except Exception:  # noqa: BLE001
            detail = response.text
        if response.status_code == 404:
            raise KeyError(str(detail))
        if response.status_code in {401, 403}:
            raise PermissionError(str(detail))
        if response.status_code in {409, 422}:
            raise ValueError(str(detail))
        response.raise_for_status()


class _FsProvider(Provider):
    """Test provider for filesystem reads."""

    category: ClassVar[str] = "filesystem"
    path: str = "/"

    def read(self, handler):
        return handler.list_dir(self.path)


@pytest.fixture
def mock_client_api():
    """Return the MockClientApi class for wrapping TestClient as ApiManager."""
    return _MockClientApi


@pytest.fixture
def api(client):
    """MockClientApi wrapping the default test client."""
    return _MockClientApi(client)


@pytest.fixture
def fs_provider():
    """Return the FsProvider class for provider tests."""
    return _FsProvider


class _SettingsStub:
    """Minimal stub for the host-app Settings object.

    Only the attributes read by ``mint_internal_worker_token`` need to be
    present. Tests that exercise the @internal dispatch path use this to
    avoid importing ``zndraw.config.Settings`` (which is a host-app concern
    and would require full app wiring).
    """

    internal_worker_email: str = "worker@internal.user"


def _build_app(
    *,
    session_maker,
    current_user,
) -> FastAPI:
    """Build a configured FastAPI app with standard dependency overrides."""
    from zndraw_auth import (
        current_active_user,
        current_superuser,
        current_user_scoped_session,
    )
    from zndraw_auth.db import get_session_maker
    from zndraw_auth.settings import AuthSettings

    app = FastAPI()
    app.include_router(router)
    app.add_exception_handler(ProblemError, problem_exception_handler)
    app.dependency_overrides[get_session_maker] = lambda: session_maker
    app.dependency_overrides[current_active_user] = current_user
    app.dependency_overrides[current_superuser] = current_user
    app.dependency_overrides[current_user_scoped_session] = current_user
    app.state.joblib_settings = JobLibSettings()
    result_backend = InMemoryResultBackend()
    app.dependency_overrides[get_result_backend] = lambda: result_backend
    # WorkerTokenDep is still used by submit_task; stub it for @global tests.
    from zndraw_joblib.dependencies import get_worker_token

    app.dependency_overrides[get_worker_token] = lambda: "test-worker-token"

    # Set up app.state attributes consumed by mint_internal_worker_token and
    # get_internal_provider_registry. Per the "always set in lifespan" contract,
    # every attribute accessed by request-time deps must be present (value or
    # None) — no getattr fallbacks.
    app.state.settings = _SettingsStub()
    app.state.auth_settings = AuthSettings()
    _worker_user = MagicMock(spec=User)
    _worker_user.id = uuid.UUID("aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa")
    app.state.internal_worker_user = _worker_user
    app.state.internal_provider_registry = None
    return app


@pytest.fixture
def app(async_session_factory, mock_current_user):
    """Create a FastAPI app with dependency overrides."""
    return _build_app(
        session_maker=async_session_factory,
        current_user=mock_current_user,
    )


@pytest.fixture
def client(app):
    """Create a test client for the app."""
    return TestClient(app)


@pytest.fixture
def seeded_client(client):
    """Client with a pre-registered job for testing."""
    resp = client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    # Store the worker_id for tests that need it
    worker_id = resp.json().get("worker_id")
    client.seeded_worker_id = worker_id
    return client


@pytest.fixture
def client_factory(async_session_factory):
    """Factory to create clients with different user identities."""

    def create_client(identity: str, is_superuser: bool = True) -> TestClient:
        user_id = uuid.uuid5(uuid.NAMESPACE_DNS, identity)

        user = MagicMock(spec=User)
        user.id = user_id
        user.email = f"{identity}@example.com"
        user.is_active = True
        user.is_superuser = is_superuser
        user.is_verified = True

        async def get_current_user():
            return user

        app = _build_app(
            session_maker=async_session_factory,
            current_user=get_current_user,
        )

        test_client = TestClient(app)
        test_client.user_id = user_id
        return test_client

    return create_client


@pytest.fixture
async def threadsafe_engine(tmp_path):
    """File-based SQLite engine for tests with background threads.

    Unlike the in-memory StaticPool engine, this gives each thread its own
    connection so SQLite's internal locking handles concurrent access.
    """
    engine = create_async_engine(
        f"sqlite+aiosqlite:///{tmp_path / 'test.db'}",
    )
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)
    yield engine
    await engine.dispose()


@pytest.fixture
def threadsafe_client(threadsafe_engine, mock_current_user):
    """TestClient safe for multi-threaded tests (background claim loop, etc.)."""
    factory = async_sessionmaker(
        threadsafe_engine, class_=AsyncSession, expire_on_commit=False
    )
    app = _build_app(session_maker=factory, current_user=mock_current_user)
    return TestClient(app)


@pytest.fixture
async def async_client(async_session_factory, mock_current_user):
    """Async HTTP client with SQLite locking for concurrent stress testing."""
    db_lock = asyncio.Lock()

    @asynccontextmanager
    async def locked_session_maker():
        async with db_lock, async_session_factory() as session:
            yield session

    app = _build_app(
        session_maker=locked_session_maker,
        current_user=mock_current_user,
    )
    async with httpx.AsyncClient(
        transport=httpx.ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client


@pytest.fixture
def unguarded_client_factory(async_session_factory):
    """Factory to create clients without the get_worker_token override.

    Unlike ``client_factory``, this fixture does NOT stub out
    ``get_worker_token``, so the real dependency runs and reads
    ``app.state.internal_worker_user``. Intended for regression tests
    that verify behaviour when that cache is absent.
    """

    def create_client(
        identity: str,
        is_superuser: bool = True,
        *,
        internal_worker_user=None,
    ) -> TestClient:
        user_id = uuid.uuid5(uuid.NAMESPACE_DNS, identity)

        user = MagicMock(spec=User)
        user.id = user_id
        user.email = f"{identity}@example.com"
        user.is_active = True
        user.is_superuser = is_superuser
        user.is_verified = True

        async def get_current_user():
            return user

        from zndraw_auth import (
            current_active_user,
            current_superuser,
            current_user_scoped_session,
        )
        from zndraw_auth.db import get_session_maker

        app = FastAPI()
        app.include_router(router)
        app.add_exception_handler(ProblemError, problem_exception_handler)
        app.dependency_overrides[get_session_maker] = lambda: async_session_factory
        app.dependency_overrides[current_active_user] = get_current_user
        app.dependency_overrides[current_superuser] = get_current_user
        app.dependency_overrides[current_user_scoped_session] = get_current_user
        app.state.joblib_settings = JobLibSettings()
        result_backend = InMemoryResultBackend()
        app.dependency_overrides[get_result_backend] = lambda: result_backend
        # Intentionally NOT overriding get_worker_token — real dep runs.
        app.state.internal_worker_user = internal_worker_user
        app.state.internal_provider_registry = None
        app.state.settings = _SettingsStub()

        test_client = TestClient(app)
        test_client.user_id = user_id
        return test_client

    return create_client
