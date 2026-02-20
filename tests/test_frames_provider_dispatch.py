"""Tests for provider-aware frame dispatch in GET /v1/rooms/{room_id}/frames/{index}."""

import asyncio
from collections.abc import AsyncIterator
from unittest.mock import AsyncMock

import msgpack
import pytest
import pytest_asyncio
from conftest import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
    decode_msgpack_response,
)
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from zndraw_auth import User
from zndraw_joblib.dependencies import request_hash
from zndraw_joblib.exceptions import ProviderTimeout
from zndraw_joblib.models import ProviderRecord, Worker

from zndraw.exceptions import FrameNotFound, ProblemDetail
from zndraw.storage import InMemoryStorage
from zndraw.storage.base import RawFrame

# =============================================================================
# In-memory ResultBackend for testing
# =============================================================================


class InMemoryResultBackend:
    """Minimal ResultBackend for testing provider dispatch."""

    def __init__(self) -> None:
        self._data: dict[str, bytes] = {}
        self._inflight: set[str] = set()
        self._waiters: dict[str, list[asyncio.Event]] = {}

    async def store(self, key: str, data: bytes, ttl: int) -> None:
        self._data[key] = data

    async def get(self, key: str) -> bytes | None:
        return self._data.get(key)

    async def delete(self, key: str) -> None:
        self._data.pop(key, None)

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        if key in self._inflight:
            return False
        self._inflight.add(key)
        return True

    async def release_inflight(self, key: str) -> None:
        self._inflight.discard(key)

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:
        cached = self._data.get(key)
        if cached is not None:
            return cached
        event = asyncio.Event()
        self._waiters.setdefault(key, []).append(event)
        try:
            await asyncio.wait_for(event.wait(), timeout=timeout)
            return self._data.get(key)
        except asyncio.TimeoutError:
            return None
        finally:
            waiters = self._waiters.get(key, [])
            if event in waiters:
                waiters.remove(event)

    async def notify_key(self, key: str) -> None:
        for event in self._waiters.pop(key, []):
            event.set()


# =============================================================================
# Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="prov_session_factory")
async def prov_session_factory_fixture() -> AsyncIterator[
    async_sessionmaker[AsyncSession]
]:
    """Create a session factory backed by a fresh in-memory database."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)
        yield async_sessionmaker(
            bind=engine, class_=AsyncSession, expire_on_commit=False
        )
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="prov_session")
async def prov_session_fixture(
    prov_session_factory: async_sessionmaker[AsyncSession],
) -> AsyncIterator[AsyncSession]:
    """Create a session from the shared factory."""
    async with prov_session_factory() as session:
        yield session


@pytest_asyncio.fixture(name="prov_storage")
async def prov_storage_fixture() -> AsyncIterator[InMemoryStorage]:
    """Create a fresh InMemoryStorage."""
    storage = InMemoryStorage()
    yield storage
    await storage.close()


@pytest.fixture(name="prov_result_backend")
def prov_result_backend_fixture() -> InMemoryResultBackend:
    """Create an InMemoryResultBackend."""
    return InMemoryResultBackend()


@pytest_asyncio.fixture(name="prov_client")
async def prov_client_fixture(
    prov_session: AsyncSession,
    prov_session_factory: async_sessionmaker[AsyncSession],
    prov_storage: InMemoryStorage,
    prov_result_backend: InMemoryResultBackend,
) -> AsyncIterator[AsyncClient]:
    """Create a test client with provider dependencies wired."""
    from zndraw_auth import get_session
    from zndraw_auth.db import get_session_maker
    from zndraw_auth.settings import AuthSettings

    from zndraw.app import app
    from zndraw.dependencies import (
        get_redis,
        get_result_backend,
        get_storage,
        get_tsio,
    )

    mock_sio = MockSioServer()
    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield prov_session

    app.state.auth_settings = AuthSettings()
    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_session_maker] = lambda: prov_session_factory
    app.dependency_overrides[get_storage] = lambda: prov_storage
    app.dependency_overrides[get_tsio] = lambda: mock_sio
    app.dependency_overrides[get_redis] = lambda: mock_redis
    app.dependency_overrides[get_result_backend] = lambda: prov_result_backend

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client

    app.dependency_overrides.clear()


# =============================================================================
# Helpers
# =============================================================================


_create_user = create_test_user_in_db
_create_room = create_test_room
_auth = auth_header


async def _create_provider(
    session: AsyncSession, room_id: str, user: User
) -> ProviderRecord:
    """Create a frames provider for the given room."""
    worker = Worker(user_id=user.id)
    session.add(worker)
    await session.flush()
    provider = ProviderRecord(
        room_id=room_id,
        category="frames",
        name=f"mount-{room_id}",
        content_type="application/x-msgpack",
        user_id=user.id,
        worker_id=worker.id,
    )
    session.add(provider)
    await session.commit()
    await session.refresh(provider)
    return provider


# =============================================================================
# Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_frame_storage_hit_ignores_provider(
    prov_client: AsyncClient,
    prov_session: AsyncSession,
    prov_storage: InMemoryStorage,
) -> None:
    """Frame in storage, provider exists -- returns frame (no provider touch)."""
    user, token = await _create_user(prov_session)
    room = await _create_room(prov_session, user)
    await _create_provider(prov_session, room.id, user)

    await prov_storage.extend(room.id, [{"a": 1}])

    response = await prov_client.get(
        f"/v1/rooms/{room.id}/frames/0", headers=_auth(token)
    )
    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert msgpack.unpackb(frames[0][b"a"]) == 1


@pytest.mark.asyncio
async def test_get_frame_provider_cache_hit(
    prov_client: AsyncClient,
    prov_session: AsyncSession,
    prov_storage: InMemoryStorage,
    prov_result_backend: InMemoryResultBackend,
) -> None:
    """Frame in provider cache, storage slot is None -- returns 200 with frame."""
    user, token = await _create_user(prov_session)
    room = await _create_room(prov_session, user)
    provider = await _create_provider(prov_session, room.id, user)

    # Reserve slots (provider has 3 frames), slot 0 is None
    await prov_storage.reserve(room.id, 3)

    # Pre-populate provider cache
    positions_packed = msgpack.packb([1, 2, 3])
    assert isinstance(positions_packed, bytes)
    frame_data: RawFrame = {b"arrays.positions": positions_packed}
    packed = msgpack.packb(frame_data, use_bin_type=True)
    assert isinstance(packed, bytes)
    params = {"index": "0"}
    rhash = request_hash(params)
    cache_key = f"provider-result:{provider.full_name}:{rhash}"
    await prov_result_backend.store(cache_key, packed, ttl=300)

    response = await prov_client.get(
        f"/v1/rooms/{room.id}/frames/0", headers=_auth(token)
    )
    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert b"arrays.positions" in frames[0]


@pytest.mark.asyncio
async def test_get_frame_provider_timeout(
    prov_client: AsyncClient,
    prov_session: AsyncSession,
    prov_storage: InMemoryStorage,
) -> None:
    """Frame not cached, provider exists -- long-poll times out → 504."""
    user, token = await _create_user(prov_session)
    room = await _create_room(prov_session, user)
    await _create_provider(prov_session, room.id, user)

    # Reserve slots, leave them empty
    await prov_storage.reserve(room.id, 5)

    response = await prov_client.get(
        f"/v1/rooms/{room.id}/frames/2", headers=_auth(token)
    )
    assert response.status_code == 504

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == ProviderTimeout.type_uri()
    assert response.headers.get("retry-after") == "2"


@pytest.mark.asyncio
async def test_get_frame_no_provider_returns_404(
    prov_client: AsyncClient,
    prov_session: AsyncSession,
    prov_storage: InMemoryStorage,
) -> None:
    """Frame missing, no provider registered -- returns 404."""
    user, token = await _create_user(prov_session)
    room = await _create_room(prov_session, user)

    # Reserve slots but no provider registered
    await prov_storage.reserve(room.id, 3)

    response = await prov_client.get(
        f"/v1/rooms/{room.id}/frames/1", headers=_auth(token)
    )
    assert response.status_code == 404
    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_frame_dispatch_acquires_inflight(
    prov_client: AsyncClient,
    prov_session: AsyncSession,
    prov_storage: InMemoryStorage,
    prov_result_backend: InMemoryResultBackend,
) -> None:
    """After dispatch, inflight lock is acquired."""
    user, token = await _create_user(prov_session)
    room = await _create_room(prov_session, user)
    provider = await _create_provider(prov_session, room.id, user)

    await prov_storage.reserve(room.id, 3)

    response = await prov_client.get(
        f"/v1/rooms/{room.id}/frames/0", headers=_auth(token)
    )
    assert response.status_code == 504  # timeout, but inflight lock was set

    # Verify inflight lock was acquired
    params = {"index": "0"}
    rhash = request_hash(params)
    inflight_key = f"provider-inflight:{provider.full_name}:{rhash}"
    assert inflight_key in prov_result_backend._inflight


@pytest.mark.asyncio
async def test_get_frame_notify_wakes_long_poll(
    prov_client: AsyncClient,
    prov_session: AsyncSession,
    prov_storage: InMemoryStorage,
    prov_result_backend: InMemoryResultBackend,
) -> None:
    """Provider uploads result mid-poll — long-poll wakes up and returns 200."""
    user, token = await _create_user(prov_session)
    room = await _create_room(prov_session, user)
    provider = await _create_provider(prov_session, room.id, user)

    await prov_storage.reserve(room.id, 3)

    # Build the cache key that _dispatch_provider_frame will wait on
    params = {"index": "1"}
    rhash = request_hash(params)
    cache_key = f"provider-result:{provider.full_name}:{rhash}"

    # Prepare the frame data that the "provider" will upload
    positions_packed = msgpack.packb([10, 20, 30])
    assert isinstance(positions_packed, bytes)
    frame_data: RawFrame = {b"arrays.positions": positions_packed}
    packed = msgpack.packb(frame_data, use_bin_type=True)
    assert isinstance(packed, bytes)

    async def _simulate_provider_upload() -> None:
        """Simulate provider uploading result after a short delay."""
        await asyncio.sleep(0.1)
        await prov_result_backend.store(cache_key, packed, ttl=300)
        await prov_result_backend.notify_key(cache_key)

    # Start the simulated provider upload concurrently with the GET request
    upload_task = asyncio.create_task(_simulate_provider_upload())

    response = await prov_client.get(
        f"/v1/rooms/{room.id}/frames/1", headers=_auth(token)
    )

    await upload_task

    # The long-poll should have woken up and returned 200
    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert b"arrays.positions" in frames[0]
    assert msgpack.unpackb(frames[0][b"arrays.positions"]) == [10, 20, 30]


@pytest.mark.asyncio
async def test_list_frames_notify_wakes_concurrent_dispatch(
    prov_client: AsyncClient,
    prov_session: AsyncSession,
    prov_storage: InMemoryStorage,
    prov_result_backend: InMemoryResultBackend,
) -> None:
    """Multiple missing frames dispatched concurrently — all wake on notify."""
    user, token = await _create_user(prov_session)
    room = await _create_room(prov_session, user)
    provider = await _create_provider(prov_session, room.id, user)

    # Reserve 3 slots, fill only index 1
    await prov_storage.reserve(room.id, 3)
    await prov_storage.set_item(room.id, 1, {"existing": 1})

    # Build cache keys for the missing frames (index 0 and 2)
    def _cache_key(idx: int) -> str:
        return (
            f"provider-result:{provider.full_name}:{request_hash({'index': str(idx)})}"
        )

    async def _simulate_provider_upload() -> None:
        """Upload results for both missing frames after a short delay."""
        await asyncio.sleep(0.1)
        for idx in (0, 2):
            idx_packed = msgpack.packb(idx)
            assert isinstance(idx_packed, bytes)
            frame: RawFrame = {b"idx": idx_packed}
            packed = msgpack.packb(frame, use_bin_type=True)
            assert isinstance(packed, bytes)
            key = _cache_key(idx)
            await prov_result_backend.store(key, packed, ttl=300)
            await prov_result_backend.notify_key(key)

    upload_task = asyncio.create_task(_simulate_provider_upload())

    response = await prov_client.get(
        f"/v1/rooms/{room.id}/frames?indices=0,1,2", headers=_auth(token)
    )

    await upload_task

    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 3
    # Index 0 and 2 came from provider, index 1 from storage
    assert msgpack.unpackb(frames[0][b"idx"]) == 0
    assert msgpack.unpackb(frames[1][b"existing"]) == 1
    assert msgpack.unpackb(frames[2][b"idx"]) == 2
