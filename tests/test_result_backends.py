"""Unit tests for StorageResultBackend, CompositeResultBackend, and RedisResultBackend."""

import asyncio
from unittest.mock import AsyncMock

import pytest
import pytest_asyncio
from redis.asyncio import Redis

from zndraw.result_backends import (
    CompositeResultBackend,
    RedisResultBackend,
    StorageResultBackend,
)
from zndraw.storage.memory import InMemoryStorage

# =============================================================================
# StorageResultBackend
# =============================================================================


@pytest.fixture
def storage_backend():
    """Fresh InMemoryStorage for each test."""
    return InMemoryStorage()


@pytest.fixture
def backend(storage_backend):
    """StorageResultBackend wrapping InMemoryStorage."""
    return StorageResultBackend(storage_backend)


@pytest.mark.asyncio
async def test_store_and_get(backend):
    """Stored data is retrievable."""
    await backend.store("key1", b"hello", ttl=60)
    assert await backend.get("key1") == b"hello"


@pytest.mark.asyncio
async def test_get_missing_returns_none(backend):
    """Missing key returns None."""
    assert await backend.get("nonexistent") is None


@pytest.mark.asyncio
async def test_delete(backend):
    """Deleted key is no longer retrievable."""
    await backend.store("key1", b"data", ttl=60)
    await backend.delete("key1")
    assert await backend.get("key1") is None


@pytest.mark.asyncio
async def test_overwrite(backend):
    """Storing to the same key overwrites the value."""
    await backend.store("key1", b"old", ttl=60)
    await backend.store("key1", b"new", ttl=60)
    assert await backend.get("key1") == b"new"


@pytest.mark.asyncio
async def test_large_binary_data(backend):
    """Large binary blobs round-trip correctly."""
    data = b"\x00\xff" * 500_000  # 1 MB
    await backend.store("big", data, ttl=60)
    assert await backend.get("big") == data


@pytest.mark.asyncio
async def test_delete_missing_key_no_error(backend):
    """Deleting a non-existent key does not raise."""
    await backend.delete("nonexistent")


@pytest.mark.asyncio
async def test_keys_are_isolated(backend):
    """Different keys don't interfere with each other."""
    await backend.store("a", b"1", ttl=60)
    await backend.store("b", b"2", ttl=60)
    assert await backend.get("a") == b"1"
    assert await backend.get("b") == b"2"
    await backend.delete("a")
    assert await backend.get("a") is None
    assert await backend.get("b") == b"2"


@pytest.mark.asyncio
async def test_inflight_raises(backend):
    """Inflight operations raise (must go through Redis)."""
    with pytest.raises(NotImplementedError):
        await backend.acquire_inflight("lock", ttl=30)
    with pytest.raises(NotImplementedError):
        await backend.release_inflight("lock")


@pytest.mark.asyncio
async def test_wait_for_key_raises(backend):
    """wait_for_key raises (must go through CompositeResultBackend)."""
    with pytest.raises(NotImplementedError):
        await backend.wait_for_key("k", timeout=1.0)


@pytest.mark.asyncio
async def test_notify_key_raises(backend):
    """notify_key raises (must go through CompositeResultBackend)."""
    with pytest.raises(NotImplementedError):
        await backend.notify_key("k")


# =============================================================================
# CompositeResultBackend
# =============================================================================


@pytest.fixture
def mock_redis():
    """AsyncMock-based RedisResultBackend."""
    mock = AsyncMock(spec=RedisResultBackend)
    mock.get.return_value = b"redis-data"
    mock.acquire_inflight.return_value = True
    return mock


@pytest.fixture
def mock_frames():
    """AsyncMock-based StorageResultBackend."""
    mock = AsyncMock(spec=StorageResultBackend)
    mock.get.return_value = b"frame-data"
    return mock


@pytest.fixture
def composite(mock_redis, mock_frames):
    """CompositeResultBackend with mocked sub-backends."""
    return CompositeResultBackend(redis=mock_redis, frames=mock_frames)


@pytest.mark.asyncio
async def test_composite_routes_frames_to_storage(composite, mock_frames):
    """Keys with :frames: route to storage backend."""
    key = "provider-result:room1:frames:source:abc123"

    await composite.store(key, b"data", ttl=300)
    mock_frames.store.assert_awaited_once_with(key, b"data", 300)

    await composite.get(key)
    mock_frames.get.assert_awaited_once_with(key)

    await composite.delete(key)
    mock_frames.delete.assert_awaited_once_with(key)


@pytest.mark.asyncio
async def test_composite_routes_filesystem_to_redis(composite, mock_redis):
    """Keys with :filesystem: route to Redis."""
    key = "provider-result:room1:filesystem:local:abc123"

    await composite.store(key, b"data", ttl=300)
    mock_redis.store.assert_awaited_once_with(key, b"data", 300)

    await composite.get(key)
    mock_redis.get.assert_awaited_once_with(key)

    await composite.delete(key)
    mock_redis.delete.assert_awaited_once_with(key)


@pytest.mark.asyncio
async def test_composite_inflight_always_uses_redis(composite, mock_redis):
    """Inflight locks always go to Redis, even for frame keys."""
    frame_key = "provider-inflight:room1:frames:source:abc123"

    await composite.acquire_inflight(frame_key, ttl=30)
    mock_redis.acquire_inflight.assert_awaited_once_with(frame_key, 30)

    await composite.release_inflight(frame_key)
    mock_redis.release_inflight.assert_awaited_once_with(frame_key)


@pytest.mark.asyncio
async def test_composite_notify_key_delegates_to_redis(composite, mock_redis):
    """notify_key always delegates to Redis backend."""
    await composite.notify_key("provider-result:room1:frames:source:abc123")
    mock_redis.notify_key.assert_awaited_once_with(
        "provider-result:room1:frames:source:abc123"
    )


# ---------------------------------------------------------------------------
# RedisResultBackend
# ---------------------------------------------------------------------------


@pytest_asyncio.fixture
async def redis_raw():
    """Raw Redis client (no decode_responses) for binary operations."""
    redis = Redis.from_url("redis://localhost")
    await redis.flushdb()
    yield redis
    await redis.flushdb()
    await redis.aclose()


@pytest_asyncio.fixture
async def redis_backend(redis_raw):
    return RedisResultBackend(redis_raw)


@pytest.mark.asyncio
async def test_redis_store_and_get(redis_backend):
    await redis_backend.store("test-key", b"hello", ttl=60)
    result = await redis_backend.get("test-key")
    assert result == b"hello"


@pytest.mark.asyncio
async def test_redis_get_missing_returns_none(redis_backend):
    result = await redis_backend.get("nonexistent")
    assert result is None


@pytest.mark.asyncio
async def test_redis_delete(redis_backend):
    await redis_backend.store("del-key", b"data", ttl=60)
    await redis_backend.delete("del-key")
    assert await redis_backend.get("del-key") is None


@pytest.mark.asyncio
async def test_redis_acquire_inflight(redis_backend):
    assert await redis_backend.acquire_inflight("lock", ttl=30) is True
    assert await redis_backend.acquire_inflight("lock", ttl=30) is False


@pytest.mark.asyncio
async def test_redis_release_inflight(redis_backend):
    await redis_backend.acquire_inflight("lock", ttl=30)
    await redis_backend.release_inflight("lock")
    assert await redis_backend.acquire_inflight("lock", ttl=30) is True


@pytest.mark.asyncio
async def test_redis_wait_for_key_cache_hit(redis_backend):
    """Returns immediately if data is already cached."""
    await redis_backend.store("k", b"data", ttl=60)
    result = await redis_backend.wait_for_key("k", timeout=1.0)
    assert result == b"data"


@pytest.mark.asyncio
async def test_redis_wait_for_key_timeout(redis_backend):
    """Returns None after timeout when no data arrives."""
    result = await redis_backend.wait_for_key("missing", timeout=0.1)
    assert result is None


@pytest.mark.asyncio
async def test_redis_wait_for_key_notify(redis_backend):
    """Wakes up when notify_key is called after store."""

    async def _upload_after_delay():
        await asyncio.sleep(0.05)
        await redis_backend.store("k", b"arrived", ttl=60)
        await redis_backend.notify_key("k")

    task = asyncio.create_task(_upload_after_delay())
    result = await redis_backend.wait_for_key("k", timeout=2.0)
    await task
    assert result == b"arrived"


@pytest.mark.asyncio
async def test_redis_wait_for_key_race_condition(redis_backend):
    """Data stored between first cache check and subscribe is still found."""
    await redis_backend.store("k", b"already-there", ttl=60)
    result = await redis_backend.wait_for_key("k", timeout=0.1)
    assert result == b"already-there"


@pytest.mark.asyncio
async def test_redis_notify_key_no_waiters(redis_backend):
    """notify_key does not raise when there are no subscribers."""
    await redis_backend.notify_key("no-one-listening")
