"""Tests for FrameStorage provider metadata (frame count, mount, room isolation)."""

import uuid
from collections.abc import AsyncGenerator

import pytest
import pytest_asyncio
from redis.asyncio import Redis as AsyncRedis

from zndraw.storage import FrameStorage


@pytest_asyncio.fixture
async def redis() -> AsyncGenerator[AsyncRedis, None]:  # type: ignore[type-arg]
    """Fresh Redis, flushed before and after."""
    r: AsyncRedis = AsyncRedis.from_url("redis://localhost:6379/1")  # type: ignore[type-arg]
    await r.flushdb()
    yield r
    await r.flushdb()
    await r.aclose()


@pytest_asyncio.fixture
async def storage(redis: AsyncRedis) -> AsyncGenerator[FrameStorage, None]:  # type: ignore[type-arg]
    """FrameStorage with real Redis backend."""
    s = FrameStorage("memory://", redis)
    yield s
    await s.close()


@pytest.fixture
def room() -> str:
    """Unique room ID per test."""
    return uuid.uuid4().hex


@pytest.mark.asyncio
async def test_has_mount_false_by_default(storage: FrameStorage, room: str) -> None:
    """Normal room has no mount."""
    assert not await storage.has_mount(room)


@pytest.mark.asyncio
async def test_set_frame_count_makes_room_mounted(
    storage: FrameStorage, room: str
) -> None:
    """set_frame_count marks room as provider-backed."""
    await storage.set_frame_count(room, 50)
    assert await storage.has_mount(room)


@pytest.mark.asyncio
async def test_get_length_returns_provider_count(
    storage: FrameStorage, room: str
) -> None:
    """get_length returns provider frame count from Redis when storage is empty."""
    await storage.set_frame_count(room, 42)
    assert await storage.get_length(room) == 42


@pytest.mark.asyncio
async def test_get_length_returns_storage_count(
    storage: FrameStorage, room: str
) -> None:
    """get_length returns storage count when storage has frames, ignores Redis."""
    await storage.set_frame_count(room, 100)
    await storage[room].extend([{b"a": b"1"}, {b"b": b"2"}])
    assert await storage.get_length(room) == 2


@pytest.mark.asyncio
async def test_clear_frame_count_removes_mount(
    storage: FrameStorage, room: str
) -> None:
    """clear_frame_count removes provider-backed status."""
    await storage.set_frame_count(room, 50)
    await storage.clear_frame_count(room)
    assert not await storage.has_mount(room)
    assert await storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_clear_removes_provider_count(storage: FrameStorage, room: str) -> None:
    """clear() removes both io data and provider frame count from Redis."""
    await storage.set_frame_count(room, 10)
    await storage.clear(room)
    assert not await storage.has_mount(room)
    assert await storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_room_isolation(storage: FrameStorage) -> None:
    """Different rooms have independent AsyncBlobIO instances."""
    room_a = uuid.uuid4().hex
    room_b = uuid.uuid4().hex

    await storage.set_frame_count(room_a, 10)
    await storage[room_b].extend([{b"x": b"1"}])

    assert await storage.has_mount(room_a)
    assert not await storage.has_mount(room_b)
    assert await storage.get_length(room_a) == 10
    assert await storage.get_length(room_b) == 1
