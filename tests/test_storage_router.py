"""Tests for StorageRouter (provider frame counts and write protection)."""

import uuid
from collections.abc import AsyncGenerator

import pytest
import pytest_asyncio
from conftest import make_raw_frame
from redis.asyncio import Redis as AsyncRedis

from zndraw.exceptions import ProblemException, RoomReadOnly
from zndraw.storage import AsebytesStorage
from zndraw.storage.router import StorageRouter


@pytest_asyncio.fixture
async def default_backend() -> AsyncGenerator[AsebytesStorage, None]:
    """Default AsebytesStorage backend."""
    s = AsebytesStorage("memory://")
    yield s
    await s.close()


@pytest_asyncio.fixture
async def redis() -> AsyncGenerator[AsyncRedis, None]:  # type: ignore[type-arg]
    """Fresh Redis, flushed before and after."""
    r: AsyncRedis = AsyncRedis()  # type: ignore[type-arg]
    await r.flushdb()
    yield r
    await r.flushdb()
    await r.aclose()


@pytest.fixture
def room() -> str:
    """Unique room ID per test to avoid memory:// backend pollution."""
    return uuid.uuid4().hex


@pytest_asyncio.fixture
async def router(
    default_backend: AsebytesStorage,
    redis: AsyncRedis,  # type: ignore[type-arg]
) -> AsyncGenerator[StorageRouter, None]:
    """StorageRouter with default backend."""
    r = StorageRouter(default=default_backend, redis=redis)
    yield r
    await r.close()


# --- Normal room routing ---


@pytest.mark.asyncio
async def test_router_get_delegates_to_default(
    router: StorageRouter, default_backend: AsebytesStorage, room: str
) -> None:
    """Normal room reads go to default backend."""
    await default_backend.extend(room, [make_raw_frame({"a": 1})])
    result = await router.get(room, 0)
    assert result == make_raw_frame({"a": 1})


@pytest.mark.asyncio
async def test_router_extend_delegates_to_default(
    router: StorageRouter, room: str
) -> None:
    """Normal room writes go to default backend."""
    count = await router.extend(room, [make_raw_frame({"a": 1})])
    assert count == 1
    assert await router.get_length(room) == 1


@pytest.mark.asyncio
async def test_router_has_mount_false(router: StorageRouter, room: str) -> None:
    """Normal room has no mount."""
    assert not await router.has_mount(room)


# --- Provider frame count ---


@pytest.mark.asyncio
async def test_set_frame_count_makes_room_mounted(
    router: StorageRouter, room: str
) -> None:
    """set_frame_count marks room as provider-backed."""
    await router.set_frame_count(room, 50)
    assert await router.has_mount(room)


@pytest.mark.asyncio
async def test_get_length_returns_provider_count(
    router: StorageRouter, room: str
) -> None:
    """get_length returns provider frame count from Redis."""
    await router.set_frame_count(room, 42)
    assert await router.get_length(room) == 42


@pytest.mark.asyncio
async def test_clear_frame_count_removes_mount(
    router: StorageRouter, room: str
) -> None:
    """clear_frame_count removes provider-backed status."""
    await router.set_frame_count(room, 50)
    await router.clear_frame_count(room)
    assert not await router.has_mount(room)
    assert await router.get_length(room) == 0


# --- Provider-backed write protection ---


@pytest.mark.asyncio
async def test_extend_raises_on_provider_room(router: StorageRouter, room: str) -> None:
    """extend on provider-backed room raises RoomReadOnly."""
    await router.set_frame_count(room, 10)
    with pytest.raises(ProblemException) as exc_info:
        await router.extend(room, [])
    assert exc_info.value.problem.status == RoomReadOnly.status


@pytest.mark.asyncio
async def test_set_item_raises_on_provider_room(
    router: StorageRouter, room: str
) -> None:
    """set_item on provider-backed room raises RoomReadOnly."""
    await router.set_frame_count(room, 10)
    with pytest.raises(ProblemException) as exc_info:
        await router.set_item(room, 0, {})
    assert exc_info.value.problem.status == RoomReadOnly.status


@pytest.mark.asyncio
async def test_delete_range_raises_on_provider_room(
    router: StorageRouter, room: str
) -> None:
    """delete_range on provider-backed room raises RoomReadOnly."""
    await router.set_frame_count(room, 10)
    with pytest.raises(ProblemException) as exc_info:
        await router.delete_range(room, 0, 1)
    assert exc_info.value.problem.status == RoomReadOnly.status


# --- Clear cleans up frame count ---


@pytest.mark.asyncio
async def test_clear_removes_provider_count(router: StorageRouter, room: str) -> None:
    """clear() removes provider frame count from Redis."""
    await router.set_frame_count(room, 10)
    await router.clear(room)
    assert not await router.has_mount(room)
    assert await router.get_length(room) == 0


# --- Provider-backed room reads return positional None for virtual frames ---


@pytest.mark.asyncio
async def test_get_range_returns_none_for_provider_frames(
    router: StorageRouter, room: str
) -> None:
    """get_range on provider-backed room returns positional None list.

    Default storage is empty but Redis reports 100 frames — get_range
    must return [None] * count, not [].
    """
    await router.set_frame_count(room, 100)
    result = await router.get_range(room, 0, 5)
    assert len(result) == 5
    assert all(f is None for f in result)


@pytest.mark.asyncio
async def test_get_range_stop_none_returns_all_provider_frames(
    router: StorageRouter, room: str
) -> None:
    """get_range with stop=None returns all provider frames as None."""
    await router.set_frame_count(room, 10)
    result = await router.get_range(room, 0, None)
    assert len(result) == 10
    assert all(f is None for f in result)


@pytest.mark.asyncio
async def test_get_many_returns_none_for_provider_frames(
    router: StorageRouter, room: str
) -> None:
    """get_many on provider-backed room returns positional None list."""
    await router.set_frame_count(room, 100)
    result = await router.get_many(room, [0, 5, 99])
    assert len(result) == 3
    assert all(f is None for f in result)


@pytest.mark.asyncio
async def test_get_range_clamps_to_provider_count(
    router: StorageRouter, room: str
) -> None:
    """get_range with stop beyond provider count is clamped."""
    await router.set_frame_count(room, 10)
    result = await router.get_range(room, 8, 20)
    assert len(result) == 2
    assert all(f is None for f in result)
