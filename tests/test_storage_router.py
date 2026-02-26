"""Tests for StorageRouter (provider frame counts and write protection)."""

from collections.abc import AsyncGenerator

import pytest
import pytest_asyncio
from conftest import make_raw_frame
from redis.asyncio import Redis as AsyncRedis

from zndraw.exceptions import ProblemException, RoomReadOnly
from zndraw.storage import InMemoryStorage
from zndraw.storage.router import StorageRouter

ROOM_NORMAL = "normal-room"
ROOM_PROVIDER = "provider-room"


@pytest_asyncio.fixture
async def default_backend() -> AsyncGenerator[InMemoryStorage, None]:
    """Default InMemoryStorage backend."""
    s = InMemoryStorage()
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


@pytest_asyncio.fixture
async def router(
    default_backend: InMemoryStorage,
    redis: AsyncRedis,  # type: ignore[type-arg]
) -> AsyncGenerator[StorageRouter, None]:
    """StorageRouter with default backend."""
    r = StorageRouter(default=default_backend, redis=redis)
    yield r
    await r.close()


# --- Normal room routing ---


@pytest.mark.asyncio
async def test_router_get_delegates_to_default(
    router: StorageRouter, default_backend: InMemoryStorage
) -> None:
    """Normal room reads go to default backend."""
    await default_backend.extend(ROOM_NORMAL, [make_raw_frame({"a": 1})])
    result = await router.get(ROOM_NORMAL, 0)
    assert result == make_raw_frame({"a": 1})


@pytest.mark.asyncio
async def test_router_extend_delegates_to_default(router: StorageRouter) -> None:
    """Normal room writes go to default backend."""
    count = await router.extend(ROOM_NORMAL, [make_raw_frame({"a": 1})])
    assert count == 1
    assert await router.get_length(ROOM_NORMAL) == 1


@pytest.mark.asyncio
async def test_router_has_mount_false(router: StorageRouter) -> None:
    """Normal room has no mount."""
    assert not await router.has_mount(ROOM_NORMAL)


# --- Provider frame count ---


@pytest.mark.asyncio
async def test_set_frame_count_makes_room_mounted(router: StorageRouter) -> None:
    """set_frame_count marks room as provider-backed."""
    await router.set_frame_count(ROOM_PROVIDER, 50)
    assert await router.has_mount(ROOM_PROVIDER)


@pytest.mark.asyncio
async def test_get_length_returns_provider_count(router: StorageRouter) -> None:
    """get_length returns provider frame count from Redis."""
    await router.set_frame_count(ROOM_PROVIDER, 42)
    assert await router.get_length(ROOM_PROVIDER) == 42


@pytest.mark.asyncio
async def test_clear_frame_count_removes_mount(router: StorageRouter) -> None:
    """clear_frame_count removes provider-backed status."""
    await router.set_frame_count(ROOM_PROVIDER, 50)
    await router.clear_frame_count(ROOM_PROVIDER)
    assert not await router.has_mount(ROOM_PROVIDER)
    assert await router.get_length(ROOM_PROVIDER) == 0


# --- Provider-backed write protection ---


@pytest.mark.asyncio
async def test_extend_raises_on_provider_room(router: StorageRouter) -> None:
    """extend on provider-backed room raises RoomReadOnly."""
    await router.set_frame_count(ROOM_PROVIDER, 10)
    with pytest.raises(ProblemException) as exc_info:
        await router.extend(ROOM_PROVIDER, [])
    assert exc_info.value.problem.status == RoomReadOnly.status


@pytest.mark.asyncio
async def test_set_item_raises_on_provider_room(router: StorageRouter) -> None:
    """set_item on provider-backed room raises RoomReadOnly."""
    await router.set_frame_count(ROOM_PROVIDER, 10)
    with pytest.raises(ProblemException) as exc_info:
        await router.set_item(ROOM_PROVIDER, 0, {})
    assert exc_info.value.problem.status == RoomReadOnly.status


@pytest.mark.asyncio
async def test_delete_range_raises_on_provider_room(router: StorageRouter) -> None:
    """delete_range on provider-backed room raises RoomReadOnly."""
    await router.set_frame_count(ROOM_PROVIDER, 10)
    with pytest.raises(ProblemException) as exc_info:
        await router.delete_range(ROOM_PROVIDER, 0, 1)
    assert exc_info.value.problem.status == RoomReadOnly.status


# --- Clear cleans up frame count ---


@pytest.mark.asyncio
async def test_clear_removes_provider_count(router: StorageRouter) -> None:
    """clear() removes provider frame count from Redis."""
    await router.set_frame_count(ROOM_PROVIDER, 10)
    await router.clear(ROOM_PROVIDER)
    assert not await router.has_mount(ROOM_PROVIDER)
    assert await router.get_length(ROOM_PROVIDER) == 0


# --- Provider-backed room reads return positional None for virtual frames ---


@pytest.mark.asyncio
async def test_get_range_returns_none_for_provider_frames(
    router: StorageRouter,
) -> None:
    """get_range on provider-backed room returns positional None list.

    Default storage is empty but Redis reports 100 frames — get_range
    must return [None] * count, not [].
    """
    await router.set_frame_count(ROOM_PROVIDER, 100)
    result = await router.get_range(ROOM_PROVIDER, 0, 5)
    assert len(result) == 5
    assert all(f is None for f in result)


@pytest.mark.asyncio
async def test_get_range_stop_none_returns_all_provider_frames(
    router: StorageRouter,
) -> None:
    """get_range with stop=None returns all provider frames as None."""
    await router.set_frame_count(ROOM_PROVIDER, 10)
    result = await router.get_range(ROOM_PROVIDER, 0, None)
    assert len(result) == 10
    assert all(f is None for f in result)


@pytest.mark.asyncio
async def test_get_many_returns_none_for_provider_frames(
    router: StorageRouter,
) -> None:
    """get_many on provider-backed room returns positional None list."""
    await router.set_frame_count(ROOM_PROVIDER, 100)
    result = await router.get_many(ROOM_PROVIDER, [0, 5, 99])
    assert len(result) == 3
    assert all(f is None for f in result)


@pytest.mark.asyncio
async def test_get_range_clamps_to_provider_count(
    router: StorageRouter,
) -> None:
    """get_range with stop beyond provider count is clamped."""
    await router.set_frame_count(ROOM_PROVIDER, 10)
    result = await router.get_range(ROOM_PROVIDER, 8, 20)
    assert len(result) == 2
    assert all(f is None for f in result)
