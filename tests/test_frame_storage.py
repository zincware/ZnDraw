"""Tests for FrameStorage with AsyncBlobIO subscript API."""

import uuid
from unittest.mock import AsyncMock

import pytest
import pytest_asyncio
from conftest import make_raw_frame

from zndraw.storage import FrameStorage


def _make_mock_redis() -> AsyncMock:
    """Create a mock Redis client for tests that don't need real Redis."""
    mock = AsyncMock()
    mock.get = AsyncMock(return_value=None)
    mock.exists = AsyncMock(return_value=0)
    mock.set = AsyncMock()
    mock.delete = AsyncMock()
    return mock


@pytest_asyncio.fixture
async def storage():
    """Fresh in-memory FrameStorage with mock Redis."""
    s = FrameStorage("memory://", _make_mock_redis())
    yield s
    await s.close()


@pytest.fixture
def room():
    """Unique room ID per test to avoid memory:// backend pollution."""
    return uuid.uuid4().hex


@pytest.mark.asyncio
async def test_empty_room_has_zero_length(storage: FrameStorage, room: str) -> None:
    assert await storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_get_raises_for_empty_room(storage: FrameStorage, room: str) -> None:
    with pytest.raises(IndexError):
        await storage[room][0]


@pytest.mark.asyncio
async def test_extend_and_get(storage: FrameStorage, room: str) -> None:
    frame = make_raw_frame({"a": 1})
    count = await storage[room].extend([frame])
    assert count == 1
    assert await storage.get_length(room) == 1
    result = await storage[room][0]
    assert result == frame


@pytest.mark.asyncio
async def test_extend_multiple(storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    count = await storage[room].extend(frames)
    assert count == 5
    assert await storage.get_length(room) == 5


@pytest.mark.asyncio
async def test_get_range(storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage[room].extend(frames)
    result = await storage[room][1:4].to_list()
    assert len(result) == 3
    assert result[0] == frames[1]
    assert result[2] == frames[3]


@pytest.mark.asyncio
async def test_get_many(storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage[room].extend(frames)
    result = await storage[room][[0, 2, 4]].to_list()
    assert len(result) == 3
    assert result[0] == frames[0]
    assert result[1] == frames[2]
    assert result[2] == frames[4]


@pytest.mark.asyncio
async def test_get_many_oob_raises(storage: FrameStorage, room: str) -> None:
    await storage[room].extend([make_raw_frame({"a": 1})])
    with pytest.raises(IndexError):
        await storage[room][[0, 99]].to_list()


@pytest.mark.asyncio
async def test_set_item(storage: FrameStorage, room: str) -> None:
    await storage[room].extend([make_raw_frame({"a": 1})])
    new_frame = make_raw_frame({"a": 2})
    await storage[room][0].set(new_frame)
    result = await storage[room][0]
    assert result == new_frame


@pytest.mark.asyncio
async def test_merge_item(storage: FrameStorage, room: str) -> None:
    original = make_raw_frame({"a": 1, "b": 2})
    await storage[room].extend([original])
    partial = make_raw_frame({"b": 99})
    await storage[room][0].update(partial)
    result = await storage[room][0]
    assert result is not None
    assert result[b"b"] == make_raw_frame({"b": 99})[b"b"]
    assert result[b"a"] == original[b"a"]


@pytest.mark.asyncio
async def test_delete_range(storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage[room].extend(frames)
    await storage[room][1:3].delete()
    assert await storage.get_length(room) == 3


@pytest.mark.asyncio
async def test_clear(storage: FrameStorage, room: str) -> None:
    await storage[room].extend([make_raw_frame({"a": 1})])
    await storage.clear(room)
    assert await storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_room_isolation(storage: FrameStorage) -> None:
    room_a = uuid.uuid4().hex
    room_b = uuid.uuid4().hex
    await storage[room_a].extend([make_raw_frame({"a": 1})])
    await storage[room_b].extend([make_raw_frame({"b": 2})])
    assert await storage.get_length(room_a) == 1
    assert await storage.get_length(room_b) == 1
    assert await storage[room_a][0] == make_raw_frame({"a": 1})
    assert await storage[room_b][0] == make_raw_frame({"b": 2})
