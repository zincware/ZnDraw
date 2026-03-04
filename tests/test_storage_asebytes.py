"""Tests for AsebytesStorage wrapper."""

import uuid

import pytest
import pytest_asyncio
from conftest import make_raw_frame

from zndraw.storage import AsebytesStorage


@pytest_asyncio.fixture
async def storage():
    """Fresh in-memory AsebytesStorage."""
    s = AsebytesStorage("memory://")
    yield s
    await s.close()


@pytest.fixture
def room():
    """Unique room ID per test to avoid memory:// backend pollution."""
    return uuid.uuid4().hex


@pytest.mark.asyncio
async def test_empty_room_has_zero_length(storage: AsebytesStorage, room: str) -> None:
    assert await storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_get_returns_none_for_empty_room(
    storage: AsebytesStorage, room: str
) -> None:
    assert await storage.get(room, 0) is None


@pytest.mark.asyncio
async def test_extend_and_get(storage: AsebytesStorage, room: str) -> None:
    frame = make_raw_frame({"a": 1})
    count = await storage.extend(room, [frame])
    assert count == 1
    assert await storage.get_length(room) == 1
    result = await storage.get(room, 0)
    assert result == frame


@pytest.mark.asyncio
async def test_extend_multiple(storage: AsebytesStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    count = await storage.extend(room, frames)
    assert count == 5
    assert await storage.get_length(room) == 5


@pytest.mark.asyncio
async def test_get_range(storage: AsebytesStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage.extend(room, frames)
    result = await storage.get_range(room, 1, 4)
    assert len(result) == 3
    assert result[0] == frames[1]
    assert result[2] == frames[3]


@pytest.mark.asyncio
async def test_get_many(storage: AsebytesStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage.extend(room, frames)
    result = await storage.get_many(room, [0, 2, 4])
    assert len(result) == 3
    assert result[0] == frames[0]
    assert result[1] == frames[2]
    assert result[2] == frames[4]


@pytest.mark.asyncio
async def test_get_many_oob_returns_none(storage: AsebytesStorage, room: str) -> None:
    await storage.extend(room, [make_raw_frame({"a": 1})])
    result = await storage.get_many(room, [0, 99])
    assert result[0] is not None
    assert result[1] is None


@pytest.mark.asyncio
async def test_set_item(storage: AsebytesStorage, room: str) -> None:
    await storage.extend(room, [make_raw_frame({"a": 1})])
    new_frame = make_raw_frame({"a": 2})
    await storage.set_item(room, 0, new_frame)
    result = await storage.get(room, 0)
    assert result == new_frame


@pytest.mark.asyncio
async def test_merge_item(storage: AsebytesStorage, room: str) -> None:
    original = make_raw_frame({"a": 1, "b": 2})
    await storage.extend(room, [original])
    partial = make_raw_frame({"b": 99})
    await storage.merge_item(room, 0, partial)
    result = await storage.get(room, 0)
    assert result is not None
    assert result[b"b"] == make_raw_frame({"b": 99})[b"b"]
    assert result[b"a"] == original[b"a"]


@pytest.mark.asyncio
async def test_delete_range(storage: AsebytesStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage.extend(room, frames)
    await storage.delete_range(room, 1, 3)
    assert await storage.get_length(room) == 3


@pytest.mark.asyncio
async def test_clear(storage: AsebytesStorage, room: str) -> None:
    await storage.extend(room, [make_raw_frame({"a": 1})])
    await storage.clear(room)
    assert await storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_reserve(storage: AsebytesStorage, room: str) -> None:
    await storage.reserve(room, 3)
    assert await storage.get_length(room) == 3
    result = await storage.get(room, 0)
    assert result is None


@pytest.mark.asyncio
async def test_remove_items(storage: AsebytesStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(3)]
    await storage.extend(room, frames)
    await storage.remove_items(room, [1])
    assert await storage.get_length(room) == 3
    assert await storage.get(room, 1) is None
    assert await storage.get(room, 0) == frames[0]
    assert await storage.get(room, 2) == frames[2]


@pytest.mark.asyncio
async def test_room_isolation(storage: AsebytesStorage) -> None:
    room_a = uuid.uuid4().hex
    room_b = uuid.uuid4().hex
    await storage.extend(room_a, [make_raw_frame({"a": 1})])
    await storage.extend(room_b, [make_raw_frame({"b": 2})])
    assert await storage.get_length(room_a) == 1
    assert await storage.get_length(room_b) == 1
    assert await storage.get(room_a, 0) == make_raw_frame({"a": 1})
    assert await storage.get(room_b, 0) == make_raw_frame({"b": 2})


@pytest.mark.asyncio
async def test_close_clears_rooms(storage: AsebytesStorage, room: str) -> None:
    await storage.extend(room, [make_raw_frame({"a": 1})])
    await storage.close()
    # After close, room data is cleared
    assert await storage.get_length(room) == 0
