"""Tests for AsebytesStorage wrapper."""

import pytest
import pytest_asyncio
from conftest import make_raw_frame

from zndraw.storage import AsebytesStorage

ROOM = "test-room"


@pytest_asyncio.fixture
async def storage():
    """Fresh in-memory AsebytesStorage."""
    s = AsebytesStorage("memory://")
    yield s
    await s.close()


@pytest.mark.asyncio
async def test_empty_room_has_zero_length(storage: AsebytesStorage) -> None:
    assert await storage.get_length(ROOM) == 0


@pytest.mark.asyncio
async def test_get_returns_none_for_empty_room(storage: AsebytesStorage) -> None:
    assert await storage.get(ROOM, 0) is None


@pytest.mark.asyncio
async def test_extend_and_get(storage: AsebytesStorage) -> None:
    frame = make_raw_frame({"a": 1})
    count = await storage.extend(ROOM, [frame])
    assert count == 1
    assert await storage.get_length(ROOM) == 1
    result = await storage.get(ROOM, 0)
    assert result == frame


@pytest.mark.asyncio
async def test_extend_multiple(storage: AsebytesStorage) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    count = await storage.extend(ROOM, frames)
    assert count == 5
    assert await storage.get_length(ROOM) == 5


@pytest.mark.asyncio
async def test_get_range(storage: AsebytesStorage) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage.extend(ROOM, frames)
    result = await storage.get_range(ROOM, 1, 4)
    assert len(result) == 3
    assert result[0] == frames[1]
    assert result[2] == frames[3]


@pytest.mark.asyncio
async def test_get_many(storage: AsebytesStorage) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage.extend(ROOM, frames)
    result = await storage.get_many(ROOM, [0, 2, 4])
    assert len(result) == 3
    assert result[0] == frames[0]
    assert result[1] == frames[2]
    assert result[2] == frames[4]


@pytest.mark.asyncio
async def test_get_many_oob_returns_none(storage: AsebytesStorage) -> None:
    await storage.extend(ROOM, [make_raw_frame({"a": 1})])
    result = await storage.get_many(ROOM, [0, 99])
    assert result[0] is not None
    assert result[1] is None


@pytest.mark.asyncio
async def test_set_item(storage: AsebytesStorage) -> None:
    await storage.extend(ROOM, [make_raw_frame({"a": 1})])
    new_frame = make_raw_frame({"a": 2})
    await storage.set_item(ROOM, 0, new_frame)
    result = await storage.get(ROOM, 0)
    assert result == new_frame


@pytest.mark.asyncio
async def test_merge_item(storage: AsebytesStorage) -> None:
    original = make_raw_frame({"a": 1, "b": 2})
    await storage.extend(ROOM, [original])
    partial = make_raw_frame({"b": 99})
    await storage.merge_item(ROOM, 0, partial)
    result = await storage.get(ROOM, 0)
    assert result is not None
    assert result[b"b"] == make_raw_frame({"b": 99})[b"b"]
    assert result[b"a"] == original[b"a"]


@pytest.mark.asyncio
async def test_delete_range(storage: AsebytesStorage) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await storage.extend(ROOM, frames)
    await storage.delete_range(ROOM, 1, 3)
    assert await storage.get_length(ROOM) == 3


@pytest.mark.asyncio
async def test_clear(storage: AsebytesStorage) -> None:
    await storage.extend(ROOM, [make_raw_frame({"a": 1})])
    await storage.clear(ROOM)
    assert await storage.get_length(ROOM) == 0


@pytest.mark.asyncio
async def test_reserve(storage: AsebytesStorage) -> None:
    await storage.reserve(ROOM, 3)
    assert await storage.get_length(ROOM) == 3
    result = await storage.get(ROOM, 0)
    assert result is None


@pytest.mark.asyncio
async def test_remove_items(storage: AsebytesStorage) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(3)]
    await storage.extend(ROOM, frames)
    await storage.remove_items(ROOM, [1])
    assert await storage.get_length(ROOM) == 3
    assert await storage.get(ROOM, 1) is None
    assert await storage.get(ROOM, 0) == frames[0]
    assert await storage.get(ROOM, 2) == frames[2]


@pytest.mark.asyncio
async def test_room_isolation(storage: AsebytesStorage) -> None:
    await storage.extend("room-a", [make_raw_frame({"a": 1})])
    await storage.extend("room-b", [make_raw_frame({"b": 2})])
    assert await storage.get_length("room-a") == 1
    assert await storage.get_length("room-b") == 1
    assert await storage.get("room-a", 0) == make_raw_frame({"a": 1})
    assert await storage.get("room-b", 0) == make_raw_frame({"b": 2})


@pytest.mark.asyncio
async def test_close_clears_rooms(storage: AsebytesStorage) -> None:
    await storage.extend(ROOM, [make_raw_frame({"a": 1})])
    await storage.close()
    # After close, rooms dict is empty — new access creates fresh backend
    assert await storage.get_length(ROOM) == 0
