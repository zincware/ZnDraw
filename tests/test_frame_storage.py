"""Tests for FrameStorage with AsyncBlobIO subscript API."""

import uuid

import pytest
from conftest import make_raw_frame

from zndraw.storage import FrameStorage


@pytest.fixture
def room():
    """Unique room ID per test to avoid memory:// backend pollution."""
    return uuid.uuid4().hex


@pytest.mark.asyncio
async def test_empty_room_has_zero_length(
    frame_storage: FrameStorage, room: str
) -> None:
    assert await frame_storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_get_raises_for_empty_room(
    frame_storage: FrameStorage, room: str
) -> None:
    with pytest.raises(IndexError):
        await frame_storage[room][0]


@pytest.mark.asyncio
async def test_extend_and_get(frame_storage: FrameStorage, room: str) -> None:
    frame = make_raw_frame({"a": 1})
    count = await frame_storage[room].extend([frame])
    assert count == 1
    assert await frame_storage.get_length(room) == 1
    result = await frame_storage[room][0]
    assert result == frame


@pytest.mark.asyncio
async def test_extend_multiple(frame_storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    count = await frame_storage[room].extend(frames)
    assert count == 5
    assert await frame_storage.get_length(room) == 5


@pytest.mark.asyncio
async def test_get_range(frame_storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await frame_storage[room].extend(frames)
    result = await frame_storage[room][1:4].to_list()
    assert len(result) == 3
    assert result[0] == frames[1]
    assert result[2] == frames[3]


@pytest.mark.asyncio
async def test_get_many(frame_storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await frame_storage[room].extend(frames)
    result = await frame_storage[room][[0, 2, 4]].to_list()
    assert len(result) == 3
    assert result[0] == frames[0]
    assert result[1] == frames[2]
    assert result[2] == frames[4]


@pytest.mark.asyncio
async def test_get_many_oob_raises(frame_storage: FrameStorage, room: str) -> None:
    await frame_storage[room].extend([make_raw_frame({"a": 1})])
    with pytest.raises(IndexError):
        await frame_storage[room][[0, 99]].to_list()


@pytest.mark.asyncio
async def test_set_item(frame_storage: FrameStorage, room: str) -> None:
    await frame_storage[room].extend([make_raw_frame({"a": 1})])
    new_frame = make_raw_frame({"a": 2})
    await frame_storage[room][0].set(new_frame)
    result = await frame_storage[room][0]
    assert result == new_frame


@pytest.mark.asyncio
async def test_merge_item(frame_storage: FrameStorage, room: str) -> None:
    original = make_raw_frame({"a": 1, "b": 2})
    await frame_storage[room].extend([original])
    partial = make_raw_frame({"b": 99})
    await frame_storage[room][0].update(partial)
    result = await frame_storage[room][0]
    assert result is not None
    assert result[b"b"] == make_raw_frame({"b": 99})[b"b"]
    assert result[b"a"] == original[b"a"]


@pytest.mark.asyncio
async def test_delete_range(frame_storage: FrameStorage, room: str) -> None:
    frames = [make_raw_frame({"i": i}) for i in range(5)]
    await frame_storage[room].extend(frames)
    await frame_storage[room][1:3].delete()
    assert await frame_storage.get_length(room) == 3


@pytest.mark.asyncio
async def test_clear(frame_storage: FrameStorage, room: str) -> None:
    await frame_storage[room].extend([make_raw_frame({"a": 1})])
    await frame_storage.clear(room)
    assert await frame_storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_room_isolation(frame_storage: FrameStorage) -> None:
    room_a = uuid.uuid4().hex
    room_b = uuid.uuid4().hex
    await frame_storage[room_a].extend([make_raw_frame({"a": 1})])
    await frame_storage[room_b].extend([make_raw_frame({"b": 2})])
    assert await frame_storage.get_length(room_a) == 1
    assert await frame_storage.get_length(room_b) == 1
    assert await frame_storage[room_a][0] == make_raw_frame({"a": 1})
    assert await frame_storage[room_b][0] == make_raw_frame({"b": 2})
