"""Tests for FrameStorage provider metadata (frame count, mount, room isolation)."""

import uuid

import pytest

from zndraw.storage import FrameStorage


@pytest.fixture
def room() -> str:
    """Unique room ID per test."""
    return uuid.uuid4().hex


@pytest.mark.asyncio
async def test_has_mount_false_by_default(
    frame_storage: FrameStorage, room: str
) -> None:
    """Normal room has no mount."""
    assert not await frame_storage.has_mount(room)


@pytest.mark.asyncio
async def test_set_frame_count_makes_room_mounted(
    frame_storage: FrameStorage, room: str
) -> None:
    """set_frame_count marks room as provider-backed."""
    await frame_storage.set_frame_count(room, 50)
    assert await frame_storage.has_mount(room)


@pytest.mark.asyncio
async def test_get_length_returns_provider_count(
    frame_storage: FrameStorage, room: str
) -> None:
    """get_length returns provider frame count from Redis when storage is empty."""
    await frame_storage.set_frame_count(room, 42)
    assert await frame_storage.get_length(room) == 42


@pytest.mark.asyncio
async def test_get_length_returns_storage_count(
    frame_storage: FrameStorage, room: str
) -> None:
    """get_length returns storage count when storage has frames, ignores Redis."""
    await frame_storage.set_frame_count(room, 100)
    await frame_storage[room].extend([{b"a": b"1"}, {b"b": b"2"}])
    assert await frame_storage.get_length(room) == 2


@pytest.mark.asyncio
async def test_clear_frame_count_removes_mount(
    frame_storage: FrameStorage, room: str
) -> None:
    """clear_frame_count removes provider-backed status."""
    await frame_storage.set_frame_count(room, 50)
    await frame_storage.clear_frame_count(room)
    assert not await frame_storage.has_mount(room)
    assert await frame_storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_clear_removes_provider_count(
    frame_storage: FrameStorage, room: str
) -> None:
    """clear() removes both io data and provider frame count from Redis."""
    await frame_storage.set_frame_count(room, 10)
    await frame_storage.clear(room)
    assert not await frame_storage.has_mount(room)
    assert await frame_storage.get_length(room) == 0


@pytest.mark.asyncio
async def test_room_isolation(frame_storage: FrameStorage) -> None:
    """Different rooms have independent AsyncBlobIO instances."""
    room_a = uuid.uuid4().hex
    room_b = uuid.uuid4().hex

    await frame_storage.set_frame_count(room_a, 10)
    await frame_storage[room_b].extend([{b"x": b"1"}])

    assert await frame_storage.has_mount(room_a)
    assert not await frame_storage.has_mount(room_b)
    assert await frame_storage.get_length(room_a) == 10
    assert await frame_storage.get_length(room_b) == 1
