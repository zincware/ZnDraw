"""Tests for InMemoryStorage backend."""

from collections.abc import AsyncGenerator

import msgpack
import pytest
import pytest_asyncio
from conftest import make_raw_frame

from zndraw.storage import InMemoryStorage, StorageBackend
from zndraw.storage.base import RawFrame


@pytest_asyncio.fixture
async def storage() -> AsyncGenerator[InMemoryStorage, None]:
    """Create a fresh InMemoryStorage instance for each test."""
    s = InMemoryStorage()
    yield s
    await s.close()


# --- Interface ---


def test_memory_is_subclass_of_storage_backend() -> None:
    """InMemoryStorage is a proper subclass of StorageBackend."""
    assert issubclass(InMemoryStorage, StorageBackend)


def test_memory_can_instantiate() -> None:
    """InMemoryStorage can be instantiated."""
    storage = InMemoryStorage()
    assert isinstance(storage, StorageBackend)


# --- Basics ---


@pytest.mark.asyncio
async def test_memory_empty_room_has_zero_length(storage: InMemoryStorage) -> None:
    """New room has no frames."""
    length = await storage.get_length(room_id="room1")
    assert length == 0


@pytest.mark.asyncio
async def test_memory_get_nonexistent_returns_none(storage: InMemoryStorage) -> None:
    """Get on empty room returns None."""
    result = await storage.get(room_id="room1", index=0)
    assert result is None


@pytest.mark.asyncio
async def test_memory_get_out_of_bounds_returns_none(storage: InMemoryStorage) -> None:
    """Get with out of bounds index returns None."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])
    result = await storage.get(room_id="room1", index=5)
    assert result is None


# --- Extend and Get ---


@pytest.mark.asyncio
async def test_memory_extend_and_get(storage: InMemoryStorage) -> None:
    """Extend then retrieve frames."""
    frames = [{"a": 1}, {"b": 2}, {"c": 3}]
    new_length = await storage.extend(room_id="room1", frames=frames)

    assert new_length == 3
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) == make_raw_frame({"b": 2})
    assert await storage.get(room_id="room1", index=2) == make_raw_frame({"c": 3})


@pytest.mark.asyncio
async def test_memory_extend_multiple_times(storage: InMemoryStorage) -> None:
    """Extend can be called multiple times."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])
    new_length = await storage.extend(room_id="room1", frames=[{"b": 2}, {"c": 3}])

    assert new_length == 3
    assert await storage.get_length(room_id="room1") == 3


@pytest.mark.asyncio
async def test_memory_get_range(storage: InMemoryStorage) -> None:
    """Range slicing works."""
    frames = [{"a": 1}, {"b": 2}, {"c": 3}, {"d": 4}]
    await storage.extend(room_id="room1", frames=frames)

    result = await storage.get_range(room_id="room1", start=1, stop=3)
    assert result == [make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]


@pytest.mark.asyncio
async def test_memory_get_range_empty(storage: InMemoryStorage) -> None:
    """Range on empty room returns empty list."""
    result = await storage.get_range(room_id="room1", start=0, stop=10)
    assert result == []


@pytest.mark.asyncio
async def test_memory_get_range_partial(storage: InMemoryStorage) -> None:
    """Range beyond bounds returns available frames only."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    result = await storage.get_range(room_id="room1", start=0, stop=10)
    assert result == [make_raw_frame({"a": 1}), make_raw_frame({"b": 2})]


# --- SetItem ---


@pytest.mark.asyncio
async def test_memory_set_item(storage: InMemoryStorage) -> None:
    """Can overwrite a frame."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    await storage.set_item(room_id="room1", index=1, frame={"updated": True})

    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) == make_raw_frame(
        {"updated": True}
    )


@pytest.mark.asyncio
async def test_memory_set_item_out_of_bounds_raises(storage: InMemoryStorage) -> None:
    """Setting frame at invalid index raises IndexError."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])

    with pytest.raises(IndexError):
        await storage.set_item(room_id="room1", index=5, frame={"bad": True})


@pytest.mark.asyncio
async def test_memory_set_item_on_empty_raises(storage: InMemoryStorage) -> None:
    """Setting frame on empty room raises IndexError."""
    with pytest.raises(IndexError):
        await storage.set_item(room_id="room1", index=0, frame={"bad": True})


# --- MergeItem ---


@pytest.mark.asyncio
async def test_memory_merge_item_updates_existing_keys(
    storage: InMemoryStorage,
) -> None:
    """Merge overwrites existing keys."""
    await storage.extend(room_id="room1", frames=[{"a": 1, "b": 2}])
    partial: RawFrame = {b"a": msgpack.packb(99)}  # type: ignore[reportAssignmentType]
    await storage.merge_item(room_id="room1", index=0, partial=partial)

    result = await storage.get(room_id="room1", index=0)
    assert result is not None
    assert msgpack.unpackb(result[b"a"]) == 99
    assert msgpack.unpackb(result[b"b"]) == 2


@pytest.mark.asyncio
async def test_memory_merge_item_adds_new_keys(storage: InMemoryStorage) -> None:
    """Merge adds keys that don't exist in the original frame."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])
    partial: RawFrame = {b"new_key": msgpack.packb("hello")}  # type: ignore[reportAssignmentType]
    await storage.merge_item(room_id="room1", index=0, partial=partial)

    result = await storage.get(room_id="room1", index=0)
    assert result is not None
    assert msgpack.unpackb(result[b"a"]) == 1
    assert msgpack.unpackb(result[b"new_key"]) == "hello"


@pytest.mark.asyncio
async def test_memory_merge_item_preserves_untouched_keys(
    storage: InMemoryStorage,
) -> None:
    """Merge does not affect keys not in the partial."""
    await storage.extend(room_id="room1", frames=[{"a": 1, "b": 2, "c": 3}])
    partial: RawFrame = {b"b": msgpack.packb(99)}  # type: ignore[reportAssignmentType]
    await storage.merge_item(room_id="room1", index=0, partial=partial)

    result = await storage.get(room_id="room1", index=0)
    assert result is not None
    assert msgpack.unpackb(result[b"a"]) == 1
    assert msgpack.unpackb(result[b"b"]) == 99
    assert msgpack.unpackb(result[b"c"]) == 3


@pytest.mark.asyncio
async def test_memory_merge_item_out_of_bounds_raises(
    storage: InMemoryStorage,
) -> None:
    """Merging at invalid index raises IndexError."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])

    with pytest.raises(IndexError):
        await storage.merge_item(
            room_id="room1",
            index=5,
            partial={b"a": msgpack.packb(1)},  # type: ignore[reportArgumentType]
        )


@pytest.mark.asyncio
async def test_memory_merge_item_on_empty_raises(storage: InMemoryStorage) -> None:
    """Merging on empty room raises IndexError."""
    with pytest.raises(IndexError):
        await storage.merge_item(
            room_id="room1",
            index=0,
            partial={b"a": msgpack.packb(1)},  # type: ignore[reportArgumentType]
        )


# --- Room Isolation ---


@pytest.mark.asyncio
async def test_memory_rooms_are_isolated(storage: InMemoryStorage) -> None:
    """Different rooms have separate storage."""
    await storage.extend(room_id="room1", frames=[{"room1": True}])
    await storage.extend(room_id="room2", frames=[{"room2": True}, {"room2b": True}])

    assert await storage.get_length(room_id="room1") == 1
    assert await storage.get_length(room_id="room2") == 2
    assert await storage.get(room_id="room1", index=0) == make_raw_frame(
        {"room1": True}
    )
    assert await storage.get(room_id="room2", index=0) == make_raw_frame(
        {"room2": True}
    )


# --- DeleteRange ---


@pytest.mark.asyncio
async def test_memory_delete_range(storage: InMemoryStorage) -> None:
    """Deleting shifts subsequent frames."""
    frames = [{"a": 1}, {"b": 2}, {"c": 3}, {"d": 4}]
    await storage.extend(room_id="room1", frames=frames)
    await storage.delete_range(room_id="room1", start=1, stop=3)

    assert await storage.get_length(room_id="room1") == 2
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) == make_raw_frame({"d": 4})


@pytest.mark.asyncio
async def test_memory_delete_range_all(storage: InMemoryStorage) -> None:
    """Can delete all frames via range."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    await storage.delete_range(room_id="room1", start=0, stop=2)

    assert await storage.get_length(room_id="room1") == 0


@pytest.mark.asyncio
async def test_memory_delete_range_empty_room(storage: InMemoryStorage) -> None:
    """Deleting from empty room is safe (no-op)."""
    await storage.delete_range(room_id="room1", start=0, stop=10)
    assert await storage.get_length(room_id="room1") == 0


# --- Clear ---


@pytest.mark.asyncio
async def test_memory_clear(storage: InMemoryStorage) -> None:
    """Clear one room doesn't affect others."""
    await storage.extend(room_id="room1", frames=[{"room1": True}])
    await storage.extend(room_id="room2", frames=[{"room2": True}])

    await storage.clear(room_id="room1")

    assert await storage.get_length(room_id="room1") == 0
    assert await storage.get_length(room_id="room2") == 1
    assert await storage.get(room_id="room2", index=0) == make_raw_frame(
        {"room2": True}
    )


@pytest.mark.asyncio
async def test_memory_clear_empty_room(storage: InMemoryStorage) -> None:
    """Clearing empty room is safe (no-op)."""
    await storage.clear(room_id="room1")
    assert await storage.get_length(room_id="room1") == 0


# --- Reserve ---


@pytest.mark.asyncio
async def test_memory_reserve_grows_length(storage: InMemoryStorage) -> None:
    """Reserve grows length to the target count."""
    await storage.reserve(room_id="room1", count=10)
    assert await storage.get_length(room_id="room1") == 10


@pytest.mark.asyncio
async def test_memory_reserve_is_grow_only(storage: InMemoryStorage) -> None:
    """Reserve with smaller count is a no-op."""
    await storage.reserve(room_id="room1", count=10)
    await storage.reserve(room_id="room1", count=5)
    assert await storage.get_length(room_id="room1") == 10


@pytest.mark.asyncio
async def test_memory_reserved_slots_return_none(storage: InMemoryStorage) -> None:
    """Reserved-but-unfilled slots return None from get()."""
    await storage.reserve(room_id="room1", count=5)
    assert await storage.get(room_id="room1", index=0) is None
    assert await storage.get(room_id="room1", index=4) is None


@pytest.mark.asyncio
async def test_memory_set_item_on_reserved_slot(storage: InMemoryStorage) -> None:
    """Can set_item on a reserved-but-unfilled slot."""
    await storage.reserve(room_id="room1", count=5)
    await storage.set_item(room_id="room1", index=3, frame={"x": 1})
    assert await storage.get(room_id="room1", index=3) == make_raw_frame({"x": 1})


@pytest.mark.asyncio
async def test_memory_get_range_preserves_none(storage: InMemoryStorage) -> None:
    """get_range returns positional list with None for empty slots."""
    await storage.reserve(room_id="room1", count=5)
    await storage.set_item(room_id="room1", index=2, frame={"x": 1})
    result = await storage.get_range(room_id="room1", start=0, stop=5)
    assert result == [None, None, make_raw_frame({"x": 1}), None, None]


@pytest.mark.asyncio
async def test_memory_get_many_preserves_none(storage: InMemoryStorage) -> None:
    """get_many returns positional list with None for empty slots."""
    await storage.reserve(room_id="room1", count=5)
    await storage.set_item(room_id="room1", index=1, frame={"y": 2})
    result = await storage.get_many(room_id="room1", indices=[0, 1, 2])
    assert result == [None, make_raw_frame({"y": 2}), None]


# --- RemoveItems ---


@pytest.mark.asyncio
async def test_memory_remove_items_makes_none(storage: InMemoryStorage) -> None:
    """Removed items return None from get()."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}, {"c": 3}])
    await storage.remove_items(room_id="room1", indices=[1])
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) is None
    assert await storage.get(room_id="room1", index=2) == make_raw_frame({"c": 3})


@pytest.mark.asyncio
async def test_memory_remove_items_preserves_length(storage: InMemoryStorage) -> None:
    """remove_items does not change get_length()."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    await storage.remove_items(room_id="room1", indices=[0])
    assert await storage.get_length(room_id="room1") == 2


@pytest.mark.asyncio
async def test_memory_remove_out_of_bounds_is_noop(storage: InMemoryStorage) -> None:
    """Removing out-of-bounds index is safe."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])
    await storage.remove_items(room_id="room1", indices=[999])
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})


# --- Close ---


@pytest.mark.asyncio
async def test_memory_close_clears_all_data(storage: InMemoryStorage) -> None:
    """Close clears all data from all rooms."""
    await storage.extend(room_id="room1", frames=[{"room1": True}])
    await storage.extend(room_id="room2", frames=[{"room2": True}])

    await storage.close()

    assert await storage.get_length(room_id="room1") == 0
    assert await storage.get_length(room_id="room2") == 0
