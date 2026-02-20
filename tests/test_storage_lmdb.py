"""Tests for LMDB storage backend with async safety."""

import asyncio
import shutil
import tempfile
import time
from collections.abc import AsyncGenerator
from pathlib import Path

import msgpack
import pytest
import pytest_asyncio
from conftest import make_raw_frame

from zndraw.storage import LMDBStorage, StorageBackend
from zndraw.storage.base import RawFrame


@pytest_asyncio.fixture
async def storage() -> AsyncGenerator[LMDBStorage, None]:
    """Create a fresh LMDBStorage instance for each test."""
    tmpdir = Path(tempfile.mkdtemp())
    s = LMDBStorage(path=tmpdir / "test.lmdb", map_size=10_000_000)
    yield s
    await s.close()
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest_asyncio.fixture
async def storage_path() -> AsyncGenerator[Path, None]:
    """Create a temp directory and yield its path for persistence tests."""
    tmpdir = Path(tempfile.mkdtemp())
    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


# --- Interface ---


def test_lmdb_is_subclass_of_storage_backend() -> None:
    """LMDBStorage is a proper subclass of StorageBackend."""
    assert issubclass(LMDBStorage, StorageBackend)


def test_lmdb_can_instantiate() -> None:
    """LMDBStorage can be instantiated."""
    tmpdir = Path(tempfile.mkdtemp())
    try:
        storage = LMDBStorage(path=tmpdir / "test.lmdb", map_size=10_000_000)
        assert isinstance(storage, StorageBackend)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# --- Basics ---


@pytest.mark.asyncio
async def test_lmdb_empty_room_has_zero_length(storage: LMDBStorage) -> None:
    """New room has no frames."""
    length = await storage.get_length(room_id="room1")
    assert length == 0


@pytest.mark.asyncio
async def test_lmdb_get_nonexistent_returns_none(storage: LMDBStorage) -> None:
    """Get on empty room returns None."""
    result = await storage.get(room_id="room1", index=0)
    assert result is None


@pytest.mark.asyncio
async def test_lmdb_get_out_of_bounds_returns_none(storage: LMDBStorage) -> None:
    """Get with out of bounds index returns None."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])
    result = await storage.get(room_id="room1", index=5)
    assert result is None


# --- Extend and Get ---


@pytest.mark.asyncio
async def test_lmdb_extend_and_get(storage: LMDBStorage) -> None:
    """Extend then retrieve frames."""
    frames = [{"a": 1}, {"b": 2}, {"c": 3}]
    new_length = await storage.extend(room_id="room1", frames=frames)

    assert new_length == 3
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) == make_raw_frame({"b": 2})
    assert await storage.get(room_id="room1", index=2) == make_raw_frame({"c": 3})


@pytest.mark.asyncio
async def test_lmdb_extend_multiple_times(storage: LMDBStorage) -> None:
    """Extend can be called multiple times."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])
    new_length = await storage.extend(room_id="room1", frames=[{"b": 2}, {"c": 3}])

    assert new_length == 3
    assert await storage.get_length(room_id="room1") == 3


@pytest.mark.asyncio
async def test_lmdb_get_range(storage: LMDBStorage) -> None:
    """Range slicing works."""
    frames = [{"a": 1}, {"b": 2}, {"c": 3}, {"d": 4}]
    await storage.extend(room_id="room1", frames=frames)

    result = await storage.get_range(room_id="room1", start=1, stop=3)
    assert result == [make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]


@pytest.mark.asyncio
async def test_lmdb_get_range_empty(storage: LMDBStorage) -> None:
    """Range on empty room returns empty list."""
    result = await storage.get_range(room_id="room1", start=0, stop=10)
    assert result == []


@pytest.mark.asyncio
async def test_lmdb_get_range_partial(storage: LMDBStorage) -> None:
    """Range beyond bounds returns available frames only."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    result = await storage.get_range(room_id="room1", start=0, stop=10)
    assert result == [make_raw_frame({"a": 1}), make_raw_frame({"b": 2})]


# --- SetItem ---


@pytest.mark.asyncio
async def test_lmdb_set_item(storage: LMDBStorage) -> None:
    """Can overwrite a frame."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    await storage.set_item(room_id="room1", index=1, frame={"updated": True})

    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) == make_raw_frame(
        {"updated": True}
    )


@pytest.mark.asyncio
async def test_lmdb_set_item_out_of_bounds_raises(storage: LMDBStorage) -> None:
    """Setting frame at invalid index raises IndexError."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])

    with pytest.raises(IndexError):
        await storage.set_item(room_id="room1", index=5, frame={"bad": True})


@pytest.mark.asyncio
async def test_lmdb_set_item_on_empty_raises(storage: LMDBStorage) -> None:
    """Setting frame on empty room raises IndexError."""
    with pytest.raises(IndexError):
        await storage.set_item(room_id="room1", index=0, frame={"bad": True})


# --- MergeItem ---


@pytest.mark.asyncio
async def test_lmdb_merge_item_updates_existing_keys(storage: LMDBStorage) -> None:
    """Merge overwrites existing keys."""
    await storage.extend(room_id="room1", frames=[{"a": 1, "b": 2}])
    partial: RawFrame = {b"a": msgpack.packb(99)}  # type: ignore[reportAssignmentType]
    await storage.merge_item(room_id="room1", index=0, partial=partial)

    result = await storage.get(room_id="room1", index=0)
    assert result is not None
    assert msgpack.unpackb(result[b"a"]) == 99
    assert msgpack.unpackb(result[b"b"]) == 2


@pytest.mark.asyncio
async def test_lmdb_merge_item_adds_new_keys(storage: LMDBStorage) -> None:
    """Merge adds keys that don't exist in the original frame."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])
    partial: RawFrame = {b"new_key": msgpack.packb("hello")}  # type: ignore[reportAssignmentType]
    await storage.merge_item(room_id="room1", index=0, partial=partial)

    result = await storage.get(room_id="room1", index=0)
    assert result is not None
    assert msgpack.unpackb(result[b"a"]) == 1
    assert msgpack.unpackb(result[b"new_key"]) == "hello"


@pytest.mark.asyncio
async def test_lmdb_merge_item_preserves_untouched_keys(
    storage: LMDBStorage,
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
async def test_lmdb_merge_item_out_of_bounds_raises(storage: LMDBStorage) -> None:
    """Merging at invalid index raises IndexError."""
    await storage.extend(room_id="room1", frames=[{"a": 1}])

    with pytest.raises(IndexError):
        await storage.merge_item(
            room_id="room1",
            index=5,
            partial={b"a": msgpack.packb(1)},  # type: ignore[reportArgumentType]
        )


@pytest.mark.asyncio
async def test_lmdb_merge_item_on_empty_raises(storage: LMDBStorage) -> None:
    """Merging on empty room raises IndexError."""
    with pytest.raises(IndexError):
        await storage.merge_item(
            room_id="room1",
            index=0,
            partial={b"a": msgpack.packb(1)},  # type: ignore[reportArgumentType]
        )


# --- Room Isolation ---


@pytest.mark.asyncio
async def test_lmdb_rooms_are_isolated(storage: LMDBStorage) -> None:
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
async def test_lmdb_delete_range(storage: LMDBStorage) -> None:
    """Deleting shifts subsequent frames."""
    frames = [{"a": 1}, {"b": 2}, {"c": 3}, {"d": 4}]
    await storage.extend(room_id="room1", frames=frames)
    await storage.delete_range(room_id="room1", start=1, stop=3)

    assert await storage.get_length(room_id="room1") == 2
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) == make_raw_frame({"d": 4})


@pytest.mark.asyncio
async def test_lmdb_delete_range_all(storage: LMDBStorage) -> None:
    """Can delete all frames via range."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    await storage.delete_range(room_id="room1", start=0, stop=2)

    assert await storage.get_length(room_id="room1") == 0


@pytest.mark.asyncio
async def test_lmdb_delete_range_empty_room(storage: LMDBStorage) -> None:
    """Deleting from empty room is safe (no-op)."""
    await storage.delete_range(room_id="room1", start=0, stop=10)
    assert await storage.get_length(room_id="room1") == 0


# --- Clear ---


@pytest.mark.asyncio
async def test_lmdb_clear(storage: LMDBStorage) -> None:
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
async def test_lmdb_clear_empty_room(storage: LMDBStorage) -> None:
    """Clearing empty room is safe (no-op)."""
    await storage.clear(room_id="room1")
    assert await storage.get_length(room_id="room1") == 0


# --- Reserve ---


@pytest.mark.asyncio
async def test_lmdb_reserve_grows_length(storage: LMDBStorage) -> None:
    """Reserve grows length to the target count."""
    await storage.reserve(room_id="room1", count=10)
    assert await storage.get_length(room_id="room1") == 10


@pytest.mark.asyncio
async def test_lmdb_reserve_is_grow_only(storage: LMDBStorage) -> None:
    """Reserve with smaller count is a no-op."""
    await storage.reserve(room_id="room1", count=10)
    await storage.reserve(room_id="room1", count=5)
    assert await storage.get_length(room_id="room1") == 10


@pytest.mark.asyncio
async def test_lmdb_reserved_slots_return_none(storage: LMDBStorage) -> None:
    """Reserved-but-unfilled slots return None from get()."""
    await storage.reserve(room_id="room1", count=5)
    assert await storage.get(room_id="room1", index=0) is None
    assert await storage.get(room_id="room1", index=4) is None


@pytest.mark.asyncio
async def test_lmdb_set_item_on_reserved_slot(storage: LMDBStorage) -> None:
    """Can set_item on a reserved-but-unfilled slot."""
    await storage.reserve(room_id="room1", count=5)
    await storage.set_item(room_id="room1", index=3, frame={"x": 1})
    assert await storage.get(room_id="room1", index=3) == make_raw_frame({"x": 1})


# --- RemoveItems ---


@pytest.mark.asyncio
async def test_lmdb_remove_items_makes_none(storage: LMDBStorage) -> None:
    """Removed items return None from get()."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}, {"c": 3}])
    await storage.remove_items(room_id="room1", indices=[1])
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"a": 1})
    assert await storage.get(room_id="room1", index=1) is None
    assert await storage.get(room_id="room1", index=2) == make_raw_frame({"c": 3})


@pytest.mark.asyncio
async def test_lmdb_remove_items_preserves_length(storage: LMDBStorage) -> None:
    """remove_items does not change get_length()."""
    await storage.extend(room_id="room1", frames=[{"a": 1}, {"b": 2}])
    await storage.remove_items(room_id="room1", indices=[0])
    assert await storage.get_length(room_id="room1") == 2


@pytest.mark.asyncio
async def test_lmdb_remove_nonexistent_is_noop(storage: LMDBStorage) -> None:
    """Removing index with no data is safe."""
    await storage.reserve(room_id="room1", count=5)
    await storage.remove_items(room_id="room1", indices=[3])
    assert await storage.get(room_id="room1", index=3) is None


# --- Persistence ---


@pytest.mark.asyncio
async def test_lmdb_persistence(storage_path: Path) -> None:
    """Data persists across storage instances."""
    db_path = storage_path / "persist.lmdb"

    # Write data with first instance
    storage1 = LMDBStorage(path=db_path, map_size=10_000_000)
    await storage1.extend(room_id="room1", frames=[{"persistent": True}, {"data": 42}])
    await storage1.close()

    # Read with new instance
    storage2 = LMDBStorage(path=db_path, map_size=10_000_000)
    try:
        assert await storage2.get_length(room_id="room1") == 2
        assert await storage2.get(room_id="room1", index=0) == make_raw_frame(
            {"persistent": True}
        )
        assert await storage2.get(room_id="room1", index=1) == make_raw_frame(
            {"data": 42}
        )
    finally:
        await storage2.close()


# --- Async Safety ---


@pytest.mark.asyncio
async def test_lmdb_async_safety_no_blocking(storage: LMDBStorage) -> None:
    """Concurrent operations don't block the event loop.

    This test verifies that LMDB operations run in a thread pool and don't
    block other async tasks from executing.
    """
    # Add some data to work with
    frames = [{"frame": i} for i in range(100)]
    await storage.extend(room_id="room1", frames=frames)

    # Track event loop responsiveness
    heartbeats: list[float] = []
    stop_heartbeat = False

    async def heartbeat_task() -> None:
        """Simulates event loop activity (like Socket.IO heartbeats)."""
        nonlocal stop_heartbeat
        while not stop_heartbeat:
            heartbeats.append(time.time())
            await asyncio.sleep(0.01)  # 10ms heartbeat

    async def storage_operations() -> None:
        """Perform multiple storage operations."""
        get_tasks = [storage.get(room_id="room1", index=i % 100) for i in range(50)]
        range_tasks = [
            storage.get_range(room_id="room1", start=0, stop=50) for _ in range(10)
        ]
        await asyncio.gather(*get_tasks, *range_tasks)

    # Run heartbeat and storage operations concurrently
    heartbeat = asyncio.create_task(heartbeat_task())
    await storage_operations()
    stop_heartbeat = True
    await heartbeat

    # Verify heartbeats continued during storage operations
    # If LMDB was blocking, heartbeats would stop or have large gaps
    assert len(heartbeats) > 0, "Heartbeat task should have run"

    # Check for gaps > 100ms (would indicate blocking)
    if len(heartbeats) > 1:
        gaps = [heartbeats[i + 1] - heartbeats[i] for i in range(len(heartbeats) - 1)]
        max_gap = max(gaps)
        assert max_gap < 0.1, f"Event loop blocked for {max_gap:.3f}s"


@pytest.mark.asyncio
async def test_lmdb_concurrent_rooms(storage: LMDBStorage) -> None:
    """Multiple rooms can be accessed concurrently."""
    # Create data in multiple rooms
    for room_id in range(1, 6):
        frames = [{"room": room_id, "frame": i} for i in range(10)]
        await storage.extend(room_id=f"room{room_id}", frames=frames)

    # Access all rooms concurrently
    async def access_room(room_id: int) -> list[RawFrame | None]:
        results: list[RawFrame | None] = []
        for i in range(10):
            frame = await storage.get(room_id=f"room{room_id}", index=i)
            results.append(frame)
        return results

    tasks = [access_room(room_id) for room_id in range(1, 6)]
    all_results = await asyncio.gather(*tasks)

    # Verify each room got correct data
    for room_idx, results in enumerate(all_results):
        room_id = room_idx + 1
        for frame_idx, frame in enumerate(results):
            assert frame == make_raw_frame({"room": room_id, "frame": frame_idx})


# --- Close ---


@pytest.mark.asyncio
async def test_lmdb_close_shuts_down_executor(storage_path: Path) -> None:
    """Close shuts down the thread pool executor."""
    storage = LMDBStorage(path=storage_path / "close.lmdb", map_size=10_000_000)
    await storage.extend(room_id="room1", frames=[{"test": True}])

    # Verify storage works before close
    assert await storage.get(room_id="room1", index=0) == make_raw_frame({"test": True})

    await storage.close()

    # Executor should be shut down
    assert storage._executor._shutdown
