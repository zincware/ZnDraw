"""LMDB storage backend with async-safe blocking calls.

Python's LMDB bindings are synchronous (blocking). Without run_in_executor,
large trajectory reads would block the event loop, causing Socket.IO heartbeat
timeouts and mass client disconnections.

Frame data is stored as msgpack with raw=True to preserve bytes keys/values:
- Keys: b"arrays.positions", b"cell", etc.
- Values: msgpack-numpy encoded bytes

This matches the reference implementation and enables efficient binary transfer.
"""

import asyncio
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, cast

import lmdb
import msgpack

from .base import RawFrame, StorageBackend, to_raw_frame


class LMDBStorage(StorageBackend):
    """LMDB storage with async-safe blocking calls.

    Uses ThreadPoolExecutor to prevent event loop blocking.
    Data is stored as msgpack-encoded dict[bytes, bytes], keyed by "frame:{room_id}:{index}".
    Room lengths are tracked in "meta:room:{room_id}:length".

    Frame data preserves raw bytes format for efficient binary transfer.
    """

    def __init__(
        self,
        path: Path,
        map_size: int = 1_073_741_824,
        max_workers: int = 4,
    ) -> None:
        """Initialize LMDB storage.

        Parameters
        ----------
        path
            Path to the LMDB database directory.
        map_size
            Maximum database size in bytes (default 1GB).
        max_workers
            Maximum threads in the executor pool.
        """
        self.path = Path(path)
        self.map_size = map_size
        self._env: lmdb.Environment | None = None
        self._executor = ThreadPoolExecutor(
            max_workers=max_workers, thread_name_prefix="lmdb"
        )

    def _get_env(self) -> lmdb.Environment:
        """Lazily initialize and return the LMDB environment.

        Returns
        -------
        lmdb.Environment
            The LMDB environment instance.

        Raises
        ------
        RuntimeError
            If LMDB environment cannot be opened.
        """
        if self._env is None:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            env = lmdb.open(
                str(self.path),
                map_size=self.map_size,
                subdir=True,
                max_dbs=0,
            )
            if env is None:
                raise RuntimeError(f"Failed to open LMDB environment at {self.path}")
            self._env = cast("lmdb.Environment", env)
        return self._env

    @staticmethod
    def _frame_key(room_id: str, index: int) -> bytes:
        """Generate the key for a frame.

        Parameters
        ----------
        room_id
            The room identifier.
        index
            The frame index.

        Returns
        -------
        bytes
            The encoded key bytes.
        """
        return f"frame:{room_id}:{index}".encode()

    @staticmethod
    def _length_key(room_id: str) -> bytes:
        """Generate the key for a room's length metadata.

        Parameters
        ----------
        room_id
            The room identifier.

        Returns
        -------
        bytes
            The encoded key bytes.
        """
        return f"meta:room:{room_id}:length".encode()

    def _sync_get(self, room_id: str, index: int) -> RawFrame | None:
        """Synchronous get - runs in thread pool.

        Parameters
        ----------
        room_id
            The room identifier.
        index
            The frame index.

        Returns
        -------
        RawFrame | None
            The frame data as dict[bytes, bytes] or None if not found.
        """
        env = self._get_env()
        with env.begin() as txn:
            data = txn.get(self._frame_key(room_id, index))
            if data is None:
                return None
            # raw=True preserves bytes keys and values
            return msgpack.unpackb(data, raw=True)

    async def get(self, room_id: str, index: int) -> RawFrame | None:
        """Get a single frame by index.

        Parameters
        ----------
        room_id
            The room identifier.
        index
            The frame index (0-based).

        Returns
        -------
        RawFrame | None
            The frame data as dict[bytes, bytes], or None if not found.
        """
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(
            self._executor, self._sync_get, room_id, index
        )

    def _sync_get_range(
        self, room_id: str, start: int, stop: int | None
    ) -> list[RawFrame | None]:
        """Synchronous get_range - runs in thread pool.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        start
            Start index (inclusive).
        stop
            Stop index (exclusive), or None for all remaining frames.

        Returns
        -------
        list[RawFrame]
            List of frames as dict[bytes, bytes] in the specified range.
        """
        env = self._get_env()
        length = self._sync_get_length(room_id)
        actual_stop = min(stop, length) if stop is not None else length

        results: list[RawFrame | None] = []
        with env.begin() as txn:
            for i in range(start, actual_stop):
                data = txn.get(self._frame_key(room_id, i))
                if data is not None:
                    results.append(msgpack.unpackb(data, raw=True))
                else:
                    results.append(None)
        return results

    async def get_range(
        self, room_id: str, start: int, stop: int | None
    ) -> list[RawFrame | None]:
        """Get a range of frames [start, stop).

        Parameters
        ----------
        room_id
            The room identifier.
        start
            Start index (inclusive).
        stop
            Stop index (exclusive).

        Returns
        -------
        list[RawFrame | None]
            Positional list. ``None`` for empty slots.
        """
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(
            self._executor, self._sync_get_range, room_id, start, stop
        )

    def _sync_get_many(self, room_id: str, indices: list[int]) -> list[RawFrame | None]:
        """Synchronous get_many - runs in thread pool.

        Parameters
        ----------
        room_id
            The room identifier.
        indices
            List of frame indices to retrieve.

        Returns
        -------
        list[RawFrame | None]
            Positional list. ``None`` for out-of-bounds or empty slots.
        """
        env = self._get_env()
        length = self._sync_get_length(room_id)

        results: list[RawFrame | None] = []
        with env.begin() as txn:
            for i in indices:
                if 0 <= i < length:
                    data = txn.get(self._frame_key(room_id, i))
                    if data is not None:
                        results.append(msgpack.unpackb(data, raw=True))
                    else:
                        results.append(None)
                else:
                    results.append(None)
        return results

    async def get_many(self, room_id: str, indices: list[int]) -> list[RawFrame | None]:
        """Get multiple frames by specific indices.

        Parameters
        ----------
        room_id
            The room identifier.
        indices
            List of frame indices to retrieve.

        Returns
        -------
        list[RawFrame | None]
            Positional list. ``None`` for out-of-bounds or empty slots.
        """
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(
            self._executor, self._sync_get_many, room_id, indices
        )

    def _sync_extend(
        self, room_id: str, frames: list[dict[str, Any]] | list[RawFrame]
    ) -> int:
        """Synchronous extend - runs in thread pool.

        Converts input frames to raw bytes format before storing.

        Parameters
        ----------
        room_id
            The room identifier.
        frames
            List of frame dictionaries to append.

        Returns
        -------
        int
            New total frame count for the room.
        """
        env = self._get_env()
        current_length = self._sync_get_length(room_id)

        with env.begin(write=True) as txn:
            for i, frame in enumerate(frames):
                key = self._frame_key(room_id, current_length + i)
                # Convert to raw bytes format and store
                raw_frame = to_raw_frame(frame)
                txn.put(key, msgpack.packb(raw_frame, use_bin_type=True))

            new_length = current_length + len(frames)
            txn.put(self._length_key(room_id), str(new_length).encode())

        return new_length

    async def extend(
        self, room_id: str, frames: list[dict[str, Any]] | list[RawFrame]
    ) -> int:
        """Append frames to storage.

        Converts input frames to raw bytes format before storing.

        Parameters
        ----------
        room_id
            The room identifier.
        frames
            List of frame dictionaries to append.

        Returns
        -------
        int
            New total frame count for the room.
        """
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(
            self._executor, self._sync_extend, room_id, frames
        )

    def _sync_set_item(
        self, room_id: str, index: int, frame: dict[str, Any] | RawFrame
    ) -> None:
        """Synchronous set_item - runs in thread pool."""
        env = self._get_env()
        length = self._sync_get_length(room_id)

        if index < 0 or index >= length:
            raise IndexError(f"Frame index {index} out of bounds for room {room_id}")

        with env.begin(write=True) as txn:
            key = self._frame_key(room_id, index)
            raw_frame = to_raw_frame(frame)
            txn.put(key, msgpack.packb(raw_frame, use_bin_type=True))

    async def set_item(
        self, room_id: str, index: int, frame: dict[str, Any] | RawFrame
    ) -> None:
        """Set a frame at a specific index.

        Raises IndexError if index is out of bounds.
        """
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(
            self._executor, self._sync_set_item, room_id, index, frame
        )

    def _sync_merge_item(self, room_id: str, index: int, partial: RawFrame) -> None:
        """Synchronous merge_item - runs in thread pool.

        Reads existing frame, merges partial data, writes back in a single
        LMDB write transaction (atomic).
        """
        env = self._get_env()
        length = self._sync_get_length(room_id)

        if index < 0 or index >= length:
            raise IndexError(f"Frame index {index} out of bounds for room {room_id}")

        with env.begin(write=True) as txn:
            key = self._frame_key(room_id, index)
            data = txn.get(key)
            existing: RawFrame = msgpack.unpackb(data, raw=True) if data else {}
            existing.update(partial)
            txn.put(key, msgpack.packb(existing, use_bin_type=True))

    async def merge_item(self, room_id: str, index: int, partial: RawFrame) -> None:
        """Merge partial frame data into existing frame at index.

        Atomic read-merge-write inside a single LMDB write transaction.
        """
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(
            self._executor, self._sync_merge_item, room_id, index, partial
        )

    def _sync_get_length(self, room_id: str) -> int:
        """Synchronous get_length - runs in thread pool.

        Parameters
        ----------
        room_id
            The room identifier.

        Returns
        -------
        int
            The number of frames stored for this room.
        """
        env = self._get_env()
        with env.begin() as txn:
            data = txn.get(self._length_key(room_id))
            if data is None:
                return 0
            return int(bytes(data).decode())

    async def get_length(self, room_id: str) -> int:
        """Get total frame count for a room.

        Parameters
        ----------
        room_id
            The room identifier.

        Returns
        -------
        int
            The number of frames stored for this room.
        """
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(
            self._executor, self._sync_get_length, room_id
        )

    def _sync_delete_range(self, room_id: str, start: int, stop: int) -> None:
        """Synchronous delete_range - runs in thread pool.

        Frames after the deleted range are shifted to fill the gap.

        Parameters
        ----------
        room_id
            The room identifier.
        start
            Start index (inclusive).
        stop
            Stop index (exclusive).
        """
        env = self._get_env()
        length = self._sync_get_length(room_id)

        if length == 0 or start >= length:
            return

        actual_stop = min(stop, length)
        delete_count = actual_stop - start

        with env.begin(write=True) as txn:
            # Shift frames after the deleted range
            for i in range(actual_stop, length):
                old_key = self._frame_key(room_id, i)
                new_key = self._frame_key(room_id, i - delete_count)
                data = txn.get(old_key)
                if data is not None:
                    txn.put(new_key, data)
                    txn.delete(old_key)

            # Delete frames in the range that weren't overwritten by shifting
            # (only needed when there are no frames to shift)
            if actual_stop >= length:
                for i in range(start, actual_stop):
                    txn.delete(self._frame_key(room_id, i))

            # Update length
            new_length = length - delete_count
            txn.put(self._length_key(room_id), str(new_length).encode())

    async def delete_range(self, room_id: str, start: int, stop: int) -> None:
        """Delete a range of frames [start, stop).

        Frames after the deleted range are shifted to fill the gap.

        Parameters
        ----------
        room_id
            The room identifier.
        start
            Start index (inclusive).
        stop
            Stop index (exclusive).
        """
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(
            self._executor, self._sync_delete_range, room_id, start, stop
        )

    def _sync_clear(self, room_id: str) -> None:
        """Synchronous clear - runs in thread pool.

        Parameters
        ----------
        room_id
            The room identifier.
        """
        env = self._get_env()
        length = self._sync_get_length(room_id)

        with env.begin(write=True) as txn:
            # Delete all frames for this room
            for i in range(length):
                txn.delete(self._frame_key(room_id, i))

            # Reset length to 0
            txn.put(self._length_key(room_id), b"0")

    async def clear(self, room_id: str) -> None:
        """Delete all frames for a room.

        Parameters
        ----------
        room_id
            The room identifier.
        """
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(self._executor, self._sync_clear, room_id)

    def _sync_reserve(self, room_id: str, count: int) -> None:
        """Synchronous reserve - runs in thread pool."""
        env = self._get_env()
        current = self._sync_get_length(room_id)
        if count > current:
            with env.begin(write=True) as txn:
                txn.put(self._length_key(room_id), str(count).encode())

    async def reserve(self, room_id: str, count: int) -> None:
        """Pre-allocate space for count frames. Grow only."""
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(self._executor, self._sync_reserve, room_id, count)

    def _sync_remove_items(self, room_id: str, indices: list[int]) -> None:
        """Synchronous remove_items - runs in thread pool."""
        env = self._get_env()
        with env.begin(write=True) as txn:
            for idx in indices:
                txn.delete(self._frame_key(room_id, idx))

    async def remove_items(self, room_id: str, indices: list[int]) -> None:
        """Remove frames at indices without shifting."""
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(
            self._executor, self._sync_remove_items, room_id, indices
        )

    def _sync_close(self) -> None:
        """Synchronous close - runs in thread pool."""
        if self._env is not None:
            self._env.close()
            self._env = None

    async def close(self) -> None:
        """Clean up resources.

        Closes the LMDB environment and shuts down the thread pool executor.
        """
        loop = asyncio.get_running_loop()
        await loop.run_in_executor(self._executor, self._sync_close)
        self._executor.shutdown(wait=True)
