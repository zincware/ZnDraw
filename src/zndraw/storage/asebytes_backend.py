"""Asebytes-backed frame storage with room-scoped AsyncBlobIO instances."""

import base64
from typing import Any

import msgpack
from asebytes import AsyncBlobIO
from asebytes._async_adapters import AsyncObjectToBlobReadWriteAdapter
from asebytes._async_backends import AsyncReadWriteBackend

# Type alias for raw frame data: dict with bytes keys and bytes values
# Keys are like b"arrays.positions", b"cell", etc.
# Values are msgpack-numpy encoded bytes
RawFrame = dict[bytes, bytes]


def to_raw_frame(frame: dict[str, Any] | RawFrame) -> RawFrame:
    """Convert input frame dict to raw bytes format.

    Handles two input formats:
    1. Already raw: dict[bytes, bytes] -- pass through
    2. Base64 encoded: dict[str, str] with keys like "b64:..." -- decode

    Parameters
    ----------
    frame
        Input frame data (dict[str, Any] or dict[bytes, bytes]).
    """
    if not frame:
        return {}

    first_key = next(iter(frame.keys()))

    if isinstance(first_key, bytes):
        return frame  # type: ignore[return-value]

    result: RawFrame = {}
    for k, v in frame.items():
        key_bytes: bytes
        if isinstance(k, str) and k.startswith("b64:"):
            key_bytes = base64.b64decode(k[4:])
        elif isinstance(k, str):
            key_bytes = k.encode()
        elif isinstance(k, bytes):
            key_bytes = k
        else:
            key_bytes = str(k).encode()

        val_bytes: bytes
        if isinstance(v, str):
            val_bytes = base64.b64decode(v)
        elif isinstance(v, bytes):
            val_bytes = v
        else:
            packed = msgpack.packb(v)
            val_bytes = packed if packed is not None else b""

        result[key_bytes] = val_bytes

    return result


def _create_blob_backend(
    uri: str, group: str | None = None
) -> AsyncReadWriteBackend[bytes, bytes]:
    """Create an async blob backend from a URI string.

    All backends from asebytes registries are object-level (str keys).
    This function wraps them to blob-level (bytes keys) via
    AsyncObjectToBlobReadWriteAdapter, converting sync backends to
    async first when needed.

    Parameters
    ----------
    uri
        Storage backend URI string.
    group
        Optional group/namespace for data isolation (e.g., room_id).
    """
    from asebytes._registry import get_async_backend_cls

    cls = get_async_backend_cls(uri, readonly=False)
    if hasattr(cls, "from_uri"):
        backend = cls.from_uri(uri, group=group)
    else:
        backend = cls(uri, group=group)

    # Ensure async
    if not isinstance(backend, AsyncReadWriteBackend):
        from asebytes._async_backends import sync_to_async

        backend = sync_to_async(backend)  # type: ignore[arg-type]

    # All registry backends are object-level (str keys) → wrap to blob
    return AsyncObjectToBlobReadWriteAdapter(backend)  # type: ignore[arg-type]


class AsebytesStorage:
    """Room-scoped frame storage backed by asebytes.

    One AsyncBlobIO instance per room_id. Backend selected by URI:
    - ``memory://``              -- in-memory (dev/testing)
    - ``path/to/data.lmdb``      -- LMDB (single-node persistence)
    - ``mongodb://host:port/db`` -- MongoDB (distributed)

    Parameters
    ----------
    uri
        Storage backend URI string.
    """

    def __init__(self, uri: str) -> None:
        self._uri = uri
        self._rooms: dict[str, AsyncBlobIO] = {}

    def _get_io(self, room_id: str) -> AsyncBlobIO:
        """Get or create an AsyncBlobIO for a room."""
        if room_id not in self._rooms:
            backend = _create_blob_backend(self._uri, group=room_id)
            self._rooms[room_id] = AsyncBlobIO(backend)
        return self._rooms[room_id]

    async def get(self, room_id: str, index: int) -> RawFrame | None:
        """Get a single frame by index. None for OOB or empty slots."""
        io = self._get_io(room_id)
        length = await io.len()
        if index < 0 or index >= length:
            return None
        return await io.get(index)

    async def get_range(
        self, room_id: str, start: int, stop: int | None
    ) -> list[RawFrame | None]:
        """Get a range of frames [start, stop)."""
        io = self._get_io(room_id)
        length = await io.len()
        actual_stop = min(stop, length) if stop is not None else length
        if start >= actual_stop:
            return []
        indices = list(range(start, actual_stop))
        return await io._backend.get_many(indices)

    async def get_many(self, room_id: str, indices: list[int]) -> list[RawFrame | None]:
        """Get multiple frames by indices. None for OOB or empty slots."""
        io = self._get_io(room_id)
        length = await io.len()
        valid_map: dict[int, int] = {}  # position → index
        for pos, idx in enumerate(indices):
            if 0 <= idx < length:
                valid_map[pos] = idx
        if not valid_map:
            return [None] * len(indices)
        valid_indices = list(valid_map.values())
        rows = await io._backend.get_many(valid_indices)
        row_lookup = dict(zip(valid_indices, rows))
        return [row_lookup.get(idx) for idx in indices]

    async def extend(
        self, room_id: str, frames: list[dict[str, Any]] | list[RawFrame]
    ) -> int:
        """Append frames to storage. Returns new total count."""
        io = self._get_io(room_id)
        raw_frames = [to_raw_frame(f) for f in frames]
        return await io.extend(raw_frames)  # type: ignore[arg-type]

    async def set_item(
        self, room_id: str, index: int, frame: dict[str, Any] | RawFrame
    ) -> None:
        """Set a frame at a specific index. Raises IndexError if OOB."""
        io = self._get_io(room_id)
        raw = to_raw_frame(frame)
        await io._backend.set(index, raw)  # type: ignore[union-attr]

    async def merge_item(self, room_id: str, index: int, partial: RawFrame) -> None:
        """Merge partial frame data into existing frame at index."""
        io = self._get_io(room_id)
        await io._backend.update(index, partial)  # type: ignore[union-attr]

    async def get_length(self, room_id: str) -> int:
        """Get total frame count for a room."""
        io = self._get_io(room_id)
        return await io.len()

    async def delete_range(self, room_id: str, start: int, stop: int) -> None:
        """Delete a range of frames [start, stop) with index shifting."""
        io = self._get_io(room_id)
        await io._backend.delete_many(start, stop)  # type: ignore[union-attr]

    async def clear(self, room_id: str) -> None:
        """Delete all frames for a room."""
        io = self._get_io(room_id)
        await io.clear()

    async def reserve(self, room_id: str, count: int) -> None:
        """Pre-allocate count placeholder slots."""
        io = self._get_io(room_id)
        await io.reserve(count)

    async def remove_items(self, room_id: str, indices: list[int]) -> None:
        """Remove frames at indices WITHOUT shifting (set to None)."""
        io = self._get_io(room_id)
        for idx in indices:
            await io._backend.set(idx, None)  # type: ignore[union-attr]

    async def close(self) -> None:
        """Release in-memory handles without deleting stored data."""
        self._rooms.clear()
