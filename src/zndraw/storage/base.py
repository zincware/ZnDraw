"""Abstract base class for frame storage backends."""

import base64
from abc import ABC, abstractmethod
from typing import Any

import msgpack

# Type alias for raw frame data: dict with bytes keys and bytes values
# Keys are like b"arrays.positions", b"cell", etc.
# Values are msgpack-numpy encoded bytes
RawFrame = dict[bytes, bytes]


def to_raw_frame(frame: dict[str, Any] | RawFrame) -> RawFrame:
    """Convert input frame dict to raw bytes format.

    Handles two input formats:
    1. Already raw: dict[bytes, bytes] — pass through
    2. Base64 encoded: dict[str, str] with keys like "b64:..." — decode

    Parameters
    ----------
    frame
        Input frame data (dict[str, Any] or dict[bytes, bytes]).
    """
    if not frame:
        return {}

    first_key = next(iter(frame.keys()))

    # Already raw bytes format
    if isinstance(first_key, bytes):
        return frame  # type: ignore[return-value]

    # Base64 encoded format (from JSON input)
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


class StorageBackend(ABC):
    """Abstract base class for frame storage backends.

    All storage backends must implement these async methods.
    Frames are stored as dict[bytes, bytes] where:
    - Keys are raw bytes (e.g., b"arrays.positions", b"cell")
    - Values are msgpack-numpy encoded bytes

    This raw bytes format matches the reference implementation and enables
    efficient binary transfer without base64 encoding overhead.

    The storage is organized by room_id, with each room having an ordered
    sequence of frames indexed from 0. Implementations must ensure thread-safety
    and handle concurrent access appropriately.
    """

    @abstractmethod
    async def get(self, room_id: str, index: int) -> RawFrame | None:
        """Get a single frame by index.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        index
            The frame index (0-based).

        Returns
        -------
        RawFrame | None
            The frame data, or None if not found.
        """

    @abstractmethod
    async def get_range(
        self, room_id: str, start: int, stop: int | None
    ) -> list[RawFrame | None]:
        """Get a range of frames [start, stop).

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
        list[RawFrame | None]
            Positional list. ``None`` for slots without data (virtual frames).
        """

    @abstractmethod
    async def get_many(self, room_id: str, indices: list[int]) -> list[RawFrame | None]:
        """Get multiple frames by specific indices.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        indices
            Frame indices to retrieve.

        Returns
        -------
        list[RawFrame | None]
            Positional list: ``result[i]`` corresponds to ``indices[i]``.
            ``None`` for out-of-bounds or empty slots.
        """

    @abstractmethod
    async def extend(
        self, room_id: str, frames: list[dict[str, Any]] | list[RawFrame]
    ) -> int:
        """Append frames to storage.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        frames
            Frame dicts to append. Accepts ``dict[str, Any]`` (JSON)
            or ``dict[bytes, bytes]`` (raw). Converted internally.

        Returns
        -------
        int
            New total frame count for the room.
        """

    @abstractmethod
    async def set_item(
        self, room_id: str, index: int, frame: dict[str, Any] | RawFrame
    ) -> None:
        """Set a frame at a specific index.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        index
            The frame index to set.
        frame
            The frame data to store.

        Raises
        ------
        IndexError
            If index is out of bounds.
        """

    @abstractmethod
    async def merge_item(self, room_id: str, index: int, partial: RawFrame) -> None:
        """Merge partial frame data into existing frame at index.

        Reads the existing frame, merges the partial data (last-write-wins),
        and writes back atomically.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        index
            The frame index to merge into.
        partial
            Partial frame data to merge (dict[bytes, bytes]).

        Raises
        ------
        IndexError
            If index is out of bounds.
        """

    @abstractmethod
    async def get_length(self, room_id: str) -> int:
        """Get total frame count for a room.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        """

    @abstractmethod
    async def delete_range(self, room_id: str, start: int, stop: int) -> None:
        """Delete a range of frames [start, stop).

        Frames after the deleted range are shifted to fill the gap.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        start
            Start index (inclusive).
        stop
            Stop index (exclusive).
        """

    @abstractmethod
    async def clear(self, room_id: str) -> None:
        """Delete all frames for a room.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        """

    @abstractmethod
    async def reserve(self, room_id: str, count: int) -> None:
        """Pre-allocate space for count frames. Idempotent for growth.

        After reserve(room_id, N), set_item(room_id, i, frame) must succeed
        for any 0 <= i < N without IndexError. get() returns None for
        reserved-but-unfilled slots.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        count
            Target capacity. Only grows — if current length >= count, no-op.
        """

    @abstractmethod
    async def remove_items(self, room_id: str, indices: list[int]) -> None:
        """Remove frames at indices WITHOUT shifting.

        After removal, get() returns None for those indices.
        Does not change get_length().

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        indices
            Frame indices to remove.
        """

    async def set_frame_count(self, room_id: str, count: int) -> None:
        """Set external frame count (e.g. from a provider mount).

        Only supported by backends with external metadata storage (e.g.
        ``StorageRouter`` backed by Redis).  Default raises
        ``NotImplementedError``.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        count
            The frame count to set.
        """
        raise NotImplementedError(
            "This storage backend does not support external frame counts"
        )

    async def clear_frame_count(self, room_id: str) -> None:
        """Clear external frame count.

        Parameters
        ----------
        room_id
            The room identifier (UUID string).
        """
        raise NotImplementedError(
            "This storage backend does not support external frame counts"
        )

    @abstractmethod
    async def close(self) -> None:
        """Clean up resources (file handles, connections, etc.)."""
