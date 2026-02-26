"""In-memory storage backend using Python lists.

Stores frame data as dict[bytes, bytes] for efficient binary transfer.
Input frames (dict[str, Any]) are converted to raw bytes on storage.
"""

from collections import defaultdict
from typing import Any

from .base import RawFrame, StorageBackend, to_raw_frame


class InMemoryStorage(StorageBackend):
    """In-memory storage using Python lists.

    Stores frames as dict[bytes, bytes] for efficient binary transfer.
    Fast but not persistent - data is lost on server restart.
    Good for development and testing.

    Slots may be None for reserved-but-unfilled or removed entries.
    """

    def __init__(self) -> None:
        """Initialize empty storage."""
        self._data: dict[str, list[RawFrame | None]] = defaultdict(list)

    async def get(self, room_id: str, index: int) -> RawFrame | None:
        """Get a single frame by index.

        Returns None for out-of-bounds or reserved-but-unfilled slots.
        """
        frames = self._data[room_id]
        if 0 <= index < len(frames):
            return frames[index]
        return None

    async def get_range(
        self, room_id: str, start: int, stop: int | None
    ) -> list[RawFrame | None]:
        """Get a range of frames [start, stop). None for empty slots."""
        return list(self._data[room_id][start:stop])

    async def get_many(self, room_id: str, indices: list[int]) -> list[RawFrame | None]:
        """Get multiple frames by specific indices. None for missing/empty."""
        frames = self._data[room_id]
        length = len(frames)
        return [frames[i] if 0 <= i < length else None for i in indices]

    async def extend(
        self, room_id: str, frames: list[dict[str, Any]] | list[RawFrame]
    ) -> int:
        """Append frames to storage."""
        raw_frames = [to_raw_frame(f) for f in frames]
        self._data[room_id].extend(raw_frames)
        return len(self._data[room_id])

    async def set_item(
        self, room_id: str, index: int, frame: dict[str, Any] | RawFrame
    ) -> None:
        """Set a frame at a specific index.

        Raises IndexError if index is out of bounds.
        """
        frames = self._data[room_id]
        if index < 0 or index >= len(frames):
            raise IndexError(f"Frame index {index} out of bounds for room {room_id}")
        frames[index] = to_raw_frame(frame)

    async def merge_item(self, room_id: str, index: int, partial: RawFrame) -> None:
        """Merge partial frame data into existing frame at index."""
        frames = self._data[room_id]
        if index < 0 or index >= len(frames):
            raise IndexError(f"Frame index {index} out of bounds for room {room_id}")
        existing = frames[index]
        if existing is None:
            raise IndexError(f"Frame index {index} is empty for room {room_id}")
        existing.update(partial)

    async def get_length(self, room_id: str) -> int:
        """Get total frame count for a room."""
        return len(self._data[room_id])

    async def delete_range(self, room_id: str, start: int, stop: int) -> None:
        """Delete a range of frames [start, stop). Shifts subsequent frames."""
        del self._data[room_id][start:stop]

    async def clear(self, room_id: str) -> None:
        """Delete all frames for a room."""
        self._data[room_id] = []

    async def reserve(self, room_id: str, count: int) -> None:
        """Pre-allocate space for count frames. Grow only."""
        frames = self._data[room_id]
        if count > len(frames):
            frames.extend([None] * (count - len(frames)))

    async def remove_items(self, room_id: str, indices: list[int]) -> None:
        """Remove frames at indices without shifting."""
        frames = self._data[room_id]
        for idx in indices:
            if 0 <= idx < len(frames):
                frames[idx] = None

    async def close(self) -> None:
        """Clean up resources."""
        self._data.clear()
