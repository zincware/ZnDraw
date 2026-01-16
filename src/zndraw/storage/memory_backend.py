"""In-memory storage backend for single-instance deployments."""

import logging

import numpy as np

from .base import StorageBackend

log = logging.getLogger(__name__)


class InMemoryStorageBackend(StorageBackend):
    """In-memory storage backend using a Python list.

    This backend stores frames entirely in memory with:
    - Zero persistence (data lost on restart)
    - Fast access (no I/O overhead)
    - Simple implementation (no external dependencies)

    Ideal for:
    - Development and testing
    - Single-instance deployments with transient data
    - Scenarios where persistence is not required

    Data format:
    - Storage: list[dict[bytes, bytes]] (msgpack-serialized values)
    - Keys use dot notation: b"arrays.positions", b"info.energy"
    """

    def __init__(self):
        """Initialize in-memory storage backend."""
        self._frames: list[dict[bytes, bytes]] = []
        log.debug("Initialized InMemoryStorageBackend")

    def get(
        self,
        index: int | list[int] | slice | np.ndarray,
        keys: list[str] | None = None,
    ) -> dict[bytes, bytes] | list[dict[bytes, bytes]]:
        """Get frame(s) with optional key filtering.

        Parameters
        ----------
        index : int | list[int] | slice | np.ndarray
            Frame index/indices to retrieve
        keys : list[str] | None
            Optional list of keys to filter (e.g., ["arrays.positions", "info.energy"])

        Returns
        -------
        dict[bytes, bytes] | list[dict[bytes, bytes]]
            For single index: dict[bytes, bytes] with msgpack-serialized values
            For multiple indices: list of dicts
        """
        # Handle numpy arrays and scalars
        if isinstance(index, np.ndarray):
            if index.ndim == 0:
                index = int(index.item())
            else:
                index = index.tolist()
        elif isinstance(index, np.integer):
            index = int(index)

        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))

        is_single = isinstance(index, int)
        if is_single:
            index = [index]

        # Validate bounds
        length = len(self)
        for i in index:
            if i < -length or i >= length:
                raise IndexError(
                    f"Index {i} is out of bounds for storage of length {length}"
                )

        # Normalize negative indices
        normalized_indices = [i if i >= 0 else length + i for i in index]

        # Convert keys to bytes for filtering
        keys_bytes = None
        if keys is not None:
            keys_bytes = set(k.encode() for k in keys)

        # Build results
        results = []
        for i in normalized_indices:
            frame = self._frames[i]

            # Filter by keys if requested
            if keys_bytes is not None:
                frame = {k: v for k, v in frame.items() if k in keys_bytes}

            results.append(frame)

        if is_single:
            return results[0]

        return results

    def extend(self, values: list[dict[bytes, bytes]]) -> None:
        """Extend storage with multiple frames (batch write).

        Parameters
        ----------
        values : list[dict[bytes, bytes]]
            List of frame dictionaries with msgpack-serialized values
            Keys should be like b"arrays.positions", b"info.energy"
        """
        if not values:
            return

        self._frames.extend(values)
        log.debug(f"Extended storage with {len(values)} frames")

    def get_available_keys(self, index: int) -> list[str]:
        """List all keys available for a frame.

        Parameters
        ----------
        index : int
            Frame index

        Returns
        -------
        list[str]
            List of available keys (e.g., ["arrays.positions", "info.energy"])
        """
        length = len(self)
        if index < -length or index >= length:
            raise IndexError(
                f"Index {index} is out of bounds for storage of length {length}"
            )

        # Normalize negative index
        if index < 0:
            index = length + index

        return [k.decode() for k in self._frames[index].keys()]

    def __len__(self) -> int:
        """Return the number of frames in storage."""
        return len(self._frames)
