"""ASE bytes-based storage backend using LMDB."""

import logging
import typing as t

import numpy as np
from asebytes import BytesIO

from .base import StorageBackend

log = logging.getLogger(__name__)


class ASEBytesStorageBackend(StorageBackend):
    """ASE bytes-based storage using LMDB via asebytes.

    This backend uses asebytes.BytesIO for efficient storage with:
    - Memory-mapped I/O (zero-copy reads)
    - ACID transactions (data integrity)
    - msgpack serialization (fast, handles numpy)
    - Partial key retrieval via get(index, keys=None)
    - Full MutableSequence support (delete, insert work!)

    Data format:
    - Storage: dict[bytes, bytes] (msgpack-serialized values)
    - Network: dict[bytes, bytes] (sent directly via WebSocket)
    - Keys use dot notation: b"arrays.positions", b"info.energy"
    """

    def __init__(self, db_path: str, map_size: int):
        """Initialize ASE bytes storage backend.

        Parameters
        ----------
        db_path : str
            Path to LMDB database file
        map_size : int
            Maximum size of the database in bytes (virtual allocation)
        """
        self.db_path = db_path
        # No prefix needed since we use one database per room
        self.io = BytesIO(db_path, prefix=b"", map_size=map_size)
        log.info(f"Initialized ASEBytesStorageBackend at '{db_path}' (map_size={map_size / 1024**3:.2f} GB)")

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

        # Convert keys to bytes for BytesIO
        keys_bytes = None
        if keys is not None:
            keys_bytes = [k.encode() for k in keys]

        # Get data from BytesIO
        results = []
        for i in index:
            # BytesIO.get(index, keys) returns dict[bytes, bytes]
            msgpack_dict = self.io.get(i, keys=keys_bytes)
            results.append(msgpack_dict)

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

        # BytesIO.extend expects list[dict[bytes, bytes]]
        self.io.extend(values)

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
        keys_bytes = self.io.get_available_keys(index)
        return [k.decode() for k in keys_bytes]

    def __len__(self) -> int:
        """Return the number of frames in storage."""
        return len(self.io)
