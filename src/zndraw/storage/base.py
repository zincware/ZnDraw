"""Base storage interface for ZnDraw trajectory storage backends."""

from abc import ABC, abstractmethod
from collections.abc import MutableSequence
import typing as t

import numpy as np


class StorageBackend(MutableSequence, ABC):
    """Abstract base class for storage backends.

    This interface defines the contract that all storage backends must implement
    to be compatible with ZnDraw's trajectory storage system.
    """

    @abstractmethod
    def get(
        self, index: int | list[int] | slice | np.ndarray, keys: list[str] | None = None
    ) -> dict:
        """Get frame(s) by index with optional key filtering.

        Parameters
        ----------
        index : int | list[int] | slice | np.ndarray
            Frame index/indices to retrieve
        keys : list[str] | None
            Optional list of keys to filter. If None, return all keys.

        Returns
        -------
        dict
            Frame data dictionary. For single index, returns single frame dict.
            For multiple indices, returns dict with concatenated arrays.
        """
        pass

    @abstractmethod
    def extend(self, values: list[dict]) -> None:
        """Extend storage with multiple frames (batch write).

        This is the primary write operation and should be optimized
        for bulk inserts.

        Parameters
        ----------
        values : list[dict]
            List of frame dictionaries to append
        """
        pass

    @abstractmethod
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
        pass

    @abstractmethod
    def __len__(self) -> int:
        """Return the number of frames in storage."""
        pass

    # MutableSequence required methods with default implementations
    def __getitem__(self, index: int | list[int] | slice | np.ndarray) -> dict:
        """Get frame(s) by index."""
        return self.get(index)

    def append(self, value: dict) -> None:
        """Append a single frame."""
        self.extend([value])

    def __setitem__(self, index: int | list[int] | slice, value: dict | list[dict]):
        """Not implemented - updates not supported."""
        raise NotImplementedError("Frame updates not supported")

    def __delitem__(self, index: int | list[int] | slice):
        """Not implemented - deletions not supported."""
        raise NotImplementedError("Frame deletions not supported")

    def insert(self, index: int, value: dict):
        """Not implemented - insertions not supported."""
        raise NotImplementedError("Frame insertions not supported")
