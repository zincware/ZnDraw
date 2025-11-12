"""Storage backend factory for creating storage instances."""

import logging
import os

from .base import StorageBackend
from .asebytes_backend import ASEBytesStorageBackend

log = logging.getLogger(__name__)


def create_storage(
    room_id: str,
    base_path: str | None = None,
    map_size: int | None = None,
    **kwargs,
) -> StorageBackend:
    """Factory for creating storage backends.

    Parameters
    ----------
    room_id : str
        Room identifier for storage isolation
    base_path : str | None
        Base directory path for file-based storage.
        If None, raises ValueError (in-memory storage not supported).
    map_size : int | None
        Maximum LMDB database size in bytes (virtual allocation).
        Must be provided from config.
    **kwargs
        Additional backend-specific configuration (unused)

    Returns
    -------
    StorageBackend
        Configured ASEBytesStorageBackend instance

    Raises
    ------
    ValueError
        If base_path is None or map_size is None
    """
    if base_path is None:
        raise ValueError("ASEBytes backend requires a base_path (in-memory storage not supported)")
    if map_size is None:
        raise ValueError("map_size must be provided from config")

    # One LMDB file per room
    db_path = f"{base_path}/{room_id}.lmdb"
    os.makedirs(os.path.dirname(db_path) if os.path.dirname(db_path) else ".", exist_ok=True)

    log.info(f"Created ASEBytes storage at '{db_path}' for room '{room_id}'")
    return ASEBytesStorageBackend(db_path, map_size=map_size)
