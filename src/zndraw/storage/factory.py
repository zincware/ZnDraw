"""Storage backend factory for creating storage instances."""

import logging
import os

from .base import StorageBackend
from .asebytes_backend import ASEBytesStorageBackend

log = logging.getLogger(__name__)


def create_storage(
    room_id: str,
    base_path: str | None = None,
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
    **kwargs
        Additional backend-specific configuration (unused)

    Returns
    -------
    StorageBackend
        Configured ASEBytesStorageBackend instance

    Raises
    ------
    ValueError
        If base_path is None
    """
    if base_path is None:
        raise ValueError("ASEBytes backend requires a base_path (in-memory storage not supported)")

    # One LMDB file per room
    db_path = f"{base_path}/{room_id}.lmdb"
    os.makedirs(os.path.dirname(db_path) if os.path.dirname(db_path) else ".", exist_ok=True)

    log.info(f"Created ASEBytes storage at '{db_path}' for room '{room_id}'")
    return ASEBytesStorageBackend(db_path)
