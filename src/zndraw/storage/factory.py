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
    mongodb_url: str | None = None,
    mongodb_database: str | None = None,
) -> StorageBackend:
    """Factory for creating storage backends.

    If mongodb_url is provided, creates a MongoDB backend.
    Otherwise, creates an ASEBytes (LMDB) backend.

    Parameters
    ----------
    room_id : str
        Room identifier for storage isolation
    base_path : str | None
        Base directory path for LMDB storage.
        Required if not using MongoDB.
    map_size : int | None
        Maximum LMDB database size in bytes (virtual allocation).
        Required if not using MongoDB.
    mongodb_url : str | None
        MongoDB connection URI. If set, uses MongoDB backend.
    mongodb_database : str | None
        MongoDB database name. Required if mongodb_url is set.

    Returns
    -------
    StorageBackend
        Configured storage backend instance

    Raises
    ------
    ValueError
        If required parameters are missing for the selected backend
    """
    # MongoDB backend
    if mongodb_url is not None:
        from .mongodb_backend import MongoDBStorageBackend

        if mongodb_database is None:
            mongodb_database = "zndraw"

        log.info(
            f"Creating MongoDB storage: database='{mongodb_database}', "
            f"collection='{room_id}'"
        )
        return MongoDBStorageBackend(
            uri=mongodb_url,
            database=mongodb_database,
            room_id=room_id,
        )

    # LMDB backend (default)
    if base_path is None:
        raise ValueError(
            "ASEBytes backend requires a base_path (in-memory storage not supported)"
        )
    if map_size is None:
        raise ValueError("map_size must be provided from config")

    # One LMDB file per room
    db_path = f"{base_path}/{room_id}.lmdb"
    os.makedirs(
        os.path.dirname(db_path) if os.path.dirname(db_path) else ".", exist_ok=True
    )

    log.info(f"Created ASEBytes storage at '{db_path}' for room '{room_id}'")
    return ASEBytesStorageBackend(db_path, map_size=map_size)
