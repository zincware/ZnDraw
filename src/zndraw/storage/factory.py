"""Storage backend factory for creating storage instances."""

import logging
import os

from zndraw.config import LMDBStorageConfig, MongoDBStorageConfig, StorageConfig

from .asebytes_backend import ASEBytesStorageBackend
from .base import StorageBackend
from .mongodb_backend import MongoDBStorageBackend

log = logging.getLogger(__name__)


def create_storage(room_id: str, config: StorageConfig) -> StorageBackend:
    """Factory for creating storage backends based on config type.

    Parameters
    ----------
    room_id : str
        Room identifier for storage isolation
    config : StorageConfig
        Storage configuration (LMDBStorageConfig or MongoDBStorageConfig)

    Returns
    -------
    StorageBackend
        Configured storage backend instance
    """
    match config:
        case MongoDBStorageConfig():
            log.info(
                f"Creating MongoDB storage: database='{config.database}', "
                f"collection='{room_id}'"
            )
            return MongoDBStorageBackend(
                uri=config.url,
                database=config.database,
                room_id=room_id,
            )

        case LMDBStorageConfig():
            db_path = f"{config.path}/{room_id}.lmdb"
            os.makedirs(
                os.path.dirname(db_path) if os.path.dirname(db_path) else ".",
                exist_ok=True,
            )
            log.info(f"Creating LMDB storage at '{db_path}' for room '{room_id}'")
            return ASEBytesStorageBackend(db_path, map_size=config.map_size)
