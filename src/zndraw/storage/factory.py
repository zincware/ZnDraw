"""Storage backend factory for creating storage instances."""

import logging
import os

from zndraw.config import (
    InMemoryStorageConfig,
    LMDBStorageConfig,
    MongoDBStorageConfig,
    StorageConfig,
)

from .asebytes_backend import ASEBytesStorageBackend
from .base import StorageBackend
from .memory_backend import InMemoryStorageBackend
from .mongodb_backend import MongoDBStorageBackend

log = logging.getLogger(__name__)


def create_storage(room_id: str, config: StorageConfig) -> StorageBackend:
    """Factory for creating storage backends based on config type.

    Parameters
    ----------
    room_id : str
        Room identifier for storage isolation
    config : StorageConfig
        Storage configuration (InMemoryStorageConfig, LMDBStorageConfig,
        or MongoDBStorageConfig)

    Returns
    -------
    StorageBackend
        Configured storage backend instance
    """
    match config:
        case InMemoryStorageConfig():
            log.debug(f"Creating in-memory storage for room '{room_id}'")
            return InMemoryStorageBackend()

        case MongoDBStorageConfig():
            log.debug(
                f"Creating MongoDB storage: uri='{config.get_masked_url()}', "
                f"database='{config.database}', collection='{room_id}'"
            )
            return MongoDBStorageBackend(
                uri=config.url,
                database=config.database,
                room_id=room_id,
            )

        case LMDBStorageConfig():
            db_path = os.path.join(config.path, f"{room_id}.lmdb")
            os.makedirs(config.path, exist_ok=True)
            log.debug(f"Creating LMDB storage at '{db_path}' for room '{room_id}'")
            return ASEBytesStorageBackend(db_path, map_size=config.map_size)

        case _:
            raise ValueError(
                f"Unsupported storage config type: {type(config).__name__}"
            )
