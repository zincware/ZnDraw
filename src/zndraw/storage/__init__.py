"""ZnDraw storage backend.

This package provides storage backends for trajectory data:
- ASEBytesStorageBackend: Local LMDB storage (default)
- MongoDBStorageBackend: Distributed MongoDB storage

Data is serialized with msgpack and sent directly over WebSockets.
"""

from .base import StorageBackend
from .asebytes_backend import ASEBytesStorageBackend
from .mongodb_backend import MongoDBStorageBackend
from .factory import create_storage

__all__ = [
    "StorageBackend",
    "ASEBytesStorageBackend",
    "MongoDBStorageBackend",
    "create_storage",
]
