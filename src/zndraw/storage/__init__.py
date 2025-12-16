"""ZnDraw storage backend.

This package provides storage backends for trajectory data:
- InMemoryStorageBackend: In-memory storage (no persistence)
- ASEBytesStorageBackend: Local LMDB storage (default)
- MongoDBStorageBackend: Distributed MongoDB storage

Data is serialized with msgpack and sent directly over WebSockets.
"""

from .asebytes_backend import ASEBytesStorageBackend
from .base import StorageBackend
from .factory import create_storage
from .memory_backend import InMemoryStorageBackend
from .mongodb_backend import MongoDBStorageBackend

__all__ = [
    "StorageBackend",
    "InMemoryStorageBackend",
    "ASEBytesStorageBackend",
    "MongoDBStorageBackend",
    "create_storage",
]
