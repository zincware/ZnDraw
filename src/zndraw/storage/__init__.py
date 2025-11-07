"""ZnDraw storage backend.

This package provides ASE bytes-based storage for trajectory data using LMDB.
Data is serialized with msgpack and sent directly over WebSockets.
"""

from .base import StorageBackend
from .asebytes_backend import ASEBytesStorageBackend
from .factory import create_storage

__all__ = [
    "StorageBackend",
    "ASEBytesStorageBackend",
    "create_storage",
]
