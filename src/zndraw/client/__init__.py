"""ZnDraw Python client package."""

# Re-export accessors previously importable via
# `from zndraw.client import Sessions`
from zndraw.accessors import Sessions
from zndraw.client.api import APIManager
from zndraw.client.core import ZnDraw
from zndraw.client.exceptions import NotConnectedError
from zndraw.client.lock import ZnDrawLock
from zndraw.client.serialization import (
    _MAX_CHUNK_FRAMES,
    _TARGET_CHUNK_BYTES,
    _estimate_frame_size,
    atoms_to_json_dict,
    json_dict_to_atoms,
    raw_frame_to_atoms,
)
from zndraw.client.socket import SocketManager

# Re-export shared exceptions so `from zndraw.client import ZnDrawError` still works
from zndraw.exceptions import RoomLockedError, ZnDrawError

__all__ = [
    "_MAX_CHUNK_FRAMES",
    "_TARGET_CHUNK_BYTES",
    "APIManager",
    "NotConnectedError",
    "RoomLockedError",
    "Sessions",
    "SocketManager",
    "ZnDraw",
    "ZnDrawError",
    "ZnDrawLock",
    "_estimate_frame_size",
    "atoms_to_json_dict",
    "json_dict_to_atoms",
    "raw_frame_to_atoms",
]
