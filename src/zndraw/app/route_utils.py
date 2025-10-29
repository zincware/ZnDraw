"""Shared utility functions for route handlers.

This module contains common functions used across multiple route blueprints,
following DRY principles and reducing code duplication.
"""

import json
import logging

from flask import current_app
from zarr.storage import MemoryStore
import zarr

from zndraw.server import socketio
from zndraw.storage import ZarrStorageSequence

from .constants import SocketEvents

log = logging.getLogger(__name__)

# Shared storage dictionary for Zarr stores
STORAGE: dict[str, MemoryStore] = {}


def get_lock_key(room: str, target: str) -> str:
    """Construct a standardized Redis key for a lock.

    Parameters
    ----------
    room : str
        Room identifier
    target : str
        Lock target (e.g., 'trajectory:meta')

    Returns
    -------
    str
        Redis key for the lock
    """
    return f"room:{room}:lock:{target}"


def get_zarr_store_path(room_id: str) -> str:
    """Return the path to the Zarr store for a given room.

    Parameters
    ----------
    room_id : str
        Room identifier

    Returns
    -------
    str
        Path to the Zarr store
    """
    storage_path = current_app.config.get("STORAGE_PATH", "./zndraw-data.zarr")
    # Remove .zarr extension if present to append room_id
    base_path = storage_path.rstrip("/").removesuffix(".zarr")
    return f"{base_path}/{room_id}.zarr"


def get_storage(room_id: str) -> ZarrStorageSequence:
    """Get or create a Zarr storage sequence for a room.

    Parameters
    ----------
    room_id : str
        Room identifier

    Returns
    -------
    ZarrStorageSequence
        Zarr storage sequence for the room
    """
    if room_id not in STORAGE:
        STORAGE[room_id] = MemoryStore()
    root = zarr.group(STORAGE[room_id])
    storage = ZarrStorageSequence(root)
    return storage


def check_room_locked(
    room_id: str, client_id: str | None = None
) -> tuple[dict[str, str], int] | None:
    """Check if a room is locked.

    This checks BOTH:
    1. Permanent room lock (room:{room_id}:locked)
    2. Temporary trajectory lock (trajectory:meta) - used by vis.lock context manager

    Parameters
    ----------
    room_id : str
        The room ID to check
    client_id : str | None, optional
        Optional client ID - if provided and this client holds the trajectory lock,
        the check will pass

    Returns
    -------
    tuple[dict[str, str], int] | None
        Error tuple if locked and client doesn't hold lock, None otherwise
    """
    redis_client = current_app.extensions["redis"]

    # Check permanent room lock
    locked = redis_client.get(f"room:{room_id}:locked")
    if locked == "1":
        return {"error": "Room is locked and cannot be modified"}, 403

    # Check trajectory lock (from vis.lock context manager)
    # This protects ALL trajectory operations: append, delete, upload, etc.
    trajectory_lock_key = get_lock_key(room_id, "trajectory:meta")
    lock_holder = redis_client.get(trajectory_lock_key)

    if lock_holder:
        log.debug(
            f"Lock check: room={room_id}, lock_holder={lock_holder}, client_id={client_id}"
        )
        # If a client_id is provided and it holds the lock, allow the operation
        if client_id and lock_holder == client_id:
            return None
        # Otherwise, the room is locked by another operation
        log.warning(
            f"Lock rejected: lock_holder={lock_holder} != client_id={client_id}"
        )
        return {"error": "Room is temporarily locked by another operation"}, 423

    return None


def shift_bookmarks_on_delete(room_id: str, deleted_indices: list[int]):
    """Shift bookmark indices after frames are deleted.

    Bookmarks at deleted indices are removed, bookmarks after are shifted down.

    Parameters
    ----------
    room_id : str
        The room ID
    deleted_indices : list[int]
        Sorted list of deleted frame indices
    """
    if not deleted_indices:
        return

    r = current_app.extensions["redis"]
    bookmarks_key = f"room:{room_id}:bookmarks"
    bookmarks_raw = r.hgetall(bookmarks_key)

    if not bookmarks_raw:
        return

    # Convert to int keys and sort deletion indices
    bookmarks = {int(k): v for k, v in bookmarks_raw.items()}
    deleted_set = set(deleted_indices)

    # Build new bookmarks dict with shifted indices
    new_bookmarks = {}
    for idx, label in bookmarks.items():
        # Remove bookmarks at deleted indices
        if idx in deleted_set:
            continue

        # Calculate how many indices below this one were deleted
        shift = sum(1 for del_idx in deleted_indices if del_idx < idx)
        new_idx = idx - shift
        new_bookmarks[new_idx] = label

    # Update Redis with shifted bookmarks
    r.delete(bookmarks_key)
    if new_bookmarks:
        r.hset(bookmarks_key, mapping={str(k): v for k, v in new_bookmarks.items()})


def shift_bookmarks_on_insert(room_id: str, insert_position: int):
    """Shift bookmark indices after a frame is inserted.

    All bookmarks at or after the insert position shift up by 1.

    Parameters
    ----------
    room_id : str
        The room ID
    insert_position : int
        The index where the frame was inserted
    """
    r = current_app.extensions["redis"]
    bookmarks_key = f"room:{room_id}:bookmarks"
    bookmarks_raw = r.hgetall(bookmarks_key)

    if not bookmarks_raw:
        return

    # Convert to int keys
    bookmarks = {int(k): v for k, v in bookmarks_raw.items()}

    # Build new bookmarks dict with shifted indices
    new_bookmarks = {}
    for idx, label in bookmarks.items():
        if idx >= insert_position:
            new_bookmarks[idx + 1] = label
        else:
            new_bookmarks[idx] = label

    # Update Redis with shifted bookmarks
    r.delete(bookmarks_key)
    if new_bookmarks:
        r.hset(bookmarks_key, mapping={str(k): v for k, v in new_bookmarks.items()})


def remove_bookmark_at_index(room_id: str, index: int):
    """Remove bookmark at a specific index (used when frame is replaced).

    Parameters
    ----------
    room_id : str
        The room ID
    index : int
        The frame index
    """
    r = current_app.extensions["redis"]
    bookmarks_key = f"room:{room_id}:bookmarks"
    r.hdel(bookmarks_key, str(index))


def emit_bookmarks_invalidate(room_id: str):
    """Emit bookmarks invalidate event to all clients after frame mapping changes.

    Clients will refetch bookmarks from the server when they receive this event.

    Parameters
    ----------
    room_id : str
        The room ID
    """
    socketio.emit(
        SocketEvents.INVALIDATE_BOOKMARK,
        {"operation": "frame_mapping_changed"},
        to=f"room:{room_id}",
    )


def emit_frames_invalidate(
    room_id: str,
    operation: str,
    affected_index: int | None = None,
    affected_from: int | None = None,
):
    """Emit frames:invalidate event to tell clients to clear their frame cache.

    Parameters
    ----------
    room_id : str
        The room identifier
    operation : str
        Type of operation - 'delete', 'insert', or 'replace'
    affected_index : int | None, optional
        For 'replace' - the specific frame index that was replaced
    affected_from : int | None, optional
        For 'delete'/'insert' - all frames from this index onward are affected
    """
    data = {
        "roomId": room_id,
        "operation": operation,
    }

    if affected_index is not None:
        data["affectedIndex"] = affected_index
    if affected_from is not None:
        data["affectedFrom"] = affected_from

    log.info(f"Emitting frames:invalidate for room '{room_id}': {data}")

    socketio.emit(
        "frames:invalidate",
        data,
        to=f"room:{room_id}",
    )


def emit_len_frames_update(room_id: str):
    """Emit len_frames event to update clients about the frame count after mutations.

    Parameters
    ----------
    room_id : str
        The room ID
    """
    room_service = current_app.extensions["room_service"]
    frame_count = room_service.get_frame_count(room_id)

    # Broadcast frame count update via room:update event
    from zndraw.app.room_manager import emit_room_update

    emit_room_update(socketio, room_id, frameCount=frame_count)


def get_metadata_lock_info(room_id: str) -> dict | None:
    """Get metadata lock information for a room.

    Checks if a metadata lock (trajectory:meta) is held for the room and returns
    lock metadata if available.

    Parameters
    ----------
    room_id : str
        The room ID

    Returns
    -------
    dict | None
        Lock metadata dictionary with 'msg', 'userName', 'timestamp' fields,
        or None if no lock is held
    """
    from .models import LockMetadata

    r = current_app.extensions["redis"]
    metadata_lock_key = get_lock_key(room_id, "trajectory:meta")

    if not r.exists(metadata_lock_key):
        return None

    # Lock exists - try to get metadata
    metadata_raw = r.get(f"{metadata_lock_key}:metadata")
    if metadata_raw:
        metadata = json.loads(metadata_raw)
        lock_metadata = LockMetadata(
            msg=metadata.get("msg"),
            userName=metadata.get("userName"),
            timestamp=metadata.get("timestamp"),
        )
        return lock_metadata.model_dump()

    # Lock exists but no metadata - get basic info from client
    lock_holder_client_id = r.get(metadata_lock_key)
    user_name = (
        r.hget(f"client:{lock_holder_client_id}", "userName")
        if lock_holder_client_id
        else None
    )
    lock_metadata = LockMetadata(
        msg=None, userName=user_name or "unknown", timestamp=None
    )
    return lock_metadata.model_dump()
