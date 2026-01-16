"""Shared utility functions for route handlers.

This module contains common functions used across multiple route blueprints,
following DRY principles and reducing code duplication.
"""

import functools
import json
import logging
import typing as t

from flask import current_app, jsonify, request

from zndraw.auth import AuthError, get_current_user
from zndraw.server import socketio
from zndraw.storage import StorageBackend, create_storage

from .constants import SocketEvents
from .redis_keys import RoomKeys, SessionKeys

log = logging.getLogger(__name__)

# Shared storage dictionary for storage backend instances
STORAGE: dict[str, StorageBackend] = {}


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
    return RoomKeys(room).lock(target)


def get_storage(room_id: str) -> StorageBackend:
    """Get or create a storage backend for a room.

    Parameters
    ----------
    room_id : str
        Room identifier

    Returns
    -------
    StorageBackend
        Storage backend instance for the room
    """
    if room_id not in STORAGE:
        config = current_app.extensions["config"]
        STORAGE[room_id] = create_storage(room_id=room_id, config=config.storage)

    return STORAGE[room_id]


LockTargetType = str | t.Callable[[t.Any, dict], str] | None


def check_lock(
    target: LockTargetType = None,
    forbid: list[str] | None = None,
):
    """Check for blocking locks before operation. Does NOT require holding lock.

    This decorator checks if an intentional lock blocks the operation, but does
    NOT require the caller to hold the lock. If no lock exists, the operation
    proceeds (FIFO ordering handles concurrency). Only blocks if ANOTHER session
    holds a lock on the target.

    Parameters
    ----------
    target : str | Callable | None, optional
        Lock target to check. Can be:
        - str: Static target like "trajectory:meta", "bookmarks"
        - Callable: Function(request, kwargs) -> str for dynamic targets
          like "geometry:{key}" where key comes from request
        - None: No target lock check (only forbid + global lock checks)

        If a lock exists for this target held by a DIFFERENT session -> 423 Locked
        If no lock exists OR same session holds it -> proceed
        Default: None

    forbid : list[str] | None, optional
        List of permanent lock states that block the operation.
        Examples: ["room:locked"]
        These return 403 Forbidden (not 423 Locked).
        Default: None

    Execution Order
    ---------------
    1. Check forbidden permanent states (e.g., room:locked) -> 403 if violated
    2. Authenticate user via JWT -> 401 if invalid
    3. Validate session ID -> 400 if missing, 401 if invalid
    4. Check global lock -> 423 if held by other session
    5. Check target lock -> 423 if held by other session
    6. Execute route handler (FIFO handles ordering)

    Returns
    -------
    function
        Decorated function that checks locks before execution

    Examples
    --------
    >>> # Proceed unless blocked by trajectory lock or room locked
    >>> @check_lock(target="trajectory:meta", forbid=["room:locked"])
    >>> def append_frame(room_id: str, session_id: str, user_id: str):
    ...     pass

    >>> # Per-geometry lock check with dynamic target
    >>> @check_lock(
    ...     target=lambda req, kw: f"geometry:{req.json.get('key')}",
    ...     forbid=["room:locked"]
    ... )
    >>> def create_geometry(room_id: str, session_id: str, user_id: str):
    ...     pass

    >>> # Only check room lock, no specific target
    >>> @check_lock(forbid=["room:locked"])
    >>> def set_bookmark(room_id: str, session_id: str, user_id: str):
    ...     pass
    """

    def decorator(f):
        @functools.wraps(f)
        def wrapped(*args, **kwargs):
            room_id = kwargs.get("room_id")
            if not room_id:
                return jsonify({"error": "room_id required"}), 400

            r = current_app.extensions["redis"]
            keys = RoomKeys(room_id)

            # ===== 1. CHECK FORBIDDEN PERMANENT STATES -> 403 =====
            if forbid:
                for forbidden_state in forbid:
                    if forbidden_state == "room:locked":
                        locked = r.get(keys.locked())
                        if locked == "1":
                            locked_by = r.get(keys.locked_by())
                            error_msg = "Room is locked and cannot be modified"
                            if locked_by:
                                error_msg = (
                                    f"Room is locked by admin '{locked_by}' "
                                    "and cannot be modified"
                                )
                            log.warning(
                                f"Forbidden state violated: {forbidden_state} "
                                f"in room {room_id}"
                            )
                            return jsonify({"error": error_msg}), 403

            # ===== 2. AUTHENTICATE USER VIA JWT -> 401 =====
            try:
                user_id = get_current_user()
            except AuthError as e:
                log.warning(f"JWT auth failed: {e.message}")
                return jsonify({"error": e.message}), e.status_code

            # ===== 3. VALIDATE SESSION ID -> 400/401 =====
            session_id = request.headers.get("X-Session-ID")
            if not session_id:
                log.warning(f"X-Session-ID header missing in room {room_id}")
                return jsonify({"error": "X-Session-ID header required"}), 400

            session_key = SessionKeys.session_data(session_id)
            session_data_str = r.get(session_key)

            if not session_data_str:
                return jsonify({"error": "Invalid or expired session"}), 401

            try:
                session_data = json.loads(session_data_str)
                session_user = session_data.get("userId")
            except json.JSONDecodeError:
                return jsonify({"error": "Invalid session data"}), 500

            if session_user != user_id:
                return jsonify({"error": "Session/user mismatch"}), 403

            # ===== 4. CHECK GLOBAL LOCK -> 423 if held by other session =====
            global_lock_str = r.get(keys.lock("global"))
            if global_lock_str:
                try:
                    global_lock_data = json.loads(global_lock_str)
                    if global_lock_data.get("sessionId") != session_id:
                        holder = global_lock_data.get("userId", "unknown")
                        log.debug(
                            f"Global lock blocks operation in room {room_id} "
                            f"(held by {holder})"
                        )
                        return (
                            jsonify(
                                {"error": "Room is globally locked", "holder": holder}
                            ),
                            423,
                        )
                except json.JSONDecodeError:
                    log.warning(f"Invalid global lock data in room {room_id}")

            # ===== 5. CHECK TARGET LOCK -> 423 if held by other session =====
            if target is not None:
                # Resolve dynamic target if callable
                if callable(target):
                    try:
                        actual_target = target(request, kwargs)
                    except (AttributeError, TypeError, ValueError, KeyError) as e:
                        log.warning(
                            f"Failed to resolve dynamic lock target in room "
                            f"{room_id}: {e}"
                        )
                        return (
                            jsonify(
                                {
                                    "error": f"Invalid or missing request body for "
                                    f"dynamic target: {e}"
                                }
                            ),
                            400,
                        )
                else:
                    actual_target = target

                lock_key = keys.lock(actual_target)
                lock_data_str = r.get(lock_key)

                if lock_data_str:
                    try:
                        lock_data = json.loads(lock_data_str)
                        if lock_data.get("sessionId") != session_id:
                            holder = lock_data.get("userId", "unknown")
                            log.debug(
                                f"Lock '{actual_target}' blocks operation in room "
                                f"{room_id} (held by {holder})"
                            )
                            return (
                                jsonify(
                                    {
                                        "error": f"Resource locked by {holder}",
                                        "target": actual_target,
                                        "holder": holder,
                                    }
                                ),
                                423,
                            )
                    except json.JSONDecodeError:
                        log.warning(
                            f"Invalid lock data for '{actual_target}' in room {room_id}"
                        )

            # ===== 6. INJECT VALIDATED PARAMETERS AND PROCEED =====
            kwargs["session_id"] = session_id
            kwargs["user_id"] = user_id

            return f(*args, **kwargs)

        return wrapped

    return decorator


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
    room_keys = RoomKeys(room_id)
    bookmarks_key = room_keys.bookmarks()
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
    room_keys = RoomKeys(room_id)
    bookmarks_key = room_keys.bookmarks()
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
    room_keys = RoomKeys(room_id)
    r.hdel(room_keys.bookmarks(), str(index))


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
    affected_keys: list[str] | None = None,
):
    """Emit frames:invalidate event to tell clients to clear their frame cache.

    Parameters
    ----------
    room_id : str
        The room identifier
    operation : str
        Type of operation - 'delete', 'insert', 'replace', or 'bulk_replace'
    affected_index : int | None, optional
        For 'replace' - the specific frame index that was replaced
    affected_from : int | None, optional
        For 'delete'/'insert'/'bulk_replace' - all frames from this index onward are affected
    affected_keys : list[str] | None, optional
        Optional list of frame data keys that changed (e.g., ["arrays.position", "arrays.radii"]).
        If None, invalidates ALL keys for affected frames (backward compatible behavior).
        If empty list, no frame data changed (only metadata/structure).
    """
    data = {
        "roomId": room_id,
        "operation": operation,
    }

    if affected_index is not None:
        data["affectedIndex"] = affected_index
    if affected_from is not None:
        data["affectedFrom"] = affected_from
    if affected_keys is not None:
        data["affectedKeys"] = affected_keys

    log.debug(f"Emitting frames:invalidate for room '{room_id}': {data}")

    socketio.emit(
        "frames:invalidate",
        data,
        to=f"room:{room_id}",
    )


def emit_lock_update(
    room_id: str,
    target: str,
    action: str,
    user_name: str | None = None,
    message: str | None = None,
    timestamp: str | None = None,
    session_id: str | None = None,
):
    """Emit lock status update to all clients in room.

    Parameters
    ----------
    room_id : str
        The room identifier
    target : str
        Lock target (e.g., "trajectory:meta")
    action : str
        Lock action - "acquired", "released", or "refreshed"
    user_name : str | None, optional
        Username of lock holder (None for released)
    message : str | None, optional
        Lock message describing the operation
    timestamp : str | None, optional
        ISO 8601 timestamp
    session_id : str | None, optional
        Session ID of the client that initiated the lock action.
        Clients can use this to filter out their own events.
    """
    payload = {
        "roomId": room_id,
        "target": target,
        "action": action,
        "holder": user_name,
        "message": message,
        "timestamp": timestamp,
        "sessionId": session_id,
    }

    socketio.emit(
        "lock:update",
        payload,
        to=f"room:{room_id}",
        skip_sid=None,  # Send to all clients; they filter by sessionId
    )

    log.debug(f"Emitted lock:update for room '{room_id}': {payload}")


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


def emit_step_update_on_delete(room_id: str, deleted_indices: list[int]):
    """Adjust and broadcast step after frames are deleted.

    When frames are deleted before the current step, the step needs to be
    decremented to point to the same logical frame data.

    Parameters
    ----------
    room_id : str
        The room ID
    deleted_indices : list[int]
        List of deleted frame indices (not necessarily sorted)
    """
    if not deleted_indices:
        return

    r = current_app.extensions["redis"]
    room_keys = RoomKeys(room_id)

    # Get current step
    current_step_raw = r.get(room_keys.current_frame())
    if current_step_raw is None:
        return

    current_step = int(current_step_raw)

    # Calculate how many frames before current step were deleted
    shift = sum(1 for idx in deleted_indices if idx < current_step)

    if shift == 0:
        return

    # Calculate new step
    new_step = current_step - shift

    # Ensure step is non-negative
    new_step = max(0, new_step)

    # Update Redis
    r.set(room_keys.current_frame(), new_step)

    # Emit frame update to all clients
    socketio.emit(
        SocketEvents.FRAME_UPDATE,
        {"frame": new_step},
        to=f"room:{room_id}",
    )

    log.debug(
        f"Step adjusted after delete: {current_step} -> {new_step} "
        f"(shift: -{shift}) in room '{room_id}'"
    )


def emit_step_update_on_insert(room_id: str, insert_position: int):
    """Adjust and broadcast step after a frame is inserted.

    When a frame is inserted at or before the current step, the step needs
    to be incremented to point to the same logical frame data.

    Parameters
    ----------
    room_id : str
        The room ID
    insert_position : int
        The index where the frame was inserted
    """
    r = current_app.extensions["redis"]
    room_keys = RoomKeys(room_id)

    # Get current step
    current_step_raw = r.get(room_keys.current_frame())
    if current_step_raw is None:
        return

    current_step = int(current_step_raw)

    # Only shift if frame was inserted at or before current step
    if insert_position > current_step:
        return

    # Calculate new step
    new_step = current_step + 1

    # Update Redis
    r.set(room_keys.current_frame(), new_step)

    # Emit frame update to all clients
    socketio.emit(
        SocketEvents.FRAME_UPDATE,
        {"frame": new_step},
        to=f"room:{room_id}",
    )

    log.debug(
        f"Step adjusted after insert: {current_step} -> {new_step} "
        f"(insert at {insert_position}) in room '{room_id}'"
    )


def parse_frame_mapping(
    mapping_entry: str | bytes, default_room: str
) -> tuple[str, int]:
    """Parse a frame mapping entry to extract source room and physical index.

    Parameters
    ----------
    mapping_entry : str | bytes
        Mapping entry in format "room_id:physical_index" or just "physical_index"
    default_room : str
        Default room ID to use if mapping entry doesn't contain room prefix

    Returns
    -------
    tuple[str, int]
        Tuple of (source_room_id, physical_index)
    """
    # Decode bytes to string if necessary
    if isinstance(mapping_entry, bytes):
        mapping_entry = mapping_entry.decode("utf-8")

    # Parse mapping entry
    if ":" in mapping_entry:
        source_room_id, physical_index_str = mapping_entry.split(":", 1)
        return source_room_id, int(physical_index_str)
    else:
        return default_room, int(mapping_entry)


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
    room_keys = RoomKeys(room_id)
    metadata_lock_key = get_lock_key(room_id, "trajectory:meta")

    if not r.exists(metadata_lock_key):
        return None

    # Lock exists - try to get metadata
    metadata_raw = r.get(room_keys.lock_metadata("trajectory:meta"))
    if metadata_raw:
        metadata = json.loads(metadata_raw)
        lock_metadata = LockMetadata(
            msg=metadata.get("msg"),
            userName=metadata.get("userName"),
            timestamp=metadata.get("timestamp"),
        )
        return lock_metadata.model_dump()

    # Lock exists but no metadata - get basic info from lock data
    lock_data_str = r.get(metadata_lock_key)
    if lock_data_str:
        lock_data = json.loads(lock_data_str)
        lock_holder_user_name = lock_data.get("userId", "unknown")
    else:
        lock_holder_user_name = "unknown"

    lock_metadata = LockMetadata(
        msg=None, userName=lock_holder_user_name, timestamp=None
    )
    return lock_metadata.model_dump()
