"""REST API routes for distributed lock management."""

import json
import logging
import uuid

from flask import Blueprint, current_app, jsonify, request

from zndraw.auth import AuthError, get_current_user
from zndraw.utils.time import utc_now_iso

from .constants import LockConfig
from .redis_keys import SessionKeys
from .route_utils import emit_lock_update, get_lock_key

log = logging.getLogger(__name__)

locks = Blueprint("locks", __name__)


@locks.route("/api/rooms/<room_id>/locks/<target>/acquire", methods=["POST"])
def acquire_lock(room_id, target):
    """Acquire a lock for a specific target in a room.

    Parameters
    ----------
    room_id : str
        Room identifier
    target : str
        Lock target (e.g., "trajectory:meta")

    Headers
    -------
    X-Session-ID : str (required)
        Session ID from /join response

    Request Body
    ------------
    msg : str, optional
        Optional message describing the lock purpose

    Returns
    -------
    dict
        {"success": true, "lockToken": "...", "ttl": 60, "refreshInterval": 30} on success
        {"success": false, "error": "..."} on failure (423 if locked)
    """
    try:
        user_name = get_current_user()
    except AuthError as e:
        return jsonify({"success": False, "error": str(e)}), 401

    # Extract session ID from header
    session_id = request.headers.get("X-Session-ID")
    if not session_id:
        return jsonify({"success": False, "error": "Session ID required"}), 400

    from .redis_keys import SessionKeys

    r = current_app.extensions["redis"]
    msg = request.json.get("msg") if request.json else None

    # Validate session exists and matches user
    session_key = SessionKeys.session_data(session_id)
    session_data_str = r.get(session_key)
    if not session_data_str:
        return jsonify({"success": False, "error": "Invalid or expired session"}), 401

    try:
        session_data = json.loads(session_data_str)
        session_user = session_data.get("userId")
    except json.JSONDecodeError:
        return jsonify({"success": False, "error": "Invalid session data"}), 500

    # Verify session belongs to authenticated user
    if session_user != user_name:
        return jsonify({"success": False, "error": "Session/user mismatch"}), 403

    # Check global lock - block acquiring any other lock if global lock is held by another session
    if target != "global":
        global_lock_key = get_lock_key(room_id, "global")
        global_lock_str = r.get(global_lock_key)
        if global_lock_str:
            try:
                global_lock_data = json.loads(global_lock_str)
                if global_lock_data.get("sessionId") != session_id:
                    holder = global_lock_data.get("userId", "unknown")
                    log.debug(
                        f"Cannot acquire lock '{target}' - room {room_id} is globally locked by {holder}"
                    )
                    return (
                        jsonify(
                            {
                                "success": False,
                                "error": f"Cannot acquire lock - room is globally locked by {holder}",
                            }
                        ),
                        423,
                    )
            except json.JSONDecodeError:
                log.warning(f"Invalid global lock data in room {room_id}")

    # Use default TTL from config
    ttl = LockConfig.DEFAULT_TTL
    refresh_interval = LockConfig.DEFAULT_REFRESH_INTERVAL

    lock_key = get_lock_key(room_id, target)

    # Generate unique lock token for this acquisition
    lock_token = str(uuid.uuid4())

    # Store lock data (sessionId + userId + token)
    lock_data = json.dumps(
        {"sessionId": session_id, "userId": user_name, "token": lock_token}
    )

    # Acquire lock (nx=True means only set if not exists)
    if r.set(lock_key, lock_data, nx=True, ex=int(ttl)):
        # Track lock in session's lock set for efficient cleanup on disconnect
        session_locks_key = SessionKeys.session_locks(session_id)
        r.sadd(session_locks_key, lock_key)
        # Set TTL on the tracking set to ensure cleanup even if disconnect handler fails
        r.expire(session_locks_key, int(ttl))

        # Store metadata if provided
        timestamp = None
        if msg:
            metadata_key = f"{lock_key}:metadata"
            timestamp = utc_now_iso()
            metadata = {
                "msg": msg,
                "userName": user_name,
                "timestamp": timestamp,
            }
            r.set(metadata_key, json.dumps(metadata), ex=int(ttl))

        # Broadcast lock acquisition event
        emit_lock_update(
            room_id=room_id,
            target=target,
            action="acquired",
            user_name=user_name,
            message=msg,
            timestamp=timestamp,
            session_id=session_id,
        )

        log.debug(
            f"Lock acquired for '{target}' in room '{room_id}' by user {user_name} with token {lock_token[:8]}... and TTL {ttl}s"
        )

        return jsonify(
            {
                "success": True,
                "lockToken": lock_token,
                "ttl": ttl,
                "refreshInterval": refresh_interval,
                "refreshed": False,
            }
        )
    else:
        # Lock acquisition failed - check if we already hold it (idempotent behavior)
        # TODO: TOCTOU race condition - the GET/EXPIRE/SET sequence below is not atomic.
        # Another client could release/re-acquire the lock between operations.
        # Fix requires atomic WATCH/MULTI/EXEC transaction.
        existing_lock_str = r.get(lock_key)
        if existing_lock_str:
            try:
                existing_lock = json.loads(existing_lock_str)
            except json.JSONDecodeError:
                log.warning(
                    f"Corrupted lock data for '{target}' in room '{room_id}': "
                    f"failed to parse JSON: {existing_lock_str!r}"
                )
                lock_holder = "corrupted lock data"
                existing_lock = None

            if existing_lock is not None:
                if not isinstance(existing_lock, dict):
                    log.warning(
                        f"Corrupted lock data for '{target}' in room '{room_id}': "
                        f"expected dict, got {type(existing_lock).__name__}: {existing_lock_str!r}"
                    )
                    lock_holder = "corrupted lock data"
                elif existing_lock.get("sessionId") == session_id:
                    # If same session holds the lock, refresh it (idempotent acquire)
                    # Refresh lock TTL
                    r.expire(lock_key, int(ttl))

                    # Refresh session_locks tracking set TTL (critical for disconnect cleanup)
                    session_locks_key = SessionKeys.session_locks(session_id)
                    r.expire(session_locks_key, int(ttl))

                    # Update metadata if msg provided
                    timestamp = None
                    if msg:
                        metadata_key = f"{lock_key}:metadata"
                        timestamp = utc_now_iso()
                        metadata = {
                            "msg": msg,
                            "userName": user_name,
                            "timestamp": timestamp,
                        }
                        r.set(metadata_key, json.dumps(metadata), ex=int(ttl))
                    else:
                        # Just refresh metadata TTL without updating
                        metadata_key = f"{lock_key}:metadata"
                        if r.exists(metadata_key):
                            r.expire(metadata_key, int(ttl))

                    # Broadcast refresh event
                    emit_lock_update(
                        room_id=room_id,
                        target=target,
                        action="refreshed",
                        user_name=user_name,
                        message=msg,
                        timestamp=timestamp,
                        session_id=session_id,
                    )

                    log.debug(
                        f"Lock refreshed (idempotent acquire) for '{target}' in room '{room_id}' "
                        f"by user {user_name} - same session already held lock"
                    )

                    return jsonify(
                        {
                            "success": True,
                            "lockToken": existing_lock.get("token"),
                            "ttl": ttl,
                            "refreshInterval": refresh_interval,
                            "refreshed": True,
                        }
                    )
                else:
                    # Different session holds the lock
                    lock_holder = existing_lock.get("userId", "unknown")
        else:
            # Lock was released between our SET attempt and this GET
            log.debug(
                f"Lock for '{target}' in room '{room_id}' was released during "
                f"acquire attempt by {user_name} (race condition)"
            )
            lock_holder = "no lock (released)"

        log.debug(
            f"Lock for '{target}' in room '{room_id}' already held by {lock_holder}, denied for {user_name}"
        )
        return (
            jsonify({"success": False, "error": f"Lock already held by {lock_holder}"}),
            423,
        )  # Locked


@locks.route("/api/rooms/<room_id>/locks/<target>/release", methods=["POST"])
def release_lock(room_id, target):
    """Release a lock.

    Parameters
    ----------
    room_id : str
        Room identifier
    target : str
        Lock target (e.g., "trajectory:meta")

    Headers
    -------
    X-Session-ID : str (required)
        Session ID from /join response

    Request Body
    ------------
    lockToken : str
        Lock token from acquire response

    Returns
    -------
    dict
        {"success": true} on success
        {"success": false, "error": "..."} on failure (403 if not lock holder or invalid token)
    """
    try:
        user_name = get_current_user()
    except AuthError as e:
        return jsonify({"success": False, "error": str(e)}), 401

    # Extract session ID from header
    session_id = request.headers.get("X-Session-ID")
    if not session_id:
        return jsonify({"success": False, "error": "Session ID required"}), 400

    r = current_app.extensions["redis"]
    request_data = request.json or {}
    lock_token = request_data.get("lockToken")

    if not lock_token:
        return jsonify({"success": False, "error": "Lock token required"}), 400

    lock_key = get_lock_key(room_id, target)
    lock_data_str = r.get(lock_key)

    if not lock_data_str:
        return jsonify({"success": False, "error": "Lock not held"}), 403

    # Parse lock data
    try:
        lock_data = json.loads(lock_data_str)
    except json.JSONDecodeError:
        return jsonify({"success": False, "error": "Invalid lock data"}), 500

    # Verify sessionId AND token match
    if lock_data.get("sessionId") != session_id or lock_data.get("token") != lock_token:
        log.warning(
            f"Failed release: Lock for '{target}' in room '{room_id}' - session or token mismatch"
        )
        return (
            jsonify(
                {"success": False, "error": "Lock not held by caller or invalid token"}
            ),
            403,
        )

    # Delete lock and metadata
    r.delete(lock_key)
    r.delete(f"{lock_key}:metadata")

    # Remove from session's lock set
    session_locks_key = SessionKeys.session_locks(session_id)
    r.srem(session_locks_key, lock_key)

    # Broadcast lock release event
    emit_lock_update(
        room_id=room_id,
        target=target,
        action="released",
        user_name=None,
        message=None,
        timestamp=None,
        session_id=session_id,
    )

    log.debug(f"Lock released for '{target}' in room '{room_id}' by user {user_name}")

    return jsonify({"success": True})
