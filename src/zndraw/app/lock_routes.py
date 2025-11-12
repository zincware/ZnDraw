"""REST API routes for distributed lock management."""

import datetime
import json
import logging
import uuid

from flask import Blueprint, current_app, jsonify, request

from zndraw.auth import AuthError, get_current_user

from .constants import LockConfig
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

    # Use default TTL from config
    ttl = LockConfig.DEFAULT_TTL
    refresh_interval = LockConfig.DEFAULT_REFRESH_INTERVAL

    lock_key = get_lock_key(room_id, target)

    # Generate unique lock token for this acquisition
    lock_token = str(uuid.uuid4())

    # Store lock data (sessionId + userId + token)
    lock_data = json.dumps({
        "sessionId": session_id,
        "userId": user_name,
        "token": lock_token
    })

    # Acquire lock (nx=True means only set if not exists)
    if r.set(lock_key, lock_data, nx=True, ex=int(ttl)):
        # Store metadata if provided
        timestamp = None
        if msg:
            metadata_key = f"{lock_key}:metadata"
            timestamp = datetime.datetime.utcnow().isoformat()
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
            session_id=session_id
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
            }
        )
    else:
        lock_data_str = r.get(lock_key)
        try:
            lock_data = json.loads(lock_data_str) if lock_data_str else {}
            lock_holder = lock_data.get("userId", "unknown")  # Changed from userName to userId
        except (json.JSONDecodeError, AttributeError):
            lock_holder = "unknown"

        log.info(
            f"Lock for '{target}' in room '{room_id}' already held by {lock_holder}, denied for {user_name}"
        )
        return (
            jsonify(
                {"success": False, "error": f"Lock already held by {lock_holder}"}
            ),
            423,
        )  # Locked


@locks.route("/api/rooms/<room_id>/locks/<target>/refresh", methods=["POST"])
def refresh_lock(room_id, target):
    """Refresh lock TTL and optionally update message.

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
    msg : str, optional
        Optional updated message (if None, keeps existing message)

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
    msg = request_data.get("msg")

    if not lock_token:
        return jsonify({"success": False, "error": "Lock token required"}), 400

    ttl = LockConfig.DEFAULT_TTL

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
            f"Failed refresh: Lock for '{target}' in room '{room_id}' - session or token mismatch"
        )
        return (
            jsonify({"success": False, "error": "Lock not held by caller or invalid token"}),
            403,
        )

    # Refresh lock TTL (keep same lock data with user+token)
    r.set(lock_key, lock_data_str, ex=int(ttl))

    # Update metadata if msg provided
    timestamp = None
    if msg:
        metadata_key = f"{lock_key}:metadata"
        timestamp = datetime.datetime.utcnow().isoformat()
        metadata = {
            "msg": msg,
            "userName": user_name,
            "timestamp": timestamp,
        }
        r.set(metadata_key, json.dumps(metadata), ex=int(ttl))

        log.debug(
            f"Lock refreshed for '{target}' in room '{room_id}' by {user_name} with updated message"
        )
    else:
        # Just refresh metadata TTL without updating - check if metadata exists
        metadata_key = f"{lock_key}:metadata"
        if r.exists(metadata_key):
            r.expire(metadata_key, int(ttl))
            # Get existing metadata for broadcast
            metadata_raw = r.get(metadata_key)
            if metadata_raw:
                existing_metadata = json.loads(metadata_raw)
                msg = existing_metadata.get("msg")
                timestamp = existing_metadata.get("timestamp")

        log.debug(
            f"Lock refreshed for '{target}' in room '{room_id}' by {user_name}"
        )

    # Broadcast lock refresh event
    emit_lock_update(
        room_id=room_id,
        target=target,
        action="refreshed",
        user_name=user_name,
        message=msg,
        timestamp=timestamp,
        session_id=session_id
    )

    return jsonify({"success": True})


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
            jsonify({"success": False, "error": "Lock not held by caller or invalid token"}),
            403,
        )

    # Delete lock and metadata
    r.delete(lock_key)
    r.delete(f"{lock_key}:metadata")

    # Broadcast lock release event
    emit_lock_update(
        room_id=room_id,
        target=target,
        action="released",
        user_name=None,
        message=None,
        timestamp=None,
        session_id=session_id
    )

    log.debug(
        f"Lock released for '{target}' in room '{room_id}' by user {user_name}"
    )

    return jsonify({"success": True})
