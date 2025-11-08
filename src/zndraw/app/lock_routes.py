"""REST API routes for distributed lock management."""

import datetime
import json
import logging

from flask import Blueprint, current_app, jsonify, request

from zndraw.auth import AuthError, get_current_user

from .constants import LockConfig
from .room_manager import emit_room_update
from .route_utils import get_lock_key

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

    Request Body
    ------------
    msg : str, optional
        Optional message describing the lock purpose

    Returns
    -------
    dict
        {"success": true, "ttl": 60, "refreshInterval": 30} on success
        {"success": false, "error": "..."} on failure (423 if locked)
    """
    try:
        user_name = get_current_user()
    except AuthError as e:
        return jsonify({"success": False, "error": str(e)}), 401

    r = current_app.extensions["redis"]
    msg = request.json.get("msg") if request.json else None

    # Use default TTL from config
    ttl = LockConfig.DEFAULT_TTL
    refresh_interval = LockConfig.DEFAULT_REFRESH_INTERVAL

    lock_key = get_lock_key(room_id, target)

    # Acquire lock (nx=True means only set if not exists)
    if r.set(lock_key, user_name, nx=True, ex=int(ttl)):
        # Store metadata if provided
        if msg:
            metadata_key = f"{lock_key}:metadata"
            metadata = {
                "msg": msg,
                "userName": user_name,
                "timestamp": datetime.datetime.utcnow().isoformat(),
            }
            r.set(metadata_key, json.dumps(metadata), ex=int(ttl))

            # Broadcast lock acquisition for trajectory:meta locks
            if target == "trajectory:meta":
                from . import events

                emit_room_update(
                    events.socketio, room_id, metadataLocked=metadata, skip_sid=None
                )

        log.debug(
            f"Lock acquired for '{target}' in room '{room_id}' by user {user_name} with TTL {ttl}s"
        )

        return jsonify(
            {"success": True, "ttl": ttl, "refreshInterval": refresh_interval}
        )
    else:
        lock_holder = r.get(lock_key)
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

    Request Body
    ------------
    msg : str, optional
        Optional updated message (if None, keeps existing message)

    Returns
    -------
    dict
        {"success": true} on success
        {"success": false, "error": "..."} on failure (403 if not lock holder)
    """
    try:
        user_name = get_current_user()
    except AuthError as e:
        return jsonify({"success": False, "error": str(e)}), 401

    r = current_app.extensions["redis"]
    msg = request.json.get("msg") if request.json else None
    ttl = LockConfig.DEFAULT_TTL

    lock_key = get_lock_key(room_id, target)
    lock_holder = r.get(lock_key)

    # Verify caller holds the lock
    if lock_holder != user_name:
        log.warning(
            f"Failed refresh: Lock for '{target}' in room '{room_id}' held by {lock_holder}, not by {user_name}"
        )
        return (
            jsonify({"success": False, "error": "Lock not held by caller"}),
            403,
        )

    # Refresh lock TTL
    r.expire(lock_key, int(ttl))

    # Update metadata if msg provided
    if msg:
        metadata_key = f"{lock_key}:metadata"
        metadata = {
            "msg": msg,
            "userName": user_name,
            "timestamp": datetime.datetime.utcnow().isoformat(),
        }
        r.set(metadata_key, json.dumps(metadata), ex=int(ttl))

        # Broadcast update for trajectory:meta locks
        if target == "trajectory:meta":
            from . import events

            emit_room_update(
                events.socketio, room_id, metadataLocked=metadata, skip_sid=None
            )

        log.debug(
            f"Lock refreshed for '{target}' in room '{room_id}' by {user_name} with updated message"
        )
    else:
        # Just refresh metadata TTL without updating
        metadata_key = f"{lock_key}:metadata"
        if r.exists(metadata_key):
            r.expire(metadata_key, int(ttl))

        log.debug(
            f"Lock refreshed for '{target}' in room '{room_id}' by {user_name}"
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

    Returns
    -------
    dict
        {"success": true} on success
        {"success": false, "error": "..."} on failure (403 if not lock holder)
    """
    try:
        user_name = get_current_user()
    except AuthError as e:
        return jsonify({"success": False, "error": str(e)}), 401

    r = current_app.extensions["redis"]
    lock_key = get_lock_key(room_id, target)
    lock_holder = r.get(lock_key)

    # Verify caller holds the lock
    if lock_holder != user_name:
        log.warning(
            f"Failed release: Lock for '{target}' in room '{room_id}' held by {lock_holder}, not by {user_name}"
        )
        return (
            jsonify({"success": False, "error": "Lock not held by caller"}),
            403,
        )

    # Delete lock and metadata
    r.delete(lock_key)
    r.delete(f"{lock_key}:metadata")

    # Broadcast lock release for trajectory:meta locks
    if target == "trajectory:meta":
        from . import events

        emit_room_update(events.socketio, room_id, metadataLocked=None, skip_sid=None)

    log.debug(
        f"Lock released for '{target}' in room '{room_id}' by user {user_name}"
    )

    return jsonify({"success": True})
