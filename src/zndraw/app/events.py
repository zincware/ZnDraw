import logging
import typing as t
import uuid

from flask import current_app, request
from flask_socketio import join_room, leave_room, rooms

from zndraw.server import socketio
from zndraw.app import tasks

log = logging.getLogger(__name__)


# --- Helper Functions ---
def get_project_room_from_session(sid: str) -> t.Optional[str]:
    """Finds the project room a client has joined."""
    for room in rooms(sid=sid):
        if room != sid:
            return room
    return None


def get_lock_key(room: str, target: str) -> str:
    """Constructs a standardized Redis key for a lock."""
    return f"room:{room}:lock:{target}"


@socketio.on("disconnect")
def handle_disconnect():
    sid = request.sid
    r = current_app.config["redis"]
    log.info(f"Client disconnected: {sid}")
    lock_keys = r.scan_iter(f"*:lock:*")
    for key in lock_keys:
        if r.get(key) == sid:
            log.warning(
                f"Cleaning up orphaned lock '{key}' held by disconnected client {sid}"
            )
            r.delete(key)


@socketio.on("join_room")
def handle_join(data):
    room = data["room"]
    sid = request.sid
    if previous_room := get_project_room_from_session(sid):
        leave_room(previous_room)
        log.info(f"Client {sid} left room: {previous_room}")

    join_room(room)
    log.info(f"Client {sid} joined room: {room}")
    tasks.get_schema.delay(sid)


@socketio.on("lock:acquire")
def acquire_lock(data):
    sid = request.sid
    r = current_app.config["redis"]
    target = data.get("target")
    room = get_project_room_from_session(sid)

    if not room or not target:
        return {"success": False, "error": "Room or target missing"}

    lock_key = get_lock_key(room, target)
    if r.set(lock_key, sid, nx=True, ex=60):
        log.info(f"Lock acquired for '{target}' in room '{room}' by {sid}")
        return {"success": True}
    else:
        log.info(
            f"Lock for '{target}' in room '{room}' already held by {r.get(lock_key)}, denied for {sid}"
        )
        return {"success": False}


@socketio.on("lock:release")
def release_lock(data):
    sid = request.sid
    r = current_app.config["redis"]
    target = data.get("target")
    room = get_project_room_from_session(sid)

    if not room or not target:
        return {"success": False, "error": "Room or target missing"}

    lock_key = get_lock_key(room, target)
    if r.get(lock_key) == sid:
        r.delete(lock_key)
        log.info(f"Lock released for '{target}' in room '{room}' by {sid}")
        return {"success": True}

    log.warning(
        f"Failed release: Lock for '{target}' in room '{room}' not held by {sid}"
    )
    return {"success": False}


@socketio.on("upload:prepare")
def handle_upload_prepare(data):
    sid = request.sid
    r = current_app.config["redis"]
    room = get_project_room_from_session(sid)
    action = data.get("action")  # Default to append for backward compatibility
    frame_id = data.get("frame_id")  # For replace operations
    insert_position = data.get("insert_position")  # For insert operations

    if action not in {"append", "replace", "insert", "extend"}:
        return {"success": False, "error": "Invalid action specified"}

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    lock_key = get_lock_key(room, "trajectory:meta")
    if r.get(lock_key) != sid:
        return {
            "success": False,
            "error": "Client does not hold the trajectory lock.",
        }

    # For replace operations, validate the frame exists
    if action == "replace":
        if frame_id is None:
            return {
                "success": False,
                "error": "frame_id is required for replace operations",
            }

        indices_key = f"room:{room}:trajectory:indices"
        # Use zcard for an efficient count of logical frames
        num_frames = r.zcard(indices_key)

        if frame_id >= num_frames:
            return {
                "success": False,
                "error": f"Frame {frame_id} does not exist. Max index is {num_frames - 1}.",
            }

    token = str(uuid.uuid4())
    token_key = f"room:{room}:upload_token:{token}"
    # Store additional metadata with the token
    token_data = {
        "sid": sid,
        "action": action,
        "frame_id": str(frame_id) if frame_id is not None else "-1",
    }

    # Add insert_position for insert operations
    if insert_position is not None:
        token_data["insert_position"] = str(insert_position)
    r.hset(token_key, mapping=token_data)
    r.expire(token_key, 60)

    log.info(
        f"Issued {action} token for room '{room}' to {sid}"
        + (f" (frame {frame_id})" if frame_id is not None else "")
    )
    return {"success": True, "token": token}


@socketio.on("frames:count")
def handle_len_frames(data):
    sid = request.sid
    r = current_app.config["redis"]
    room = get_project_room_from_session(sid)

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    try:
        indices_key = f"room:{room}:trajectory:indices"
        # Count is the number of entries in the mapping (logical frames)
        frame_count = r.zcard(indices_key)
        return {"success": True, "count": frame_count}
    except Exception as e:
        log.error(f"Failed to get frame count: {e}")
        return {"success": False, "error": "Failed to get frame count"}


@socketio.on("frame:delete")
def handle_delete_frame(data):
    sid = request.sid
    r = current_app.config["redis"]
    room = get_project_room_from_session(sid)
    frame_id = data.get("frame_id")

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    if frame_id is None:
        return {"success": False, "error": "frame_id is required"}

    lock_key = get_lock_key(room, "trajectory:meta")
    if r.get(lock_key) != sid:
        return {
            "success": False,
            "error": "Client does not hold the trajectory lock.",
        }

    try:
        indices_key = f"room:{room}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {"success": False, "error": "No frames found in room"}

        if frame_id >= len(frame_mapping):
            return {
                "success": False,
                "error": f"Frame {frame_id} not found, max frame: {len(frame_mapping) - 1}",
            }

        # Get the physical index that we're "deleting" (just removing from mapping)
        physical_index_to_remove = int(frame_mapping[frame_id])

        # Remove the mapping entry for this logical position
        # We need to rebuild the mapping without this entry
        remaining_physical_indices = (
            frame_mapping[:frame_id] + frame_mapping[frame_id + 1 :]
        )

        # Clear and rebuild the Redis mapping
        r.delete(indices_key)
        for logical_pos, physical_idx_str in enumerate(remaining_physical_indices):
            r.zadd(indices_key, {physical_idx_str: logical_pos})

        log.info(
            f"Deleted logical frame {frame_id} (physical: {physical_index_to_remove}) from room '{room}'. Physical data preserved."
        )
        return {
            "success": True,
            "deleted_frame": frame_id,
            "physical_preserved": physical_index_to_remove,
        }
    except Exception as e:
        log.error(f"Failed to delete frame: {e}")
        return {"success": False, "error": "Failed to delete frame"}
