import logging
import typing as t
import uuid

from flask import current_app, request
from flask_socketio import join_room, leave_room, rooms, emit

from zndraw.server import socketio

log = logging.getLogger(__name__)


TOKEN_EXPIRY_SECONDS = 5000 # Short expiry for auto-cleanup


def _get_len() -> dict:
    sid = request.sid
    r = current_app.extensions["redis"]
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


# --- Helper Functions ---
def get_project_room_from_session(sid: str) -> t.Optional[str]:
    """Finds the project room a client has joined from Redis."""
    r = current_app.extensions["redis"]
    room_name = r.get(f"sid:{sid}")
    return room_name if room_name else None


def get_lock_key(room: str, target: str) -> str:
    """Constructs a standardized Redis key for a lock."""
    return f"room:{room}:lock:{target}"

@socketio.on("connect")
def handle_connect():
    sid = request.sid
    log.info(f"Client connected: {sid}")

@socketio.on("disconnect")
def handle_disconnect():
    sid = request.sid
    r = current_app.extensions["redis"]
    
    # --- New Redis-based Cleanup ---
    room_name = r.get(f"sid:{sid}")
    if room_name:
        # Get username before deleting for logging purposes
        user = r.hget(f"room:{room_name}:users", sid)
        log.info(f"Client disconnected: {sid} ({user}) from room {room_name}")
        
        # Clean up Redis entries
        r.hdel(f"room:{room_name}:users", sid)
        r.delete(f"sid:{sid}")

        # (Optional) Notify room that a user has left
        users_in_room = r.hgetall(f"room:{room_name}:users")
        emit('room_users_update', users_in_room, to=f"room:{room_name}")
    else:
        log.info(f"Client disconnected: {sid} (was not in a Redis-managed room)")

    # --- Existing Lock Cleanup Logic ---
    lock_keys = r.scan_iter("*:lock:*")
    for key in lock_keys:
        if r.get(key) == sid:
            log.warning(f"Cleaning up orphaned lock '{key}' held by disconnected client {sid}")
            r.delete(key)
            
    if room_name:
        lock_key = f"room:{room_name}:presenter_lock"
        presenter_sid = r.get(lock_key)
        if presenter_sid and presenter_sid == sid:
            r.delete(lock_key)
            # Inform everyone that the presenter left
            emit('presenter_update', {'presenterSid': None}, to=f"room:{room_name}")


@socketio.on("join_room")
def on_join(data):
    room = data["room"]
    user = data.get("userId", "Anonymous")
    sid = request.sid
    r = current_app.extensions["redis"]

    # --- New Redis-based Logic ---
    # Leave previous room if any
    if previous_room_name := r.get(f"sid:{sid}"):
        leave_room(f"room:{previous_room_name}")
        # Remove user from the old room's user list in Redis
        r.hdel(f"room:{previous_room_name}:users", sid)
        log.info(f"Client {sid} ({user}) removed from Redis room: {previous_room_name}")

    # Join the new room (flask-socketio still needs this for the message queue)
    join_room(f"room:{room}")
    # joint the room for this user, e.g. one user can have multiple connections
    join_room(f"user:{user}")    

    # Store the new room membership and user info in Redis
    r.set(f"sid:{sid}", room)
    r.hset(f"room:{room}:users", sid, user)
    log.info(f"Client {sid} ({user}) joined room: {room} (and stored in Redis)")

    # --- Existing Logic ---
    emit("len_frames", _get_len(), to=sid)
    presenter_sid = r.get(f"room:{room}:presenter_lock")
    if presenter_sid:
        emit('presenter_update', {'presenterSid': presenter_sid}, to=sid)
    current_frame = r.get(f"room:{room}:current_frame")
    if current_frame is not None:
        emit('frame_update', {'frame': int(current_frame)}, to=sid)
    
    # (Optional but recommended) Notify everyone in the room about the updated user list
    users_in_room = r.hgetall(f"room:{room}:users")
    emit('room_users_update', users_in_room, to=f"room:{room}")


@socketio.on('request_presenter_token')
def handle_request_presenter_token():
    sid = request.sid
    room = get_project_room_from_session(sid)
    lock_key = f"room:{room}:presenter_lock"

    r = current_app.extensions["redis"]
    
    # Try to acquire the lock using SETNX (set if not exists)
    # This is an atomic operation.
    was_set = r.set(lock_key, request.sid, nx=True, ex=TOKEN_EXPIRY_SECONDS)

    if was_set:
        # Success! You are the presenter.
        emit('presenter_token_granted')
        # Inform everyone else in the room who the new presenter is
        emit('presenter_update', {'presenterSid': request.sid}, to=f"room:{room}")
    else:
        # Failure, someone else holds the lock.
        emit('presenter_token_denied')


@socketio.on('set_frame_atomic')
def handle_set_frame_atomic(data):
    """
    Handles a single frame jump. REJECTED if a presenter is active.
    """
    room = get_project_room_from_session(request.sid)
    lock_key = f"room:{room}:presenter_lock"
    redis_client = current_app.extensions["redis"]

    if redis_client.exists(lock_key):
        # Presenter is active, return error
        return {"success": False, "error": "LockError", "message": "Cannot set frame while presenter is active"}

    frame = data.get('frame')
    if frame is not None:
        redis_client.set(f"room:{room}:current_frame", frame)
        emit('frame_update', {'frame': frame}, to=f"room:{room}", skip_sid=request.sid)
        return {"success": True}

@socketio.on('set_frame_continuous')
def handle_set_frame_continuous(data):
    """
    Handles continuous frame updates. REQUIRES sender to be the presenter.
    """
    room = get_project_room_from_session(request.sid)
    lock_key = f"room:{room}:presenter_lock"
    redis_client = current_app.extensions["redis"]
    
    presenter_sid = redis_client.get(lock_key)
    if presenter_sid and presenter_sid == request.sid:
        frame = data.get('frame')
        if frame is not None:
            redis_client.set(f"room:{room}:current_frame", frame)
            emit('frame_update', {'frame': frame}, to=f"room:{room}", skip_sid=request.sid)

@socketio.on('renew_presenter_token')
def handle_renew_presenter_token():
    room = get_project_room_from_session(request.sid)
    lock_key = f"room:{room}:presenter_lock"
    r = current_app.extensions["redis"]
    
    # Verify ownership before renewing
    presenter_sid = r.get(lock_key)
    if presenter_sid and presenter_sid == request.sid:
        # Refresh the expiry time
        r.expire(lock_key, TOKEN_EXPIRY_SECONDS)


@socketio.on('release_presenter_token')
def handle_release_presenter_token():
    room = get_project_room_from_session(request.sid)
    lock_key = f"room:{room}:presenter_lock"
    r = current_app.extensions["redis"]
    
    # Verify ownership before deleting
    presenter_sid = r.get(lock_key)
    if presenter_sid and presenter_sid == request.sid:
        r.delete(lock_key)
        # Inform everyone that scrubbing is over
        emit('presenter_update', {'presenterSid': None}, to=f"room:{room}")

@socketio.on("lock:acquire")
def acquire_lock(data):
    sid = request.sid
    r = current_app.extensions["redis"]
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
    r = current_app.extensions["redis"]
    target = data.get("target")
    room = get_project_room_from_session(sid)

    if not room or not target:
        return {"success": False, "error": "Room or target missing"}

    lock_key = get_lock_key(room, target)
    if r.get(lock_key) == sid:
        r.delete(lock_key)

        # gathering the lock means typically, making updates, so update frame count
        # TODO: later move this to the specific frames changed, because
        # we need to send a frames_changed event anyway
        emit("len_frames", _get_len(), to=f"room:{room}")

        log.info(f"Lock released for '{target}' in room '{room}' by {sid}")
        return {"success": True}

    log.warning(
        f"Failed release: Lock for '{target}' in room '{room}' not held by {sid}"
    )
    return {"success": False}


@socketio.on("upload:prepare")
def handle_upload_prepare(data):
    sid = request.sid
    r = current_app.extensions["redis"]
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



@socketio.on("frame:delete")
def handle_delete_frame(data):
    sid = request.sid
    r = current_app.extensions["redis"]
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
