import json
import logging
import typing as t
import uuid

from flask import current_app, request
from flask_socketio import emit, join_room, leave_room

from zndraw.server import socketio

from .constants import SocketEvents
from .redis_keys import ExtensionKeys
from .worker_dispatcher import dispatch_next_task

log = logging.getLogger(__name__)

# for crash handling.
TOKEN_EXPIRY_SECONDS = 10


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


def get_client_id_from_sid(sid: str) -> t.Optional[str]:
    """Gets the client_id for a given Socket.IO sid."""
    r = current_app.extensions["redis"]
    client_id = r.get(f"sid:{sid}:client_id")
    return client_id if client_id else None


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

    # Get client_id before cleanup
    client_id = r.get(f"sid:{sid}:client_id")

    # --- New Redis-based Cleanup ---
    room_name = r.get(f"sid:{sid}")
    if room_name:
        # Get username before deleting for logging purposes
        user = r.hget(f"room:{room_name}:users", sid)
        log.info(
            f"Client disconnected: sid={sid}, client_id={client_id}, user={user}, room={room_name}"
        )

        # Clean up Redis entries
        r.hdel(f"room:{room_name}:users", sid)
        r.delete(f"sid:{sid}")

        # Clean up client_id mappings
        if client_id:
            r.delete(f"sid:{sid}:client_id")
            r.delete(f"client_id:{client_id}:sid")

        # (Optional) Notify room that a user has left
        users_in_room = r.hgetall(f"room:{room_name}:users")
        emit("room_users_update", users_in_room, to=f"room:{room_name}")
    else:
        log.info(f"Client disconnected: {sid} (was not in a Redis-managed room)")

    # --- Existing Lock Cleanup Logic ---
    lock_keys = r.scan_iter("*:lock:*")
    for key in lock_keys:
        if r.get(key) == sid:
            log.warning(
                f"Cleaning up orphaned lock '{key}' held by disconnected client {sid}"
            )
            r.delete(key)

    if room_name:
        lock_key = f"room:{room_name}:presenter_lock"
        presenter_sid = r.get(lock_key)
        if presenter_sid and presenter_sid == sid:
            r.delete(lock_key)
            # Inform everyone that the presenter left
            emit("presenter_update", {"presenterSid": None}, to=f"room:{room_name}")

    # Extension cleanup - use client_id if available, otherwise fall back to sid for backwards compatibility
    worker_id = client_id if client_id else sid
    extension_categories = ["modifiers", "selections", "analyses"]
    log.info(
        f"Cleaning up extensions for worker_id={worker_id} in room '{room_name}'..."
    )

    for category in extension_categories:
        user_extensions_key = f"room:{room_name}:extensions:{category}:{worker_id}"
        # This key tells us which extensions this worker_id was providing
        user_extensions = r.smembers(user_extensions_key)

        if not user_extensions:
            continue

        log.info(
            f"Worker {worker_id} was providing extensions in '{category}': {user_extensions}"
        )

        extensions_to_delete = []

        with r.pipeline() as pipe:
            for ext_name in user_extensions:
                keys = ExtensionKeys.for_extension(room_name, category, ext_name)
                idle_key = keys.idle_workers
                progressing_key = keys.progressing_workers

                # Remove the worker_id from both possible state sets
                pipe.srem(idle_key, worker_id)
                pipe.srem(progressing_key, worker_id)

                # Check the combined cardinality of both sets to see if the extension is orphaned
                pipe.scard(idle_key)
                pipe.scard(progressing_key)

            # Each extension now produces 4 results in the pipeline
            results = pipe.execute()

        # Iterate through the results to decide which extensions to delete
        for i, ext_name in enumerate(user_extensions):
            # Get the scard results for this extension
            remaining_idle = results[i * 4 + 2]
            remaining_progressing = results[i * 4 + 3]
            total_remaining = remaining_idle + remaining_progressing

            log.info(
                f"Extension '{ext_name}': {total_remaining} workers remaining after removing {worker_id}."
            )

            # Only delete extension if no workers AND no jobs in queue
            if total_remaining == 0:
                keys = ExtensionKeys.for_extension(room_name, category, ext_name)
                queue_length = r.llen(keys.queue)

                if queue_length == 0:
                    extensions_to_delete.append(ext_name)
                    log.info(
                        f"Extension '{ext_name}' marked for deletion: no workers, no queued jobs"
                    )
                else:
                    log.info(
                        f"Extension '{ext_name}' kept despite no workers: {queue_length} jobs in queue"
                    )

        # If any extensions are now orphaned, delete them and their state sets
        if extensions_to_delete:
            print(
                f"Deleting orphaned extensions in '{category}': {extensions_to_delete}"
            )
            with r.pipeline() as pipe:
                for ext_name in extensions_to_delete:
                    keys = ExtensionKeys.for_extension(room_name, category, ext_name)
                    # Delete the state sets
                    pipe.delete(keys.idle_workers)
                    pipe.delete(keys.progressing_workers)
                    # Delete the schema from the main hash
                    pipe.hdel(keys.schema, ext_name)
                pipe.execute()

        # Clean up the user-specific reverse-lookup key
        r.delete(user_extensions_key)
        print(f"Cleaned up user-specific extension list: {user_extensions_key}")

        # Notify clients about worker count changes
        # We always invalidate if this worker had any extensions, not just when deleting
        if user_extensions:
            print(
                f"Invalidating schema for category '{category}' in room '{room_name}' "
                f"due to worker disconnect."
            )
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"roomId": room_name, "category": category},
                to=f"room:{room_name}",
            )


@socketio.on("join_room")
def on_join(data):
    room = data["room"]
    user = data.get("userId", "Anonymous")
    client_id = data.get("clientId")
    sid = request.sid
    r = current_app.extensions["redis"]

    # Leave previous room if any
    if previous_room_name := r.get(f"sid:{sid}"):
        leave_room(f"room:{previous_room_name}")
        # Remove user from the old room's user list in Redis
        r.hdel(f"room:{previous_room_name}:users", sid)
        log.info(f"Client {sid} ({user}) removed from Redis room: {previous_room_name}")

    room_exists = r.exists(f"room:{room}:template")
    if not room_exists:
        log.error(f"Room '{room}' does not exist. Cannot join non-existent room.")
        return

    # Join the new room (flask-socketio still needs this for the message queue)
    join_room(f"room:{room}")
    # joint the room for this user, e.g. one user can have multiple connections
    join_room(f"user:{user}")

    # Store the new room membership and user info in Redis
    r.set(f"sid:{sid}", room)
    r.hset(f"room:{room}:users", sid, user)

    # Store bidirectional mapping between Socket.IO sid and client_id
    if client_id:
        r.set(f"sid:{sid}:client_id", client_id)
        r.set(f"client_id:{client_id}:sid", sid)
        log.info(
            f"Client {sid} ({user}) with client_id {client_id} joined room: {room}"
        )
    else:
        log.info(f"Client {sid} ({user}) joined room: {room} (no client_id provided)")


@socketio.on("selection:set")
def handle_selection_set(data):
    room = get_project_room_from_session(request.sid)
    redis_client = current_app.extensions["redis"]

    indices = data.get("indices", [])
    if not isinstance(indices, list):
        return {
            "success": False,
            "messagemsg": "Indices must be a list",
            "error": "TypeError",
        }

    # Validate and store the selection in Redis
    # valid_indices = [idx for idx in indices if isinstance(idx, int) and 0 <= idx < _get_len()]
    if any(not isinstance(idx, int) or idx < 0 for idx in indices):
        return {
            "success": False,
            "message": f"No valid indices provided {indices}",
            "error": "ValueError",
        }

    redis_client.set(f"room:{room}:selection:default", json.dumps(indices))
    emit(
        "selection:update",
        {"indices": indices},
        to=f"room:{room}",
        skip_sid=request.sid,
    )
    return {"success": True}


@socketio.on("frame_selection:set")
def handle_frame_selection_set(data):
    room = get_project_room_from_session(request.sid)
    redis_client = current_app.extensions["redis"]

    indices = data.get("indices", [])
    if not isinstance(indices, list):
        return {
            "success": False,
            "messagemsg": "Indices must be a list",
            "error": "TypeError",
        }

    # Validate and store the selection in Redis
    # valid_indices = [idx for idx in indices if isinstance(idx, int) and 0 <= idx < _get_len()]
    if any(not isinstance(idx, int) or idx < 0 for idx in indices):
        return {
            "success": False,
            "message": f"No valid indices provided {indices}",
            "error": "ValueError",
        }

    redis_client.set(f"room:{room}:frame_selection:default", json.dumps(indices))
    emit(
        "frame_selection:update",
        {"indices": indices},
        to=f"room:{room}",
        skip_sid=request.sid,
    )
    return {"success": True}


@socketio.on("bookmarks:set")
def handle_bookmarks_set(data):
    """
    Set bookmarks for frames. Bookmarks are stored using physical indices
    so they persist through frame deletions/reordering.

    data format: {"bookmarks": {logical_index: "label", ...}}
    """
    room = get_project_room_from_session(request.sid)
    redis_client = current_app.extensions["redis"]

    bookmarks = data.get("bookmarks", {})
    if not isinstance(bookmarks, dict):
        return {
            "success": False,
            "message": "Bookmarks must be a dictionary",
            "error": "TypeError",
        }

    # Get frame mapping to convert logical -> physical indices
    indices_key = f"room:{room}:trajectory:indices"
    frame_mapping = redis_client.zrange(indices_key, 0, -1)

    if not frame_mapping:
        max_frame = -1
    else:
        max_frame = len(frame_mapping) - 1

    # Convert logical indices to physical indices and validate
    physical_bookmarks = {}
    for logical_idx_str, label in bookmarks.items():
        try:
            logical_idx = int(logical_idx_str)
        except (ValueError, TypeError):
            return {
                "success": False,
                "message": f"Invalid bookmark index: {logical_idx_str}",
                "error": "ValueError",
            }

        if logical_idx < 0 or logical_idx > max_frame:
            return {
                "success": False,
                "message": f"Bookmark index {logical_idx} out of range [0, {max_frame}]",
                "error": "IndexError",
            }

        if not isinstance(label, str):
            return {
                "success": False,
                "message": f"Bookmark label must be a string, got {type(label).__name__}",
                "error": "TypeError",
            }

        # Get physical index from mapping
        physical_key = frame_mapping[logical_idx]
        physical_bookmarks[physical_key] = label

    # Store bookmarks using physical indices
    bookmarks_key = f"room:{room}:bookmarks"

    # Clear existing bookmarks and set new ones
    redis_client.delete(bookmarks_key)
    if physical_bookmarks:
        redis_client.hset(bookmarks_key, mapping=physical_bookmarks)

    # Emit update with logical indices for clients
    emit(
        "bookmarks:update",
        {"bookmarks": bookmarks},  # Send back the logical indices
        to=f"room:{room}",
        skip_sid=request.sid,
    )
    return {"success": True}


@socketio.on("set_frame_atomic")
def handle_set_frame_atomic(data):
    """
    Handles a single frame jump. REJECTED if a presenter is active.
    """
    room = get_project_room_from_session(request.sid)
    lock_key = f"room:{room}:presenter_lock"
    redis_client = current_app.extensions["redis"]

    if redis_client.get(lock_key) not in [request.sid, None]:
        return {
            "success": False,
            "error": "LockError",
            "message": "Cannot set frame while presenter is active",
        }

    frame = data.get("frame")
    if frame is not None:
        redis_client.set(f"room:{room}:current_frame", frame)
        emit("frame_update", {"frame": frame}, to=f"room:{room}", skip_sid=request.sid)
        return {"success": True}

    return {"success": False, "error": "Invalid frame"}


@socketio.on("set_frame_continuous")
def handle_set_frame_continuous(data):
    """
    Handles continuous frame updates. REQUIRES sender to be the presenter.
    """
    room = get_project_room_from_session(request.sid)
    lock_key = f"room:{room}:presenter_lock"
    redis_client = current_app.extensions["redis"]

    presenter_sid = redis_client.get(lock_key)

    if presenter_sid and presenter_sid == request.sid:
        frame = data.get("frame")
        if frame is not None:
            redis_client.set(f"room:{room}:current_frame", frame)
            emit(
                "frame_update",
                {"frame": frame},
                to=f"room:{room}",
                skip_sid=request.sid,
            )


@socketio.on("request_presenter_token")
def handle_request_presenter_token():
    sid = request.sid
    room = get_project_room_from_session(sid)
    print(f"Presenter token requested by {sid} in room {room}")
    if not room:
        return {"success": False, "reason": "Not in a valid room"}

    lock_key = f"room:{room}:presenter_lock"
    r = current_app.extensions["redis"]

    # --- UPDATED LOGIC ---
    # Get the current holder of the lock
    current_holder = r.get(lock_key)

    # Case 1: No one has the lock, or the requester already has it (renewal)
    if current_holder is None or current_holder == sid:
        # Set (or reset) the lock with the new expiry
        r.set(lock_key, sid, ex=TOKEN_EXPIRY_SECONDS)

        # If this is a brand new presenter, inform the room
        if current_holder is None:
            emit(
                "presenter_update",
                {"presenterSid": sid},
                to=f"room:{room}",
                skip_sid=sid,
            )

        return {"success": True}
    else:
        # Case 2: Someone else has the lock
        return {"success": False, "reason": "Presenter lock is held by another user"}


@socketio.on("release_presenter_token")
def handle_release_presenter_token():
    room = get_project_room_from_session(request.sid)
    lock_key = f"room:{room}:presenter_lock"
    print(f"Presenter token release requested by {request.sid} in room {room}")
    r = current_app.extensions["redis"]

    presenter_sid = r.get(lock_key)

    if presenter_sid and presenter_sid == request.sid:
        r.delete(lock_key)
        emit(
            "presenter_update",
            {"presenterSid": None},
            to=f"room:{room}",
            skip_sid=request.sid,
        )
        return {"success": True}
    else:
        return {"success": False, "error": "Not the current presenter"}


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
                "error_type": "IndexError",
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
                "error_type": "IndexError",
            }

        # Parse the frame mapping entry to check if it's a template frame
        mapping_entry = frame_mapping[frame_id]
        if isinstance(mapping_entry, bytes):
            mapping_entry = mapping_entry.decode()

        # Check if this frame belongs to a different room (e.g., a template)
        if ":" in mapping_entry:
            source_room_id = mapping_entry.split(":", 1)[0]
            if source_room_id != room:
                # This is a template frame or from another room
                return {
                    "success": False,
                    "error": f"Cannot delete template frame. This frame belongs to '{source_room_id}'",
                    "error_type": "PermissionError",
                }

        # Get the physical index that we're "deleting" (just removing from mapping)
        # For backward compatibility, handle both old format (just index) and new format (room:index)
        if ":" in mapping_entry:
            physical_index_to_remove = int(mapping_entry.split(":", 1)[1])
        else:
            physical_index_to_remove = int(mapping_entry)

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

        # Emit bookmarks update to all clients (uses helper from routes)
        from .routes import emit_bookmarks_update, emit_frames_invalidate
        emit_bookmarks_update(room)
        # Invalidate all frames from deleted position onward (they all shift down)
        emit_frames_invalidate(room, operation="delete", affected_from=frame_id)

        return {
            "success": True,
            "deleted_frame": frame_id,
            "physical_preserved": physical_index_to_remove,
        }
    except Exception as e:
        log.error(f"Failed to delete frame: {e}")
        return {"success": False, "error": "Failed to delete frame"}


@socketio.on("register:extension")
def register_extension(data: dict):
    """Registers a new extension for the room and sets its initial state to idle."""
    sid = request.sid
    client_id = get_client_id_from_sid(sid)

    if not client_id:
        log.error(f"Cannot register extension: no client_id found for sid {sid}")
        return {"error": "Client ID not found. Please reconnect."}

    name = data.get("name")
    category = data.get("category")
    schema = data.get("schema")
    public = data.get("public", False)
    if public:
        return {"error": "Cannot register public extensions via this endpoint"}
    room_id = get_project_room_from_session(sid)

    if not name or not category or not schema:
        return {"error": "name, category, and schema are required"}

    # Security check: prevent registration of extensions with server-side names
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
    }

    if category in category_map and name in category_map[category]:
        log.warning(
            f"Blocked attempt to register extension '{name}' in category '{category}' "
            f"- name conflicts with server-side extension (security violation)"
        )
        return {
            "error": f"Cannot register extension '{name}': name is reserved for server-side extensions"
        }

    log.info(
        f"Registering extension for room {room_id}: name={name}, category={category}, "
        f"sid={sid}, client_id={client_id}"
    )

    # store in redis
    redis_client = current_app.extensions["redis"]
    keys = ExtensionKeys.for_extension(room_id, category, name)
    user_extensions_key = ExtensionKeys.user_extensions_key(
        room_id, category, client_id
    )
    existing_schema = redis_client.hget(keys.schema, name)

    if existing_schema is not None:
        existing_schema = json.loads(existing_schema)
        if existing_schema != schema:
            return {
                "error": "Extension with this name already exists with a different schema"
            }
        else:
            # Worker is re-registering. Add it to the idle set.
            redis_client.sadd(keys.idle_workers, client_id)
            redis_client.sadd(
                user_extensions_key, name
            )  # Keep this for disconnect cleanup

            # Notify clients about worker count change
            log.info(
                f"Worker {client_id} re-registered for extension '{name}' "
                f"in category '{category}', invalidating schema"
            )
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"roomId": room_id, "category": category},
                to=f"room:{room_id}",
            )

            # Check if there are queued tasks that this worker can handle
            if room_id:  # Null check for type safety
                dispatch_next_task(redis_client, socketio, client_id, room_id, category)

            return {
                "status": "success",
                "message": "Extension already registered with same schema. Worker marked as idle.",
            }
    else:
        # This is a brand new extension for the room
        with redis_client.pipeline() as pipe:
            # Set the schema
            pipe.hset(keys.schema, name, json.dumps(schema))
            # Add the current worker to the set of idle workers for this extension
            pipe.sadd(keys.idle_workers, client_id)
            # Maintain the reverse mapping for easy cleanup on disconnect
            pipe.sadd(user_extensions_key, name)
            pipe.execute()

        # Invalidate schema on all clients so they can see the new extension
        socketio.emit(
            SocketEvents.INVALIDATE_SCHEMA,
            {"roomId": room_id, "category": category},
            to=f"room:{room_id}",
        )

        # Check if there are queued tasks that this worker can handle
        if room_id:  # Null check for type safety
            dispatch_next_task(redis_client, socketio, client_id, room_id, category)

    return {"status": "success"}
