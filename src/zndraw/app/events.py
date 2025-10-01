import logging
import typing as t
import uuid

from flask import current_app, request
from flask_socketio import join_room, leave_room, rooms, emit
import json

from zndraw.server import socketio

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

    extension_categories = ["modifiers", "selections", "analyses"]
    print(f"Cleaning up extensions for disconnected SID {sid} in room '{room_name}'...")

    for category in extension_categories:
        user_extensions_key = f"room:{room_name}:extensions:{category}:{sid}"
        # This key tells us which extensions this SID was providing
        user_extensions = r.smembers(user_extensions_key)

        if not user_extensions:
            continue

        print(f"User {sid} was providing extensions in '{category}': {user_extensions}")

        schema_invalidated = False
        extensions_to_delete = []

        with r.pipeline() as pipe:
            for ext_name in user_extensions:
                idle_key = (
                    f"room:{room_name}:extensions:{category}:{ext_name}:idle_workers"
                )
                progressing_key = f"room:{room_name}:extensions:{category}:{ext_name}:progressing_workers"

                # Remove the SID from both possible state sets
                pipe.srem(idle_key, sid)
                pipe.srem(progressing_key, sid)

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

            print(
                f"Extension '{ext_name}': {total_remaining} workers remaining after removing {sid}."
            )

            if total_remaining == 0:
                extensions_to_delete.append(ext_name)
                schema_invalidated = True

        # If any extensions are now orphaned, delete them and their state sets
        if extensions_to_delete:
            print(
                f"Deleting orphaned extensions in '{category}': {extensions_to_delete}"
            )
            with r.pipeline() as pipe:
                for ext_name in extensions_to_delete:
                    # Delete the state sets
                    pipe.delete(
                        f"room:{room_name}:extensions:{category}:{ext_name}:idle_workers"
                    )
                    pipe.delete(
                        f"room:{room_name}:extensions:{category}:{ext_name}:progressing_workers"
                    )
                    # Delete the schema from the main hash
                    pipe.hdel(f"room:{room_name}:extensions:{category}", ext_name)
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
                "invalidate:schema",
                {"roomId": room_name, "category": category},
                to=f"room:{room_name}",
            )


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
    # presenter_sid = r.get(f"room:{room}:presenter_lock")
    # if presenter_sid:
    #     emit('presenter_update', {'presenterSid': presenter_sid}, to=sid)
    # current_frame = r.get(f"room:{room}:current_frame")
    # if current_frame is not None:
    #     emit('frame_update', {'frame': int(current_frame)}, to=sid)

    # # (Optional but recommended) Notify everyone in the room about the updated user list
    # users_in_room = r.hgetall(f"room:{room}:users")
    # emit('room_users_update', users_in_room, to=f"room:{room}")

    # check if this user has any extensions registered


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


@socketio.on("register:extension")
def register_extension(data: dict):
    """Registers a new extension for the room and sets its initial state to idle."""
    name = data.get("name")
    category = data.get("category")
    schema = data.get("schema")
    public = data.get("public", False)
    if public:
        return {"error": "Cannot register public extensions via this endpoint"}
    room_id = get_project_room_from_session(request.sid)

    if not name or not category or not schema:
        return {"error": "name, category, and schema are required"}

    print(
        f"Registering extension for room {room_id}: name={name}, category={category}, sid={request.sid}"
    )

    # store in redis
    redis_client = current_app.extensions["redis"]
    schema_key = f"room:{room_id}:extensions:{category}"
    existing_schema = redis_client.hget(schema_key, name)

    # Define keys for the new state model
    idle_workers_key = f"room:{room_id}:extensions:{category}:{name}:idle_workers"
    user_extensions_key = f"room:{room_id}:extensions:{category}:{request.sid}"

    if existing_schema is not None:
        existing_schema = json.loads(existing_schema)
        if existing_schema != schema:
            return {
                "error": "Extension with this name already exists with a different schema"
            }
        else:
            # Worker is re-registering. Add it to the idle set.
            redis_client.sadd(idle_workers_key, request.sid)
            redis_client.sadd(
                user_extensions_key, name
            )  # Keep this for disconnect cleanup

            # Notify clients about worker count change
            log.info(
                f"Worker {request.sid} re-registered for extension '{name}' "
                f"in category '{category}', invalidating schema"
            )
            socketio.emit(
                "invalidate:schema",
                {"roomId": room_id, "category": category},
                to=f"room:{room_id}",
            )

            # Check if there are queued tasks that this worker can handle
            _process_worker_queues(request.sid, room_id, category, redis_client)

            return {
                "status": "success",
                "message": "Extension already registered with same schema. Worker marked as idle.",
            }
    else:
        # This is a brand new extension for the room
        with redis_client.pipeline() as pipe:
            # Set the schema
            pipe.hset(schema_key, name, json.dumps(schema))
            # Add the current worker to the set of idle workers for this extension
            pipe.sadd(idle_workers_key, request.sid)
            # Maintain the reverse mapping for easy cleanup on disconnect
            pipe.sadd(user_extensions_key, name)
            pipe.execute()

        # Invalidate schema on all clients so they can see the new extension
        socketio.emit(
            "invalidate:schema",
            {"roomId": room_id, "category": category},
            to=f"room:{room_id}",
        )

        # Check if there are queued tasks that this worker can handle
        _process_worker_queues(request.sid, room_id, category, redis_client)

    return {"status": "success"}


def _process_worker_queues(sid: str, room_id: str, category: str, redis_client) -> None:
    """Process queued tasks for all extensions that a worker can handle.

    This function:
    1. Retrieves all extensions the worker is registered for
    2. For each extension, checks if there's a queued task
    3. Dispatches the first available queued task to the worker
    4. Uses FIFO ordering within each extension's queue

    Args:
        sid: The worker's session ID
        room_id: The room ID
        category: The extension category (e.g., 'modifiers', 'selections')
        redis_client: The Redis client instance
    """
    # Get all extensions this worker is registered for in this category
    # TODO: consider all categories?
    user_extensions_key = f"room:{room_id}:extensions:{category}:{sid}"
    registered_extensions = redis_client.smembers(user_extensions_key)

    if not registered_extensions:
        print(f"Worker {sid} has no registered extensions in category '{category}'")
        return

    print(f"Worker {sid} can handle extensions: {registered_extensions}")

    # Check each extension's queue for pending tasks
    for extension_name in registered_extensions:
        queue_key = f"room:{room_id}:extensions:{category}:{extension_name}:queue"
        idle_key = f"room:{room_id}:extensions:{category}:{extension_name}:idle_workers"
        progressing_key = f"room:{room_id}:extensions:{category}:{extension_name}:progressing_workers"

        # Check if worker is still idle for this extension (might have been assigned already)
        is_idle = redis_client.sismember(idle_key, sid)
        if not is_idle:
            continue

        # Try to pop a task from the queue (LPOP for FIFO)
        queued_task_json = redis_client.lpop(queue_key)

        if queued_task_json:
            # Parse the queued task
            print(f"Found queued task for extension '{extension_name}' for worker {sid}")
            try:
                queued_task = json.loads(queued_task_json)
                task_data = queued_task.get("data")
                # user_id can be used for future analytics/logging
                _user_id = queued_task.get("user_id")

                log.info(
                    f"Dispatching queued task for extension '{extension_name}' "
                    f"to worker {sid} in room '{room_id}'"
                )

                # Move worker from idle to progressing atomically
                moved = redis_client.smove(idle_key, progressing_key, sid)

                if moved:
                    # Emit the task to the worker
                    emit(
                        "task:run",
                        {
                            "data": task_data,
                            "extension": extension_name,
                            "category": category,
                        },
                        to=sid,
                    )

                    print(
                        f"Successfully dispatched task from queue to worker {sid} "
                        f"for extension '{extension_name}'"
                    )

                    # Notify all clients in room about queue update and worker state change
                    new_queue_length = redis_client.llen(queue_key)
                    idle_count = redis_client.scard(idle_key)
                    progressing_count = redis_client.scard(progressing_key)
                    emit(
                        "queue:update",
                        {
                            "roomId": room_id,
                            "category": category,
                            "extension": extension_name,
                            "queueLength": new_queue_length,
                            "idleWorkers": idle_count,
                            "progressingWorkers": progressing_count
                        },
                        to=f"room:{room_id}",
                    )

                    # The worker is now busy, so we stop checking other queues
                    return
                else:
                    # Worker was removed from idle set by another process
                    # Put the task back in the queue
                    redis_client.lpush(queue_key, queued_task_json)
                    log.warning(
                        f"Worker {sid} no longer idle, task returned to queue "
                        f"for '{extension_name}'"
                    )
                    return

            except json.JSONDecodeError as e:
                log.error(f"Failed to parse queued task: {e}")
                continue

    print(f"No queued tasks found for worker {sid} in category '{category}'")


@socketio.on("task:finished")
def handle_task_finished(data: dict):
    """Handles the completion of a task and sets the worker state back to idle.

    After marking the worker as idle, checks for queued tasks across ALL extensions
    this worker can handle and dispatches the next task if available.
    """
    sid = request.sid
    room_id = get_project_room_from_session(sid)
    category = data.get("category")
    extension = data.get("extension")

    if not all([room_id, category, extension]):
        log.warning(f"Received incomplete task_finished event from {sid}: {data}")
        return

    redis_client = current_app.extensions["redis"]
    idle_key = f"room:{room_id}:extensions:{category}:{extension}:idle_workers"
    progressing_key = f"room:{room_id}:extensions:{category}:{extension}:progressing_workers"

    # Atomically move the worker from the progressing set back to the idle set
    moved = redis_client.smove(progressing_key, idle_key, sid)
    if moved:
        log.info(f"Worker {sid} finished task '{extension}' and is now idle.")

        # Notify all clients about worker state change
        idle_count = redis_client.scard(idle_key)
        progressing_count = redis_client.scard(progressing_key)
        queue_key = f"room:{room_id}:extensions:{category}:{extension}:queue"
        queue_length = redis_client.llen(queue_key)

        emit(
            "queue:update",
            {
                "roomId": room_id,
                "category": category,
                "extension": extension,
                "queueLength": queue_length,
                "idleWorkers": idle_count,
                "progressingWorkers": progressing_count
            },
            to=f"room:{room_id}",
        )
    else:
        # This could happen if the worker disconnected and was already cleaned up.
        log.warning(
            f"Worker {sid} finished task '{extension}', but was not in the progressing set."
        )
        return

    # Process queues for all extensions this worker can handle
    _process_worker_queues(sid, room_id, category, redis_client)