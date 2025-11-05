import datetime
import json
import logging
import time
import typing as t
import uuid

from flask import current_app, request
from flask_socketio import emit

from zndraw.server import socketio

from .constants import SocketEvents
from .models import LockMetadata
from .redis_keys import ExtensionKeys, FilesystemKeys, RoomKeys, SessionKeys, UserKeys
from .room_manager import emit_room_update

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
        keys = RoomKeys(room)
        # Count is the number of entries in the mapping (logical frames)
        frame_count = r.zcard(keys.trajectory_indices())
        return {"success": True, "count": frame_count}
    except Exception as e:
        log.error(f"Failed to get frame count: {e}")
        return {"success": False, "error": "Failed to get frame count"}


# --- Helper Functions ---
def get_project_room_from_session(sid: str) -> t.Optional[str]:
    """Finds the project room a user has joined from Redis."""
    r = current_app.extensions["redis"]
    # Get userName from sid
    session_keys = SessionKeys(sid)
    user_name = r.get(session_keys.username())
    if not user_name:
        return None
    # Get room from user data
    user_keys = UserKeys(user_name)
    room_name = r.hget(user_keys.hash_key(), "currentRoom")
    return room_name if room_name else None


def get_user_name_from_sid(sid: str) -> t.Optional[str]:
    """Gets the userName for a given Socket.IO sid."""
    r = current_app.extensions["redis"]
    session_keys = SessionKeys(sid)
    user_name = r.get(session_keys.username())
    return user_name if user_name else None


def get_lock_key(room: str, target: str) -> str:
    """Constructs a standardized Redis key for a lock."""
    keys = RoomKeys(room)
    return keys.lock(target)


@socketio.on("connect")
def handle_connect(auth):
    """Handle socket connection with JWT authentication.

    Auth payload
    ------------
    {
        "token": "jwt-token-string"
    }
    """
    from flask_socketio import ConnectionRefusedError, join_room

    from zndraw.auth import AuthError, decode_jwt_token

    sid = request.sid
    r = current_app.extensions["redis"]

    # Get JWT token from auth
    token = auth.get("token") if auth else None

    if not token:
        log.warning(f"Client {sid} connected without JWT token")
        raise ConnectionRefusedError("Authentication token required")

    # Validate JWT token
    try:
        payload = decode_jwt_token(token)
        user_name = payload["sub"]  # userName is the primary identifier
        user_role = payload.get("role", "guest")  # Get role from JWT
    except AuthError as e:
        log.error(f"Client {sid} authentication failed: {e.message}")
        raise ConnectionRefusedError(e.message)

    # Verify user exists in Redis
    user_keys = UserKeys(user_name)
    if not r.exists(user_keys.hash_key()):
        log.error(f"User {user_name} not found in Redis")
        raise ConnectionRefusedError("User not registered. Call /api/login first.")

    # Update user's current SID
    r.hset(
        user_keys.hash_key(),
        mapping={
            "currentSid": sid,
            "lastActivity": datetime.datetime.utcnow().isoformat(),
        },
    )

    # Register SID → userName mapping and role mapping
    session_keys = SessionKeys(sid)
    r.set(session_keys.username(), user_name)
    r.set(session_keys.role(), user_role)  # Store role for this session

    # Get user's current room (if any)
    current_room = r.hget(user_keys.hash_key(), "currentRoom")

    if current_room:
        # Rejoin room after reconnection
        join_room(f"room:{current_room}")
        join_room(f"user:{user_name}")
        log.info(f"User {user_name} reconnected to room {current_room} (sid: {sid})")
    else:
        log.info(f"User {user_name} connected but not in any room yet (sid: {sid})")

    return {"status": "ok", "userName": user_name}


@socketio.on("disconnect")
def handle_disconnect():
    """Handle client disconnect.

    Note: No parameters needed - Flask-SocketIO provides request.sid automatically.
    The framework may pass arguments but we don't use them.
    """
    sid = request.sid
    r = current_app.extensions["redis"]

    # Get userName from connection lookup
    session_keys = SessionKeys(sid)
    user_name = r.get(session_keys.username())

    if not user_name:
        log.info(f"Client disconnected: {sid} (no userName found)")
        return

    # Get user data
    user_keys = UserKeys(user_name)
    room_name = r.hget(user_keys.hash_key(), "currentRoom")

    log.info(f"User disconnected: sid={sid}, user={user_name}, room={room_name}")

    # Clean up connection lookup
    r.delete(session_keys.username())

    # Update user's currentSid to empty (user still exists but disconnected)
    r.hset(user_keys.hash_key(), "currentSid", "")

    # Note: We don't remove the user from the room or delete user data
    # The user may reconnect and rejoin the same room
    # Only when a user explicitly joins a different room do we remove them from the old room

    if room_name:
        # Notify room that a user has disconnected (but not left)
        client_service = current_app.extensions["client_service"]
        users_in_room = client_service.get_room_users(room_name)
        emit(
            "room_clients_update",
            {"clients": list(users_in_room)},
            to=f"room:{room_name}",
        )
    else:
        log.info(f"User {user_name} disconnected (was not in a room)")

    # --- Existing Lock Cleanup Logic ---
    lock_keys = r.scan_iter("*:lock:*")
    for key in lock_keys:
        if r.get(key) == sid:
            log.warning(
                f"Cleaning up orphaned lock '{key}' held by disconnected user {sid}"
            )
            r.delete(key)

    if room_name:
        room_keys = RoomKeys(room_name)
        presenter_sid = r.get(room_keys.presenter_lock())
        if presenter_sid and presenter_sid == sid:
            r.delete(room_keys.presenter_lock())
            # Inform everyone that the presenter left via room:update
            from zndraw.app.room_manager import emit_room_update

            emit_room_update(socketio, room_name, skip_sid=sid, presenterSid=None)

    # Extension cleanup - workers are tracked by server's sid
    if not sid:
        log.error(f"No sid during disconnect cleanup")
        return

    worker_id = sid  # Workers are tracked by server's socket sid
    extension_categories = ["modifiers", "selections", "analysis"]
    log.info(
        f"Cleaning up extensions for worker_id={worker_id} (user={user_name}) in room '{room_name}'..."
    )

    for category in extension_categories:
        worker_extensions_key = ExtensionKeys.user_extensions_key(
            room_name, category, worker_id
        )
        # This key tells us which extensions this worker_id (sid) was providing
        worker_extensions = list(r.smembers(worker_extensions_key))

        if not worker_extensions:
            continue

        log.info(
            f"Worker {worker_id} was providing extensions in '{category}': {worker_extensions}"
        )

        extensions_to_delete = []

        with r.pipeline() as pipe:
            for ext_name in worker_extensions:
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
        for i, ext_name in enumerate(worker_extensions):
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
            log.info(
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

        # Clean up the worker-specific reverse-lookup key
        r.delete(worker_extensions_key)
        print(f"Cleaned up worker-specific extension list: {worker_extensions_key}")

        # Notify clients about worker count changes
        # We always invalidate if this worker had any extensions, not just when deleting
        if worker_extensions:
            print(
                f"Invalidating schema for category '{category}' in room '{room_name}' "
                f"due to worker disconnect."
            )
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
                to=f"room:{room_name}",
            )

    # Clean up global (public) extensions
    log.info(f"Checking for global extensions from worker_id={worker_id}...")
    for category in extension_categories:
        # Check if this worker registered any global extensions
        global_worker_extensions_key = ExtensionKeys.global_user_extensions_key(
            category, worker_id
        )
        global_worker_extensions = list(r.smembers(global_worker_extensions_key))

        if not global_worker_extensions:
            continue

        log.info(
            f"Worker {worker_id} was providing global extensions in '{category}': {global_worker_extensions}"
        )

        global_extensions_to_delete = []

        with r.pipeline() as pipe:
            for ext_name in global_worker_extensions:
                keys = ExtensionKeys.for_global_extension(category, ext_name)
                idle_key = keys.idle_workers
                progressing_key = keys.progressing_workers

                # Remove the worker_id from both possible state sets
                pipe.srem(idle_key, worker_id)
                pipe.srem(progressing_key, worker_id)

                # Check the combined cardinality to see if the extension is orphaned
                pipe.scard(idle_key)
                pipe.scard(progressing_key)

            # Each extension produces 4 results
            results = pipe.execute()

        # Iterate through results to decide which extensions to delete
        for i, ext_name in enumerate(global_worker_extensions):
            remaining_idle = results[i * 4 + 2]
            remaining_progressing = results[i * 4 + 3]
            total_remaining = remaining_idle + remaining_progressing

            log.info(
                f"Global extension '{ext_name}': {total_remaining} workers remaining after removing {worker_id}."
            )

            # Only delete if no workers AND no jobs in queue
            if total_remaining == 0:
                keys = ExtensionKeys.for_global_extension(category, ext_name)
                queue_length = r.llen(keys.queue)

                if queue_length == 0:
                    global_extensions_to_delete.append(ext_name)
                    log.info(
                        f"Global extension '{ext_name}' marked for deletion: no workers, no queued jobs"
                    )
                else:
                    log.info(
                        f"Global extension '{ext_name}' kept despite no workers: {queue_length} jobs in queue"
                    )

        # Delete orphaned global extensions
        if global_extensions_to_delete:
            log.info(
                f"Deleting orphaned global extensions in '{category}': {global_extensions_to_delete}"
            )
            with r.pipeline() as pipe:
                for ext_name in global_extensions_to_delete:
                    keys = ExtensionKeys.for_global_extension(category, ext_name)
                    # Delete the state sets
                    pipe.delete(keys.idle_workers)
                    pipe.delete(keys.progressing_workers)
                    # Delete the schema from the main hash
                    pipe.hdel(keys.schema, ext_name)
                pipe.execute()

        # Clean up the worker-specific reverse-lookup key for global extensions
        r.delete(global_worker_extensions_key)
        log.info(
            f"Cleaned up global worker-specific extension list: {global_worker_extensions_key}"
        )

        # Notify ALL clients about global extension removal (broadcast globally)
        if global_worker_extensions:
            log.info(
                f"Broadcasting global schema invalidation for category '{category}' due to worker disconnect"
            )
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category, "global": True},
            )

    # Filesystem cleanup - similar to extensions
    if room_name:
        worker_filesystems_key = FilesystemKeys.user_filesystems_key(room_name, worker_id)
        worker_filesystems = list(r.smembers(worker_filesystems_key))

        if worker_filesystems:
            log.info(
                f"Worker {worker_id} was providing filesystems: {worker_filesystems}"
            )
            for fs_name in worker_filesystems:
                keys = FilesystemKeys.for_filesystem(room_name, fs_name)
                # Delete filesystem metadata and worker reference
                r.delete(keys.metadata)
                r.delete(keys.worker)
                log.info(f"Deleted room-scoped filesystem '{fs_name}' for worker {worker_id}")

            # Clean up worker-specific filesystem list
            r.delete(worker_filesystems_key)
            log.info(f"Cleaned up worker-specific filesystem list: {worker_filesystems_key}")

            # Notify clients in the room that filesystems changed
            socketio.emit(
                SocketEvents.FILESYSTEMS_UPDATE,
                {"scope": "room"},
                to=f"room:{room_name}",
            )

    # Clean up global (public) filesystems
    global_worker_filesystems_key = FilesystemKeys.global_user_filesystems_key(worker_id)
    global_worker_filesystems = list(r.smembers(global_worker_filesystems_key))

    if global_worker_filesystems:
        log.info(
            f"Worker {worker_id} was providing global filesystems: {global_worker_filesystems}"
        )
        for fs_name in global_worker_filesystems:
            keys = FilesystemKeys.for_global_filesystem(fs_name)
            # Delete filesystem metadata and worker reference
            r.delete(keys.metadata)
            r.delete(keys.worker)
            log.info(f"Deleted global filesystem '{fs_name}' for worker {worker_id}")

        # Clean up worker-specific filesystem list
        r.delete(global_worker_filesystems_key)
        log.info(f"Cleaned up global worker-specific filesystem list: {global_worker_filesystems_key}")

        # Notify all clients that global filesystems changed
        socketio.emit(SocketEvents.FILESYSTEMS_UPDATE, {"scope": "global"})


@socketio.on("extension:register")
def handle_extension_register(data):
    """Register a worker extension via Socket.IO.

    Using Socket.IO ensures request.sid is consistent for both registration
    and disconnect cleanup.

    Parameters
    ----------
    data : dict
        {
            "roomId": str,
            "name": str,
            "category": str,
            "schema": dict,
            "public": bool (optional, defaults to False)
        }
    """
    sid = request.sid
    r = current_app.extensions["redis"]

    try:
        room_id = data["roomId"]
        name = data["name"]
        category = data["category"]
        schema = data["schema"]
        public = data.get("public", False)  # Default to False if not specified
    except KeyError as e:
        return {"success": False, "error": f"Missing required field: {e}"}

    # Get user name from sid
    user_name = get_user_name_from_sid(sid)
    if not user_name:
        return {"success": False, "error": "Not authenticated"}

    # Check admin privileges for public extensions
    if public:
        try:
            # Get role from Redis (stored during socket connection)
            session_keys = SessionKeys(sid)
            role = r.get(session_keys.role())
            if role is None:
                log.error(f"Role not found for sid {sid}")
                return {"success": False, "error": "Failed to verify admin privileges"}

            if role != "admin":
                log.warning(
                    f"User {user_name} (role: {role}) attempted to register public extension '{name}'"
                )
                return {
                    "success": False,
                    "error": "Only admin users can register public extensions",
                }
        except Exception as e:
            log.error(f"Failed to check user role: {e}")
            return {"success": False, "error": "Failed to verify admin privileges"}

    # Use request.sid as the worker_id - this is consistent across registration and disconnect
    worker_id = sid

    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
        "analysis": analysis,
    }

    # Prevent overriding server-side extensions
    if category in category_map and name in category_map[category]:
        log.warning(
            f"Blocked attempt to register extension '{name}' in category '{category}' "
            f"- name conflicts with server-side extension (security violation)"
        )
        return {
            "success": False,
            "error": f"Cannot register extension '{name}': name is reserved for server-side extensions",
        }

    scope = "global" if public else f"room {room_id}"
    log.info(
        f"Registering {'global' if public else 'room-scoped'} extension: name={name}, category={category}, worker_id={worker_id}, scope={scope}"
    )

    # Use global or room-scoped keys based on public flag
    if public:
        keys = ExtensionKeys.for_global_extension(category, name)
        worker_extensions_key = ExtensionKeys.global_user_extensions_key(
            category, worker_id
        )
    else:
        keys = ExtensionKeys.for_extension(room_id, category, name)
        worker_extensions_key = ExtensionKeys.user_extensions_key(
            room_id, category, worker_id
        )

    existing_schema = r.hget(keys.schema, name)

    if existing_schema is not None:
        existing_schema = json.loads(existing_schema)
        if existing_schema != schema:
            return {
                "success": False,
                "error": "Extension with this name already exists with a different schema",
            }
        r.sadd(keys.idle_workers, worker_id)
        r.sadd(worker_extensions_key, name)

        log.info(
            f"Worker {worker_id} re-registered for extension '{name}' "
            f"in category '{category}', invalidating schema"
        )

        # For global extensions, broadcast to all rooms; for room-scoped, broadcast to specific room
        if public:
            # Broadcast to all connected clients (global extension)
            # Note: Flask-SocketIO broadcasts by not specifying a room
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
            )
        else:
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
                to=f"room:{room_id}",
            )

        return {
            "success": True,
            "workerId": worker_id,  # Return server's sid to client
            "message": "Extension already registered with same schema. Worker marked as idle.",
        }
    else:
        # Brand new extension
        with r.pipeline() as pipe:
            pipe.hset(keys.schema, name, json.dumps(schema))
            pipe.sadd(keys.idle_workers, worker_id)
            pipe.sadd(worker_extensions_key, name)
            pipe.execute()

        # For global extensions, broadcast to all rooms; for room-scoped, broadcast to specific room
        if public:
            # Broadcast to all connected clients (global extension)
            # Note: Flask-SocketIO broadcasts by not specifying a room
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
            )
        else:
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
                to=f"room:{room_id}",
            )

        return {"success": True, "workerId": worker_id}  # Return server's sid to client


@socketio.on("filesystem:register")
def handle_filesystem_register(data):
    """Register a filesystem via Socket.IO.

    Using Socket.IO ensures request.sid is consistent for both registration
    and disconnect cleanup.

    Parameters
    ----------
    data : dict
        {
            "roomId": str,
            "name": str,
            "fsType": str,
            "public": bool (optional, defaults to False)
        }
    """
    sid = request.sid
    r = current_app.extensions["redis"]

    try:
        room_id = data["roomId"]
        name = data["name"]
        fs_type = data["fsType"]
        public = data.get("public", False)
    except KeyError as e:
        return {"success": False, "error": f"Missing required field: {e}"}

    # Get user name from sid
    user_name = get_user_name_from_sid(sid)
    if not user_name:
        return {"success": False, "error": "Not authenticated"}

    # Check admin privileges for public filesystems
    if public:
        # Celery filesystem workers are automatically granted admin privileges
        is_celery_worker = user_name.startswith("celery-fs:")

        if not is_celery_worker:
            try:
                session_keys = SessionKeys(sid)
                role = r.get(session_keys.role())
                if role is None:
                    log.error(f"Role not found for sid {sid}")
                    return {"success": False, "error": "Failed to verify admin privileges"}

                if role != "admin":
                    log.warning(
                        f"User {user_name} (role: {role}) attempted to register public filesystem '{name}'"
                    )
                    return {
                        "success": False,
                        "error": "Only admin users can register public filesystems",
                    }
            except Exception as e:
                log.error(f"Failed to check user role: {e}")
                return {"success": False, "error": "Failed to verify admin privileges"}
        else:
            log.info(f"Celery filesystem worker {user_name} granted admin privileges automatically")

    # Use request.sid as the worker_id
    worker_id = sid

    scope = "global" if public else f"room {room_id}"
    log.info(
        f"Registering {'global' if public else 'room-scoped'} filesystem: "
        f"name={name}, type={fs_type}, worker_id={worker_id}, scope={scope}"
    )

    # Use global or room-scoped keys based on public flag
    if public:
        keys = FilesystemKeys.for_global_filesystem(name)
        worker_filesystems_key = FilesystemKeys.global_user_filesystems_key(worker_id)
    else:
        keys = FilesystemKeys.for_filesystem(room_id, name)
        worker_filesystems_key = FilesystemKeys.user_filesystems_key(room_id, worker_id)

    # Store filesystem metadata in Redis
    # Convert boolean to string for Redis compatibility
    fs_metadata = {
        "name": name,
        "fsType": fs_type,
        "workerId": worker_id,
        "userName": user_name,
        "public": "true" if public else "false",
    }

    with r.pipeline() as pipe:
        # Store filesystem metadata as a hash
        pipe.hset(keys.metadata, mapping=fs_metadata)
        # Store worker ID
        pipe.set(keys.worker, worker_id)
        # Add filesystem name to worker's filesystem set
        pipe.sadd(worker_filesystems_key, name)
        pipe.execute()

    log.info(f"Filesystem '{name}' registered successfully by worker {worker_id}")

    # Notify clients about filesystem availability change
    # Emit to room if room-scoped, broadcast globally if public
    if public:
        # Global filesystem - notify all clients
        socketio.emit(SocketEvents.FILESYSTEMS_UPDATE, {"scope": "global"})
    else:
        # Room-scoped filesystem - notify clients in that room
        socketio.emit(
            SocketEvents.FILESYSTEMS_UPDATE,
            {"scope": "room"},
            to=f"room:{room_id}",
        )

    return {"success": True, "workerId": worker_id}


@socketio.on("set_frame_atomic")
def handle_set_frame_atomic(data):
    """
    Handles a single frame jump. REJECTED if a presenter is active.
    """
    room = get_project_room_from_session(request.sid)
    room_keys = RoomKeys(room)
    redis_client = current_app.extensions["redis"]

    if redis_client.get(room_keys.presenter_lock()) not in [request.sid, None]:
        return {
            "success": False,
            "error": "LockError",
            "message": "Cannot set frame while presenter is active",
        }

    frame = data.get("frame")
    if frame is not None:
        try:
            # Validate and convert to int (Plotly may send float)
            frame_int = int(frame)
            if frame_int < 0:
                return {"success": False, "error": "Frame must be non-negative"}
            redis_client.set(room_keys.current_frame(), frame_int)
            emit(
                "frame_update",
                {"frame": frame_int},
                to=f"room:{room}",
                skip_sid=request.sid,
            )
            return {"success": True}
        except (ValueError, TypeError) as e:
            log.error(f"Invalid frame value: {frame} - {e}")
            return {"success": False, "error": f"Invalid frame value: {frame}"}

    return {"success": False, "error": "Frame parameter missing"}


@socketio.on("set_frame_continuous")
def handle_set_frame_continuous(data):
    """
    Handles continuous frame updates. REQUIRES sender to be the presenter.
    """
    room = get_project_room_from_session(request.sid)
    room_keys = RoomKeys(room)
    redis_client = current_app.extensions["redis"]

    presenter_sid = redis_client.get(room_keys.presenter_lock())

    if presenter_sid and presenter_sid == request.sid:
        frame = data.get("frame")
        if frame is not None:
            try:
                # Validate and convert to int (Plotly may send float)
                frame_int = int(frame)
                if frame_int < 0:
                    log.warning(f"Negative frame rejected: {frame_int}")
                    return {"success": False, "error": "Frame must be non-negative"}
                redis_client.set(room_keys.current_frame(), frame_int)
                emit(
                    "frame_update",
                    {"frame": frame_int},
                    to=f"room:{room}",
                    skip_sid=request.sid,
                )
                return {"success": True}
            except (ValueError, TypeError) as e:
                log.error(f"Invalid frame value in continuous: {frame} - {e}")
                return {"success": False, "error": f"Invalid frame value: {frame}"}


@socketio.on("frame_selection:set")
def handle_frame_selection_set(data):
    """
    Handles setting the frame selection for a room.
    """
    room = get_project_room_from_session(request.sid)
    if not room:
        return {"success": False, "error": "Client has not joined a room"}

    redis_client = current_app.extensions["redis"]
    indices = data.get("indices", [])

    if not isinstance(indices, list):
        return {"success": False, "error": "indices must be a list"}

    if not all(isinstance(idx, int) and idx >= 0 for idx in indices):
        return {"success": False, "error": "All indices must be non-negative integers"}

    # Store frame selection in Redis
    room_keys = RoomKeys(room)
    redis_client.set(room_keys.frame_selection(), json.dumps(indices))

    # Emit update to all clients in the room
    emit(
        "frame_selection:update",
        {"indices": indices},
        to=f"room:{room}",
        skip_sid=request.sid,
    )

    return {"success": True}


@socketio.on("request_presenter_token")
def handle_request_presenter_token():
    sid = request.sid
    room = get_project_room_from_session(sid)
    print(f"Presenter token requested by {sid} in room {room}")
    if not room:
        return {"success": False, "reason": "Not in a valid room"}

    room_keys = RoomKeys(room)
    r = current_app.extensions["redis"]

    # --- UPDATED LOGIC ---
    # Get the current holder of the lock
    current_holder = r.get(room_keys.presenter_lock())

    # Case 1: No one has the lock, or the requester already has it (renewal)
    if current_holder is None or current_holder == sid:
        # Set (or reset) the lock with the new expiry
        r.set(room_keys.presenter_lock(), sid, ex=TOKEN_EXPIRY_SECONDS)

        # If this is a brand new presenter, inform the room
        if current_holder is None:
            from zndraw.app.room_manager import emit_room_update

            emit_room_update(socketio, room, skip_sid=sid, presenterSid=sid)

        return {"success": True}
    else:
        # Case 2: Someone else has the lock
        return {"success": False, "reason": "Presenter lock is held by another user"}


@socketio.on("release_presenter_token")
def handle_release_presenter_token():
    room = get_project_room_from_session(request.sid)
    room_keys = RoomKeys(room)
    print(f"Presenter token release requested by {request.sid} in room {room}")
    r = current_app.extensions["redis"]

    presenter_sid = r.get(room_keys.presenter_lock())

    if presenter_sid and presenter_sid == request.sid:
        r.delete(room_keys.presenter_lock())
        from zndraw.app.room_manager import emit_room_update

        emit_room_update(socketio, room, skip_sid=request.sid, presenterSid=None)
        return {"success": True}
    else:
        return {"success": False, "error": "Not the current presenter"}


@socketio.on("lock:acquire")
def acquire_lock(data):
    sid = request.sid
    r = current_app.extensions["redis"]
    target = data.get("target")
    ttl = data.get("ttl", 60)  # Default to 60 seconds if not specified
    room = get_project_room_from_session(sid)
    user_name = get_user_name_from_sid(sid)

    if not room or not target or not user_name:
        return {"success": False, "error": "Room, target, or userName missing"}

    # Validate TTL - must not exceed 300 seconds (5 minutes)
    if not isinstance(ttl, (int, float)) or ttl <= 0:
        return {"success": False, "error": "TTL must be a positive number"}
    if ttl > 300:
        return {"success": False, "error": "TTL cannot exceed 300 seconds (5 minutes)"}

    lock_key = get_lock_key(room, target)
    # Store userName in lock (not sid) so HTTP endpoints can verify
    if r.set(lock_key, user_name, nx=True, ex=int(ttl)):
        log.debug(
            f"Lock acquired for '{target}' in room '{room}' by user {user_name} (sid:{sid}) with TTL {ttl}s"
        )

        # Broadcast lock acquisition to all clients
        # Only broadcast for trajectory:meta locks (the one used by vis.lock())
        if target == "trajectory:meta":
            lock_metadata = LockMetadata(
                msg=None,  # No message yet (will be updated by lock:msg if provided)
                userName=user_name,
                timestamp=datetime.datetime.utcnow().isoformat(),
            )
            emit_room_update(
                socketio, room, skip_sid=sid, metadataLocked=lock_metadata.model_dump()
            )

        return {"success": True}
    else:
        lock_holder = r.get(lock_key)
        log.info(
            f"Lock for '{target}' in room '{room}' already held by {lock_holder}, denied for {user_name} (sid:{sid})"
        )
        return {"success": False}


@socketio.on("lock:release")
def release_lock(data):
    sid = request.sid
    r = current_app.extensions["redis"]
    target = data.get("target")
    room = get_project_room_from_session(sid)
    user_name = get_user_name_from_sid(sid)

    if not room or not target or not user_name:
        return {"success": False, "error": "Room, target, or userName missing"}

    lock_key = get_lock_key(room, target)
    lock_holder = r.get(lock_key)
    # Compare with userName (not sid) since that's what we store
    if lock_holder == user_name:
        # Delete lock AND metadata
        r.delete(lock_key)
        r.delete(f"{lock_key}:metadata")

        # Broadcast lock release via room:update event
        emit_room_update(socketio, room, skip_sid=sid, metadataLocked=None)

        log.debug(
            f"Lock released for '{target}' in room '{room}' by user {user_name} (sid:{sid})"
        )
        return {"success": True}

    log.warning(
        f"Failed release: Lock for '{target}' in room '{room}' held by {lock_holder}, not by {user_name} (sid:{sid})"
    )
    return {"success": False}


@socketio.on("lock:refresh")
def refresh_lock(data):
    """Refresh the TTL of an existing lock to prevent expiration during long operations."""
    sid = request.sid
    r = current_app.extensions["redis"]
    target = data.get("target")
    ttl = data.get("ttl", 60)  # Default to 60 seconds if not specified
    room = get_project_room_from_session(sid)
    user_name = get_user_name_from_sid(sid)

    if not room or not target or not user_name:
        return {"success": False, "error": "Room, target, or userName missing"}

    # Validate TTL - must not exceed 300 seconds (5 minutes)
    if not isinstance(ttl, (int, float)) or ttl <= 0:
        return {"success": False, "error": "TTL must be a positive number"}
    if ttl > 300:
        return {"success": False, "error": "TTL cannot exceed 300 seconds (5 minutes)"}

    lock_key = get_lock_key(room, target)
    lock_holder = r.get(lock_key)

    # Only refresh if the lock is held by this client (compare with userName)
    if lock_holder == user_name:
        # Reset the TTL
        r.expire(lock_key, int(ttl))
        log.debug(
            f"Lock refreshed for '{target}' in room '{room}' by user {user_name} (sid:{sid}) with TTL {ttl}s"
        )
        return {"success": True}

    log.warning(
        f"Failed refresh: Lock for '{target}' in room '{room}' not held by {sid}"
    )
    return {"success": False}


@socketio.on("lock:msg")
def update_lock_message(data):
    """Update metadata for an active lock.

    This endpoint allows the lock holder to send descriptive metadata
    about what they're doing with the lock (e.g., "Uploading 1000 frames").

    IMPORTANT: Server validates that the requesting client actually holds the lock.

    Parameters
    ----------
    data : dict
        Contains 'target' and 'metadata' fields

    Returns
    -------
    dict
        Success response or error
    """
    sid = request.sid
    r = current_app.extensions["redis"]
    target = data.get("target")
    metadata = data.get("metadata", {})
    room = get_project_room_from_session(sid)
    user_name = get_user_name_from_sid(sid)

    if not room or not target or not user_name:
        return {"success": False, "error": "Missing required fields"}

    lock_key = get_lock_key(room, target)
    lock_holder = r.get(lock_key)

    # Validate that this user actually holds the lock
    if lock_holder != user_name:
        log.warning(
            f"Rejected lock:msg from {user_name}: lock for '{target}' held by {lock_holder}"
        )
        return {"success": False, "error": "Lock not held by this user"}

    # Store metadata with same TTL as lock
    metadata_key = f"{lock_key}:metadata"

    # Store metadata
    metadata_with_user = {
        "userName": user_name,
        "timestamp": datetime.datetime.utcnow().isoformat(),
        **metadata,
    }

    # Serialize to JSON for storage
    r.set(metadata_key, json.dumps(metadata_with_user))

    # Match lock TTL
    lock_ttl = r.ttl(lock_key)
    if lock_ttl > 0:
        r.expire(metadata_key, lock_ttl)

    # Broadcast to room via room:update event
    lock_metadata = LockMetadata(
        msg=metadata.get("msg"),
        userName=user_name,
        timestamp=metadata_with_user["timestamp"],
    )
    emit_room_update(
        socketio, room, skip_sid=sid, metadataLocked=lock_metadata.model_dump()
    )

    log.debug(f"Lock metadata updated for '{target}' by {user_name}: {metadata}")
    return {"success": True}


@socketio.on("join:overview")
def handle_join_overview():
    """Client joining /rooms page - join overview:public room."""
    from flask_socketio import join_room

    sid = request.sid
    join_room("overview:public")
    log.debug(f"Client {sid} joined overview:public")
    return {"status": "joined", "room": "overview:public"}


@socketio.on("leave:overview")
def handle_leave_overview():
    """Client leaving /rooms page - leave overview:public room."""
    from flask_socketio import leave_room

    sid = request.sid
    leave_room("overview:public")
    log.debug(f"Client {sid} left overview:public")
    return {"status": "left", "room": "overview:public"}


@socketio.on("join:room")
def handle_join_room(data):
    """Client joining specific room page - join room:<room_id> and leave overview:public."""
    from flask_socketio import join_room, leave_room

    sid = request.sid
    room_id = data.get("roomId")

    if not room_id:
        return {"status": "error", "message": "roomId required"}

    # Get userName from sid
    r = current_app.extensions["redis"]
    user_name = get_user_name_from_sid(sid)

    if not user_name:
        log.error(f"Cannot join room: userName not found for sid {sid}")
        return {"status": "error", "message": "User not found"}

    # Leave overview if joined
    leave_room("overview:public")

    # Join specific room (Flask-SocketIO level)
    join_room(f"room:{room_id}")

    # Update Redis to track which room this user is in
    user_keys = UserKeys(user_name)
    r.hset(user_keys.hash_key(), "currentRoom", room_id)

    log.debug(f"User {sid} (userName: {user_name}) joined room:{room_id}")
    return {"status": "joined", "room": f"room:{room_id}"}


@socketio.on("leave:room")
def handle_leave_room(data):
    """Client leaving specific room page - leave room:<room_id>."""
    from flask_socketio import leave_room

    sid = request.sid
    room_id = data.get("roomId")

    if not room_id:
        return {"status": "error", "message": "roomId required"}

    leave_room(f"room:{room_id}")
    log.debug(f"Client {sid} left room:{room_id}")
    return {"status": "left", "room": f"room:{room_id}"}


@socketio.on("chat:message:create")
def handle_chat_message_create(data):
    """
    Create a new chat message.
    Payload: { "content": "message text" }
    Returns: { "success": bool, "message": Message | None, "error": str | None }
    """
    from .chat_utils import create_message

    sid = request.sid
    r = current_app.extensions["redis"]
    room = get_project_room_from_session(sid)

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    content = data.get("content")
    if not content or not isinstance(content, str):
        return {"success": False, "error": "Message content is required"}

    # Get userName from session using new schema: sid -> userName
    user_name = get_user_name_from_sid(sid)
    if not user_name:
        return {"success": False, "error": "User not found"}

    try:
        # Create message using helper function
        message = create_message(r, room, user_name, content)

        # Emit to room (excluding sender)
        emit("chat:message:new", message, to=f"room:{room}", include_self=True)

        return {"success": True, "message": message}
    except Exception as e:
        log.error(f"Failed to create chat message: {e}")
        return {"success": False, "error": str(e)}


@socketio.on("chat:message:edit")
def handle_chat_message_edit(data):
    """
    Edit an existing chat message.
    Payload: { "messageId": "msg_room_42", "content": "new text" }
    Returns: { "success": bool, "message": Message | None, "error": str | None }
    """
    from .chat_utils import get_message, update_message

    sid = request.sid
    r = current_app.extensions["redis"]
    room = get_project_room_from_session(sid)

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    message_id = data.get("messageId")
    content = data.get("content")

    if not message_id or not isinstance(message_id, str):
        return {"success": False, "error": "Message ID is required"}

    if not content or not isinstance(content, str):
        return {"success": False, "error": "Message content is required"}

    # Get userName from session using new schema: sid -> userName
    user_name = get_user_name_from_sid(sid)
    if not user_name:
        return {"success": False, "error": "User not found"}

    try:
        # Fetch existing message
        existing_message = get_message(r, room, message_id)
        if not existing_message:
            return {"success": False, "error": "Message not found"}

        # Authorization check: verify user owns the message
        if existing_message["author"]["id"] != user_name:
            return {
                "success": False,
                "error": "You can only edit your own messages",
            }

        # Update message
        updated_message = update_message(r, room, message_id, content)

        # Emit to room
        emit("chat:message:updated", updated_message, to=f"room:{room}")

        return {"success": True, "message": updated_message}
    except Exception as e:
        log.error(f"Failed to edit chat message: {e}")
        return {"success": False, "error": str(e)}
