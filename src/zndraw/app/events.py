import datetime
import json
import logging
import time
import typing as t
import uuid

from flask import current_app, request
from flask_socketio import emit

from zndraw.analytics import connected_users
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
    if "sessionId" not in auth:
        log.critical(f"Client {sid} connected without sessionId in auth")
    else:
        session_id = auth.get("sessionId")
        if session_id is None:
            log.critical(f"Client {sid} connected with null sessionId")
        else:
            # Bidirectional mapping for session cleanup:
            # session_id → sid (used for skip_sid in emit)
            r.set(SessionKeys.session_to_sid(session_id), sid)
            # sid → session_id (used for cleanup on disconnect)
            session_keys = SessionKeys(sid)
            r.set(session_keys.session_id(), session_id)
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

    # Increment connected users metric
    connected_users.inc()

    return {"status": "ok", "userName": user_name}


@socketio.on("disconnect")
def handle_disconnect(*args, **kwargs):
    """Handle client disconnect.

    Note: No parameters needed - Flask-SocketIO provides request.sid automatically.
    The framework may pass arguments but we don't use them.
    """
    sid = request.sid
    r = current_app.extensions["redis"]

    # Get userName from connection lookup
    session_keys = SessionKeys(sid)
    user_name = r.get(session_keys.username())

    # Decrement connected users metric
    connected_users.dec()

    if not user_name:
        log.info(f"Client disconnected: {sid} (no userName found)")
        return

    # Get user data
    user_keys = UserKeys(user_name)
    room_name = r.hget(user_keys.hash_key(), "currentRoom")

    log.info(f"User disconnected: sid={sid}, user={user_name}, room={room_name}")

    # Clean up session→sid bidirectional mapping (prevents memory leak)
    session_id = r.get(session_keys.session_id())
    if session_id:
        # Remove session_id→sid mapping
        r.delete(SessionKeys.session_to_sid(session_id))
        # Remove sid→session_id mapping
        r.delete(session_keys.session_id())

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

    # --- Lock Cleanup Logic ---
    # Scan for locks held by this session
    lock_keys = r.scan_iter("*:lock:*")
    for key in lock_keys:
        # Skip metadata keys
        if key.endswith(":metadata"):
            continue

        lock_data_str = r.get(key)
        if not lock_data_str:
            continue

        try:
            lock_data = json.loads(lock_data_str)
            lock_session_id = lock_data.get("sessionId")

            # If this lock is held by the disconnecting session, clean it up
            if lock_session_id == sid:
                log.warning(
                    f"Cleaning up orphaned lock '{key}' held by disconnected session {sid}"
                )
                r.delete(key)
                # Also delete associated metadata if it exists
                metadata_key = f"{key}:metadata" if isinstance(key, str) else key + b":metadata"
                r.delete(metadata_key)
        except (json.JSONDecodeError, TypeError, AttributeError):
            # If lock data is not JSON (old format), skip it
            continue

    if room_name:
        room_keys = RoomKeys(room_name)
        presenter_sid = r.get(room_keys.presenter_lock())
        if presenter_sid and presenter_sid == sid:
            r.delete(room_keys.presenter_lock())
            # Inform everyone that the presenter left via room:update
            from zndraw.app.room_manager import emit_room_update

            emit_room_update(socketio, room_name, skip_sid=sid, presenterSid=None)

    # --- Job Cleanup Logic ---
    # Handle jobs assigned to or being processed by this worker
    if not sid:
        log.error(f"No sid during disconnect cleanup")
        return

    worker_id = sid  # Workers are tracked by server's socket sid

    # Get all jobs assigned to this worker using the reverse mapping
    from .redis_keys import WorkerKeys
    worker_keys = WorkerKeys(worker_id)
    worker_job_ids = list(r.smembers(worker_keys.active_jobs()))

    if worker_job_ids:
        log.warning(
            f"Worker {worker_id} disconnected with {len(worker_job_ids)} active job(s): {worker_job_ids}"
        )

        from .job_manager import JobManager, JobStatus

        for job_id in worker_job_ids:
            try:
                # Get job details
                job_data = JobManager.get_job(r, job_id)
                if not job_data:
                    log.error(f"Job {job_id} not found during disconnect cleanup")
                    r.srem(worker_keys.active_jobs(), job_id)
                    continue

                category = job_data.get("category")
                extension = job_data.get("extension")
                job_room = job_data.get("room")
                current_status = job_data.get("status")

                log.info(
                    f"Failing job {job_id} ({category}/{extension} in room {job_room}, status: {current_status})"
                )

                # Only fail jobs that are assigned or processing (not already completed/failed)
                if current_status in [JobStatus.ASSIGNED, JobStatus.PROCESSING]:
                    # Fail the job (automatically emits job:state_changed with error metadata)
                    JobManager.fail_job(
                        r,
                        job_id,
                        f"Worker {worker_id} disconnected while processing job",
                        socketio=socketio
                    )
                    log.info(f"Marked job {job_id} as failed due to worker disconnect")
                else:
                    log.info(f"Job {job_id} has status {current_status}, not failing")

                # Job will be removed from worker set by fail_job()

            except Exception as e:
                log.error(
                    f"Failed to fail job {job_id} after worker disconnect: {e}"
                )
                # Clean up manually if fail_job errored
                r.srem(worker_keys.active_jobs(), job_id)

        # Clean up worker's job set (should be empty after fail_job removes entries)
        r.delete(worker_keys.active_jobs())
        log.info(f"Cleaned up worker job set for {worker_id}")

    # Clean up worker capacity key
    capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
    r.delete(capacity_key)
    log.info(f"Cleaned up worker capacity for {worker_id}")

    # Extension cleanup - workers are tracked by server's sid
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

                # Remove the worker_id from the workers registry
                pipe.hdel(keys.workers, worker_id)

                # Check how many workers remain
                pipe.hlen(keys.workers)

            # Each extension now produces 2 results in the pipeline
            results = pipe.execute()

        # Iterate through the results to decide which extensions to delete
        for i, ext_name in enumerate(worker_extensions):
            # Get the hlen result for this extension
            total_remaining = results[i * 2 + 1]

            log.info(
                f"Extension '{ext_name}': {total_remaining} workers remaining after removing {worker_id}."
            )

            # Only delete extension if no workers AND no pending jobs
            if total_remaining == 0:
                keys = ExtensionKeys.for_extension(room_name, category, ext_name)
                pending_jobs_count = r.zcard(keys.pending_jobs)

                if pending_jobs_count == 0:
                    extensions_to_delete.append(ext_name)
                    log.info(
                        f"Extension '{ext_name}' marked for deletion: no workers, no pending jobs"
                    )
                else:
                    log.info(
                        f"Extension '{ext_name}' kept despite no workers: {pending_jobs_count} pending jobs"
                    )

        # If any extensions are now orphaned, delete them and their state sets
        if extensions_to_delete:
            log.info(
                f"Deleting orphaned extensions in '{category}': {extensions_to_delete}"
            )
            with r.pipeline() as pipe:
                for ext_name in extensions_to_delete:
                    keys = ExtensionKeys.for_extension(room_name, category, ext_name)
                    # Delete the workers registry
                    pipe.delete(keys.workers)
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

                # Remove the worker_id from the workers registry
                pipe.hdel(keys.workers, worker_id)

                # Check how many workers remain
                pipe.hlen(keys.workers)

            # Each extension produces 2 results
            results = pipe.execute()

        # Iterate through results to decide which extensions to delete
        for i, ext_name in enumerate(global_worker_extensions):
            total_remaining = results[i * 2 + 1]

            log.info(
                f"Global extension '{ext_name}': {total_remaining} workers remaining after removing {worker_id}."
            )

            # Only delete if no workers AND no pending jobs
            if total_remaining == 0:
                keys = ExtensionKeys.for_global_extension(category, ext_name)
                pending_jobs_count = r.zcard(keys.pending_jobs)

                if pending_jobs_count == 0:
                    global_extensions_to_delete.append(ext_name)
                    log.info(
                        f"Global extension '{ext_name}' marked for deletion: no workers, no pending jobs"
                    )
                else:
                    log.info(
                        f"Global extension '{ext_name}' kept despite no workers: {pending_jobs_count} pending jobs"
                    )

        # Delete orphaned global extensions
        if global_extensions_to_delete:
            log.info(
                f"Deleting orphaned global extensions in '{category}': {global_extensions_to_delete}"
            )
            with r.pipeline() as pipe:
                for ext_name in global_extensions_to_delete:
                    keys = ExtensionKeys.for_global_extension(category, ext_name)
                    # Delete the workers registry
                    pipe.delete(keys.workers)
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

    # Clean up connection lookup - MUST BE LAST to prevent double-disconnect issues
    # If this handler is called twice, the second call will find no username and return early
    r.delete(session_keys.username())
    log.info(f"Cleaned up connection lookup for session {sid}")


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

    # Send current progress trackers to joining client
    room_keys = RoomKeys(room_id)
    progress_data = r.hgetall(room_keys.progress())

    # Convert Redis data to client format
    progress_trackers = {}
    for progress_id, progress_json in progress_data.items():
        progress_dict = json.loads(progress_json)
        progress_trackers[progress_id] = progress_dict

    # Emit progress trackers to joining client only
    if progress_trackers:
        emit("progress:initial", {"progressTrackers": progress_trackers}, to=sid)

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
        emit(SocketEvents.CHAT_MESSAGE_NEW, message, to=f"room:{room}", include_self=True)

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
        emit(SocketEvents.CHAT_MESSAGE_UPDATED, updated_message, to=f"room:{room}")

        return {"success": True, "message": updated_message}
    except Exception as e:
        log.error(f"Failed to edit chat message: {e}")
        return {"success": False, "error": str(e)}

