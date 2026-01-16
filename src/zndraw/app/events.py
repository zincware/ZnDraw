import json
import logging
import typing as t

from flask import current_app, request
from flask_socketio import emit

from zndraw.analytics import connected_users
from zndraw.geometries.camera import Camera
from zndraw.server import socketio
from zndraw.utils.time import utc_now_iso

from .constants import SocketEvents
from .redis_keys import ExtensionKeys, FilesystemKeys, RoomKeys, SessionKeys, UserKeys

log = logging.getLogger(__name__)

# for crash handling.
TOKEN_EXPIRY_SECONDS = 10


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


def get_session_camera_key(session_id: str) -> str:
    """Get the geometry key for a session camera."""
    return f"cam:session:{session_id}"


def create_session_camera(r, room: str, session_id: str) -> None:
    """Create a Camera geometry for a session.

    The session camera is the user's viewport into the scene.
    It's stored in geometries with helper_visible=False (not cluttering the scene).
    """
    camera_key = get_session_camera_key(session_id)
    keys = RoomKeys(room)

    # Check if camera already exists (reconnection case)
    existing = r.hget(keys.geometries(), camera_key)
    if existing:
        log.debug(f"Session camera {camera_key} already exists, skipping creation")
        return

    # Create camera with default position, helper_visible=False, protected=True
    camera = Camera(helper_visible=False, protected=True)
    geometry_data = {"type": "Camera", "data": camera.model_dump()}

    r.hset(keys.geometries(), camera_key, json.dumps(geometry_data))
    log.debug(f"Created session camera: {camera_key}")

    # Emit geometry invalidation so clients know about the new camera
    socketio.emit(
        SocketEvents.INVALIDATE_GEOMETRY,
        {"key": camera_key, "operation": "set"},
        to=f"room:{room}",
    )


def delete_session_camera(r, room: str, session_id: str) -> None:
    """Delete a session's Camera geometry on disconnect."""
    camera_key = get_session_camera_key(session_id)
    keys = RoomKeys(room)

    result = r.hdel(keys.geometries(), camera_key)
    if result > 0:
        log.debug(f"Deleted session camera: {camera_key}")
        # Emit geometry invalidation
        socketio.emit(
            SocketEvents.INVALIDATE_GEOMETRY,
            {"key": camera_key, "operation": "delete"},
            to=f"room:{room}",
        )
    else:
        log.debug(f"Session camera {camera_key} not found (already deleted)")


@socketio.on("connect")
def handle_connect(auth):
    """Handle socket connection with JWT authentication."""
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

    # Validate user exists (must be created via /api/user/register)
    # Users are ONLY created there - socket connect never creates users
    user_service = current_app.extensions["user_service"]
    if not user_service.username_exists(user_name):
        log.warning(f"Socket connect failed: user '{user_name}' not found")
        raise ConnectionRefusedError("User not found. Please register first.")

    # Update user's current SID
    user_keys = UserKeys(user_name)
    r.hset(
        user_keys.hash_key(),
        mapping={
            "currentSid": sid,
            "lastActivity": utc_now_iso(),
        },
    )
    # Register SID → userName mapping and role mapping
    session_keys = SessionKeys(sid)
    r.set(session_keys.username(), user_name)
    r.set(session_keys.role(), user_role)  # Store role for this session

    # Join user-specific room for user-targeted events (e.g., extension notifications)
    join_room(f"user:{user_name}")
    log.debug(f"User {user_name} connected (sid: {sid})")

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
        log.debug(f"Client disconnected: {sid} (no userName found)")
        return

    # Get session_id and room from session data (before cleanup)
    session_id = r.get(session_keys.session_id())
    room_name = None
    if session_id:
        session_data_raw = r.get(SessionKeys.session_data(session_id))
        if session_data_raw:
            session_data = json.loads(session_data_raw)
            room_name = session_data.get("roomId")

    user_keys = UserKeys(user_name)
    log.debug(f"User disconnected: sid={sid}, user={user_name}, room={room_name}")
    if session_id:
        # Remove session_id→sid mapping
        r.delete(SessionKeys.session_to_sid(session_id))
        # Remove sid→session_id mapping
        r.delete(session_keys.session_id())

        # Clean up frontend session data if in a room
        if room_name:
            room_keys = RoomKeys(room_name)

            # Remove from frontend sessions set
            r.srem(room_keys.frontend_sessions(), session_id)

            # Delete session camera geometry
            delete_session_camera(r, room_name, session_id)

            # Delete session settings
            r.delete(room_keys.session_settings(session_id))

            # Delete active camera key
            r.delete(room_keys.session_active_camera(session_id))

            log.debug(f"Cleaned up frontend session data for {session_id}")

    # Update user's currentSid to empty (user still exists but disconnected)
    r.hset(user_keys.hash_key(), "currentSid", "")

    # Note: We don't remove the user from the room or delete user data
    # The user may reconnect and rejoin the same room
    # Only when a user explicitly joins a different room do we remove them from the old room

    if not room_name:
        log.debug(f"User {user_name} disconnected (was not in a room)")

    # --- Lock Cleanup Logic ---
    if session_id:
        from .route_utils import emit_lock_update

        session_locks_key = SessionKeys.session_locks(session_id)
        lock_keys = r.smembers(session_locks_key)

        if lock_keys:
            log.debug(
                f"Cleaning up {len(lock_keys)} orphaned lock(s) for session {session_id}"
            )
            for lock_key in lock_keys:
                # Delete the lock and its metadata
                r.delete(lock_key)
                metadata_key = (
                    f"{lock_key}:metadata"
                    if isinstance(lock_key, str)
                    else lock_key + b":metadata"
                )
                r.delete(metadata_key)

                # Parse lock_key to extract room_id and target
                # Format: room:{room_id}:lock:{target}
                try:
                    parts = (
                        lock_key.split(":")
                        if isinstance(lock_key, str)
                        else lock_key.decode().split(":")
                    )
                    if len(parts) >= 4 and parts[0] == "room" and parts[2] == "lock":
                        lock_room_id = parts[1]
                        lock_target = ":".join(
                            parts[3:]
                        )  # Handle targets like "trajectory:meta"

                        # Emit lock release event so frontends update their UI
                        emit_lock_update(
                            room_id=lock_room_id,
                            target=lock_target,
                            action="released",
                            user_name=None,
                            message=None,
                            timestamp=None,
                            session_id=session_id,
                        )
                        log.debug(
                            f"Emitted lock:update (released) for '{lock_target}' in room '{lock_room_id}'"
                        )
                except Exception as e:
                    log.warning(
                        f"Failed to parse/emit lock release for '{lock_key}': {e}"
                    )

                log.debug(f"Cleaned up orphaned lock '{lock_key}'")

            # Delete the session locks set itself
            r.delete(session_locks_key)

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
        log.error("No sid during disconnect cleanup")
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

                log.debug(
                    f"Failing job {job_id} ({category}/{extension} in room {job_room}, status: {current_status})"
                )

                # Only fail jobs that are assigned or processing (not already completed/failed)
                if current_status in [JobStatus.ASSIGNED, JobStatus.PROCESSING]:
                    # Fail the job (automatically emits job:update with error metadata)
                    JobManager.fail_job(
                        r,
                        job_id,
                        f"Worker {worker_id} disconnected while processing job",
                        socketio=socketio,
                    )
                    log.debug(f"Marked job {job_id} as failed due to worker disconnect")
                else:
                    log.debug(f"Job {job_id} has status {current_status}, not failing")

                # Job will be removed from worker set by fail_job()

            except Exception as e:
                log.error(f"Failed to fail job {job_id} after worker disconnect: {e}")
                # Clean up manually if fail_job errored
                r.srem(worker_keys.active_jobs(), job_id)

        # Clean up worker's job set (should be empty after fail_job removes entries)
        r.delete(worker_keys.active_jobs())
        log.debug(f"Cleaned up worker job set for {worker_id}")

    # Clean up worker capacity key
    capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
    r.delete(capacity_key)
    log.debug(f"Cleaned up worker capacity for {worker_id}")

    # Extension cleanup - workers are tracked by server's sid
    extension_categories = ["modifiers", "selections", "analysis"]
    log.debug(
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

        log.debug(
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

            log.debug(
                f"Extension '{ext_name}': {total_remaining} workers remaining after removing {worker_id}."
            )

            # Only delete extension if no workers AND no pending jobs
            if total_remaining == 0:
                keys = ExtensionKeys.for_extension(room_name, category, ext_name)
                pending_jobs_count = r.zcard(keys.pending_jobs)

                if pending_jobs_count == 0:
                    extensions_to_delete.append(ext_name)
                    log.debug(
                        f"Extension '{ext_name}' marked for deletion: no workers, no pending jobs"
                    )
                else:
                    log.debug(
                        f"Extension '{ext_name}' kept despite no workers: {pending_jobs_count} pending jobs"
                    )

        # If any extensions are now orphaned, delete them and their state sets
        if extensions_to_delete:
            log.debug(
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
        log.debug(f"Cleaned up worker-specific extension list: {worker_extensions_key}")

        # Notify clients about worker count changes
        # We always invalidate if this worker had any extensions, not just when deleting
        if worker_extensions:
            log.debug(
                f"Invalidating schema for category '{category}' in room '{room_name}' "
                "due to worker disconnect."
            )
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
                to=f"room:{room_name}",
            )

    # Clean up global (public) extensions
    log.debug(f"Checking for global extensions from worker_id={worker_id}...")
    for category in extension_categories:
        # Check if this worker registered any global extensions
        global_worker_extensions_key = ExtensionKeys.global_user_extensions_key(
            category, worker_id
        )
        global_worker_extensions = list(r.smembers(global_worker_extensions_key))

        if not global_worker_extensions:
            continue

        log.debug(
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

            log.debug(
                f"Global extension '{ext_name}': {total_remaining} workers remaining after removing {worker_id}."
            )

            # Only delete if no workers AND no pending jobs
            if total_remaining == 0:
                keys = ExtensionKeys.for_global_extension(category, ext_name)
                pending_jobs_count = r.zcard(keys.pending_jobs)

                if pending_jobs_count == 0:
                    global_extensions_to_delete.append(ext_name)
                    log.debug(
                        f"Global extension '{ext_name}' marked for deletion: no workers, no pending jobs"
                    )
                else:
                    log.debug(
                        f"Global extension '{ext_name}' kept despite no workers: {pending_jobs_count} pending jobs"
                    )

        # Delete orphaned global extensions
        if global_extensions_to_delete:
            log.debug(
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
        log.debug(
            f"Cleaned up global worker-specific extension list: {global_worker_extensions_key}"
        )

        # Notify ALL clients about global extension removal (broadcast globally)
        if global_worker_extensions:
            log.debug(
                f"Broadcasting global schema invalidation for category '{category}' due to worker disconnect"
            )
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category, "global": True},
            )

    # Filesystem cleanup - similar to extensions
    if room_name:
        worker_filesystems_key = FilesystemKeys.user_filesystems_key(
            room_name, worker_id
        )
        worker_filesystems = list(r.smembers(worker_filesystems_key))

        if worker_filesystems:
            log.debug(
                f"Worker {worker_id} was providing filesystems: {worker_filesystems}"
            )
            for fs_name in worker_filesystems:
                keys = FilesystemKeys.for_filesystem(room_name, fs_name)
                # Delete filesystem metadata and worker reference
                r.delete(keys.metadata)
                r.delete(keys.worker)
                log.debug(
                    f"Deleted room-scoped filesystem '{fs_name}' for worker {worker_id}"
                )

            # Clean up worker-specific filesystem list
            r.delete(worker_filesystems_key)
            log.debug(
                f"Cleaned up worker-specific filesystem list: {worker_filesystems_key}"
            )

            # Notify clients in the room that filesystems changed
            socketio.emit(
                SocketEvents.FILESYSTEMS_UPDATE,
                {"scope": "room"},
                to=f"room:{room_name}",
            )

    # Clean up global (public) filesystems
    global_worker_filesystems_key = FilesystemKeys.global_user_filesystems_key(
        worker_id
    )
    global_worker_filesystems = list(r.smembers(global_worker_filesystems_key))

    if global_worker_filesystems:
        log.debug(
            f"Worker {worker_id} was providing global filesystems: {global_worker_filesystems}"
        )
        for fs_name in global_worker_filesystems:
            keys = FilesystemKeys.for_global_filesystem(fs_name)
            # Delete filesystem metadata and worker reference
            r.delete(keys.metadata)
            r.delete(keys.worker)
            log.debug(f"Deleted global filesystem '{fs_name}' for worker {worker_id}")

        # Clean up worker-specific filesystem list
        r.delete(global_worker_filesystems_key)
        log.debug(
            f"Cleaned up global worker-specific filesystem list: {global_worker_filesystems_key}"
        )

        # Notify all clients that global filesystems changed
        socketio.emit(SocketEvents.FILESYSTEMS_UPDATE, {"scope": "global"})

    # Clean up connection lookup - MUST BE LAST to prevent double-disconnect issues
    # If this handler is called twice, the second call will find no username and return early
    log.debug("Cleaning up session username mapping: %s", session_keys.username())
    r.delete(session_keys.username())


@socketio.on("overview:join")
def handle_join_overview():
    """Client joining /rooms page - join overview:public room."""
    from flask_socketio import join_room

    sid = request.sid
    join_room("overview:public")
    log.debug(f"Client {sid} joined overview:public")
    return {"status": "joined", "room": "overview:public"}


@socketio.on("overview:leave")
def handle_leave_overview():
    """Client leaving /rooms page - leave overview:public room."""
    from flask_socketio import leave_room

    sid = request.sid
    leave_room("overview:public")
    log.debug(f"Client {sid} left overview:public")
    return {"status": "left", "room": "overview:public"}


@socketio.on("room:join")
def handle_room_join(data):
    """Join room and get minimal initialization data.

    This is the primary entry point for joining a room. It creates a session,
    joins socket rooms, and returns minimal state. Additional data (geometries,
    settings, selections, etc.) is fetched via REST endpoints.

    Room creation is handled separately via POST /api/rooms.

    Request
    -------
    {
        "roomId": str,  # Room to join
        "clientType": str  # "frontend" or "python"
    }

    Response (success)
    ------------------
    {
        "status": "ok",
        "sessionId": str,
        "cameraKey": str,  # cam:session:<sessionId>
        "step": int,       # Current frame index
        "frameCount": int, # Total number of frames
        "locked": bool     # Whether room is locked
    }

    Response (error)
    ----------------
    {
        "status": "error",
        "code": 404,
        "message": "Room not found"
    }
    """
    import uuid

    from flask_socketio import join_room, leave_room

    sid = request.sid
    r = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    room_id = data.get("roomId")
    client_type = data.get("clientType")  # "frontend" or "python"

    if not room_id:
        return {"status": "error", "code": 400, "message": "roomId required"}

    if not client_type or client_type not in ("frontend", "python"):
        return {
            "status": "error",
            "code": 400,
            "message": "clientType required ('frontend' or 'python')",
        }

    # Get userName from sid (must be authenticated via handle_connect)
    user_name = get_user_name_from_sid(sid)
    if not user_name:
        log.error(f"Cannot join room: userName not found for sid {sid}")
        return {"status": "error", "code": 401, "message": "User not found"}

    # Check room exists (creation is separate via POST /api/rooms)
    if not room_service.room_exists(room_id):
        log.debug(f"Room {room_id} not found for user {user_name}")
        return {"status": "error", "code": 404, "message": "Room not found"}

    # 1. Create sessionId
    session_id = str(uuid.uuid4())

    # 2. Store SID ↔ sessionId bidirectional mapping
    r.set(SessionKeys.session_to_sid(session_id), sid)
    session_keys = SessionKeys(sid)
    r.set(session_keys.session_id(), session_id)

    # 3. Store session data for route validation (lock routes, etc.)
    session_data = {
        "userId": user_name,
        "roomId": room_id,
        "createdAt": utc_now_iso(),
    }
    r.set(SessionKeys.session_data(session_id), json.dumps(session_data))

    # 4. Leave overview if joined, then join room-specific rooms
    leave_room("overview:public")
    join_room(f"room:{room_id}")
    join_room(f"session:{session_id}")

    # 5. Update user's current room, add to room membership, and track visit
    client_service = current_app.extensions["client_service"]
    client_service.update_user_and_room_membership(user_name, room_id)

    # 6. Register as frontend session and create camera (only for browser clients)
    room_keys = RoomKeys(room_id)
    if client_type == "frontend":
        r.sadd(room_keys.frontend_sessions(), session_id)
        create_session_camera(r, room_id, session_id)
        # Set initial active camera to session's own camera
        session_camera_key = get_session_camera_key(session_id)
        r.set(room_keys.session_active_camera(session_id), session_camera_key)

    # 7. Get minimal room state (geometries, settings, etc. fetched via REST)
    frame_count = r.zcard(room_keys.trajectory_indices())
    current_step = r.get(room_keys.current_frame())
    locked = r.exists(room_keys.locked())

    # 8. Send current progress trackers to joining client
    progress_data = r.hgetall(room_keys.progress())
    if progress_data:
        progress_trackers = {
            progress_id: json.loads(progress_json)
            for progress_id, progress_json in progress_data.items()
        }
        emit("progress:init", {"progressTrackers": progress_trackers}, to=sid)

    log.debug(
        f"User {user_name} joined room {room_id} with session {session_id} (sid: {sid})"
    )

    return {
        "status": "ok",
        "sessionId": session_id,
        "cameraKey": get_session_camera_key(session_id),
        "step": int(current_step) if current_step else 0,
        "frameCount": frame_count,
        "locked": bool(locked),
    }


@socketio.on("room:leave")
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
