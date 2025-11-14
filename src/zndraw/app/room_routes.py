"""Room management routes.

Handles room creation, listing, updates, metadata, locking, duplication, and schema management.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio
from zndraw.app.constants import SocketEvents
from zndraw.auth import require_auth

from .frame_index_manager import FrameIndexManager
from .redis_keys import RoomKeys, SessionKeys
from .room_manager import emit_room_update
from .route_utils import (
    emit_bookmarks_invalidate,
    get_lock_key,
    get_metadata_lock_info,
    requires_lock,
)

log = logging.getLogger(__name__)

rooms = Blueprint("rooms", __name__)


@rooms.route("/api/rooms", methods=["GET"])
def list_rooms():
    """List all active rooms with metadata.

    Query Parameters:
        search: Optional regex pattern to search in metadata values

    Returns:
        [{
            "id": "room1",
            "description": "My room",
            "frameCount": 42,
            "locked": false,
            "metadataLocked": false,
            "hidden": false,
            "isDefault": false,
            "metadata": {"relative_file_path": "...", ...}
        }]
    """
    import re

    from zndraw.app.room_data_fetcher import BatchedRoomDataFetcher

    redis_client = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    search_pattern = request.args.get("search")

    # Scan for all room keys to find unique room IDs
    room_ids = set()
    for key in redis_client.scan_iter(match="room:*"):
        # Extract room ID from keys like "room:{room_id}:..."
        parts = key.split(":")
        if len(parts) >= 2:
            room_ids.add(parts[1])

    # Get default room
    default_room = redis_client.get("default_room")

    # OPTIMIZATION: Batch fetch all room data to avoid N+1 queries
    fetcher = BatchedRoomDataFetcher()
    rooms_data = fetcher.fetch_rooms_data(redis_client, room_service, list(room_ids))

    # Build detailed room objects using pre-fetched data
    rooms = []
    for room_id in sorted(room_ids):
        data = rooms_data[room_id]

        # Filter by search pattern if provided
        if search_pattern:
            try:
                pattern = re.compile(search_pattern, re.IGNORECASE)
                # Search in metadata values and room ID
                search_targets = list(data["metadata"].values()) + [room_id]
                if data["description"]:
                    search_targets.append(data["description"])
                if not any(pattern.search(str(v)) for v in search_targets):
                    continue
            except re.error:
                # Invalid regex, skip filtering for this room
                pass

        rooms.append(
            {
                "id": room_id,
                "description": data["description"] if data["description"] else None,
                "frameCount": data["frameCount"],
                "locked": data["locked"],
                "metadataLocked": data["metadataLocked"],
                "hidden": data["hidden"],
                "isDefault": default_room == room_id,
                "metadata": data["metadata"],
            }
        )

    return rooms, 200


@rooms.route("/api/rooms/<string:room_id>", methods=["GET"])
@require_auth
def get_room(room_id):
    """Get details for a specific room.

    Returns:
        {
            "id": "room1",
            "description": "My room",
            "frameCount": 42,
            "locked": false,
            "hidden": false,
            "metadata": {"relative_file_path": "...", ...}
        }
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    redis_client = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    keys = RoomKeys(room_id)

    # Check if room exists
    room_exists = False
    for key in redis_client.scan_iter(match=keys.all_keys_pattern(), count=1):
        room_exists = True
        break

    if not room_exists:
        return {"error": "Room not found"}, 404

    # Get frame count using service
    frame_count = room_service.get_frame_count(room_id)

    # Get metadata
    description = redis_client.get(keys.description())
    locked = redis_client.get(keys.locked()) == "1"
    hidden = redis_client.get(keys.hidden()) == "1"

    # Get file metadata
    metadata_manager = RoomMetadataManager(redis_client, room_id)
    file_metadata = metadata_manager.get_all()

    return {
        "id": room_id,
        "description": description if description else None,
        "frameCount": frame_count,
        "locked": locked,
        "hidden": hidden,
        "metadata": file_metadata,
    }, 200


@rooms.route("/api/rooms/<string:room_id>/join", methods=["POST"])
def join_room(room_id):
    """Join a room (requires JWT authentication).

    Generates a unique sessionId for this client connection.
    The sessionId identifies this specific browser tab/client instance
    and is used for per-connection lock ownership.

    Headers
    -------
    Authorization: Bearer <jwt-token> (required)

    Request
    -------
    {
        "description": "optional room description",
        "copyFrom": "optional-source-room",
        "allowCreate": true
    }

    Response
    --------
    {
        "status": "ok",
        "roomId": "room-name",
        "userName": "username",
        "sessionId": "unique-session-uuid",
        "created": true
    }

    Note
    ----
    Room data (frameCount, selections, geometries, bookmarks, settings, etc.)
    should be fetched via separate REST endpoints after joining.
    """
    import datetime
    import uuid

    from zndraw.auth import AuthError, get_current_user

    data = request.get_json() or {}
    if ":" in room_id:
        return {"error": "Room ID cannot contain ':' character"}, 400

    # Authenticate request (JWT required)
    try:
        user_name = get_current_user()
    except AuthError as e:
        return {"error": e.message}, e.status_code

    # Generate unique session ID for this client connection
    session_id = str(uuid.uuid4())
    description = data.get("description")
    copy_from = data.get("copyFrom")
    allow_create = data.get("allowCreate", True)
    r = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    client_service = current_app.extensions["client_service"]
    settings_service = current_app.extensions["settings_service"]

    # Store session mapping in Redis (sessionId → userId)
    # TTL: 24 hours (session expires if no activity)
    session_key = SessionKeys.session_data(session_id)
    session_data = json.dumps({
        "userId": user_name,
        "roomId": room_id,
        "connectedAt": datetime.datetime.utcnow().isoformat(),
        "lastActivity": datetime.datetime.utcnow().isoformat()
    })
    r.set(session_key, session_data, ex=86400)  # 24 hour TTL

    log.info(f"User {user_name} joined room {room_id} with session {session_id}")

    # Check if room already exists
    room_exists = room_service.room_exists(room_id)

    # If allowCreate is False and room doesn't exist, return 404
    if not allow_create and not room_exists:
        return {
            "status": "not_found",
            "message": f"Room '{room_id}' does not exist yet. It may still be loading.",
        }, 404

    # Update client room membership atomically (using userName as identifier)
    client_service.update_user_and_room_membership(user_name, room_id)

    response = {
        "status": "ok",
        "userName": user_name,
        "sessionId": session_id,  # Return session ID to client
        "roomId": room_id,
        "created": not room_exists,
    }

    if not room_exists:
        # Create new room using RoomService
        try:
            result = room_service.create_room(
                room_id, user_name, description, copy_from
            )
            frame_count = result["frameCount"]
        except ValueError as e:
            return {"error": str(e)}, 404

        # Include frameCount in response when creating room (especially for copyFrom)
        response["frameCount"] = frame_count

        # Initialize default settings for new user in room
        settings_service.initialize_defaults(room_id, user_name)

        # Broadcast room creation to all connected clients
        emit_room_update(
            socketio,
            room_id,
            created=True,
            description=description,
            frameCount=frame_count,
            locked=False,
            hidden=False,
            isDefault=False,
        )

    return response


@rooms.route("/api/rooms/<string:room_id>", methods=["PATCH"])
def update_room(room_id):
    """Update room metadata (description, locked, hidden).

    Request body:
        {
            "description": "My custom description",  // Optional
            "locked": true,                          // Optional
            "hidden": false                          // Optional
        }
    """
    redis_client = current_app.extensions["redis"]
    data = request.get_json() or {}
    keys = RoomKeys(room_id)

    # Check if room exists
    room_exists = False
    for key in redis_client.scan_iter(match=keys.all_keys_pattern(), count=1):
        room_exists = True
        break

    if not room_exists:
        return {"error": "Room not found"}, 404

    # Track what changed for socket event
    changes = {}

    # Update description
    if "description" in data:
        if data["description"] is None:
            redis_client.delete(keys.description())
            changes["description"] = None
        else:
            redis_client.set(keys.description(), data["description"])
            changes["description"] = data["description"]

    # Update locked status
    if "locked" in data:
        redis_client.set(keys.locked(), "1" if data["locked"] else "0")
        changes["locked"] = bool(data["locked"])

    # Update hidden status
    if "hidden" in data:
        redis_client.set(keys.hidden(), "1" if data["hidden"] else "0")
        changes["hidden"] = bool(data["hidden"])

    # Emit socket event for real-time updates
    if changes:
        from zndraw.app.room_manager import emit_room_update

        emit_room_update(socketio, room_id, **changes)
        log.debug(f"Emitted room:update for room '{room_id}': {changes}")

    log.info(f"Updated room '{room_id}' metadata: {data}")
    return {"status": "ok"}, 200


@rooms.route("/api/rooms/default", methods=["GET"])
def get_default_room():
    """Get the default room ID.

    Returns:
        {"roomId": "room1"} or {"roomId": null}
    """
    redis_client = current_app.extensions["redis"]
    default_room = redis_client.get("default_room")
    return {"roomId": default_room if default_room else None}, 200


@rooms.route("/api/rooms/default", methods=["PUT"])
def set_default_room():
    """Set the default room.

    Request body:
        {"roomId": "room1"}  // or null to unset
    """
    redis_client = current_app.extensions["redis"]
    data = request.get_json() or {}

    room_id = data.get("roomId")

    # Get previous default room
    previous_default = redis_client.get("default_room")

    if room_id is None:
        # Unset default room
        redis_client.delete("default_room")
        log.info("Unset default room")

        # Update previous default room
        if previous_default:
            from zndraw.app.room_manager import emit_room_update

            emit_room_update(socketio, previous_default, isDefault=False)
    else:
        # Verify room exists
        room_keys = RoomKeys(room_id)
        room_exists = False
        for key in redis_client.scan_iter(match=room_keys.all_keys_pattern(), count=1):
            room_exists = True
            break

        if not room_exists:
            return {"error": "Room not found"}, 404

        # Set default room
        redis_client.set("default_room", room_id)
        log.info(f"Set default room to '{room_id}'")

        # Update previous default room if different
        if previous_default and previous_default != room_id:
            from zndraw.app.room_manager import emit_room_update

            emit_room_update(socketio, previous_default, isDefault=False)

        # Update new default room
        from zndraw.app.room_manager import emit_room_update

        emit_room_update(socketio, room_id, isDefault=True)

    log.debug(f"Updated default room from {previous_default} to {room_id}")

    return {"status": "ok"}, 200


@rooms.route("/api/rooms/<string:room_id>/metadata", methods=["GET"])
def get_room_metadata(room_id: str):
    """Get all metadata for a room.

    Returns:
        {"metadata": {"relative_file_path": "...", ...}}
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    redis_client = current_app.extensions["redis"]
    manager = RoomMetadataManager(redis_client, room_id)
    metadata = manager.get_all()

    return {"metadata": metadata}, 200


@rooms.route("/api/rooms/<string:room_id>/metadata", methods=["POST"])
def update_room_metadata(room_id: str):
    """Update room metadata. Respects room lock.

    Request body:
        {"relative_file_path": "data/file.xyz", "file_size": "12345", ...}

    Returns:
        {"success": true, "metadata": {...}}
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    # Check permanent room lock only (metadata doesn't use trajectory lock)
    redis_client = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    locked = redis_client.get(keys.locked())
    if locked == "1":
        return {"error": "Room is locked and cannot be modified"}, 403

    data = request.get_json() or {}

    # Validate all values are strings
    for key, value in data.items():
        if not isinstance(value, str):
            return {
                "error": f"All values must be strings. Field '{key}' has type {type(value).__name__}"
            }, 400

    # Update metadata
    manager = RoomMetadataManager(redis_client, room_id)
    try:
        manager.update(data)
        metadata = manager.get_all()
        return {"success": True, "metadata": metadata}, 200
    except Exception as e:
        log.error(f"Error updating metadata for room '{room_id}': {e}")
        return {"error": str(e)}, 500


@rooms.route("/api/rooms/<string:room_id>/metadata/<string:field>", methods=["DELETE"])
def delete_room_metadata_field(room_id: str, field: str):
    """Delete specific metadata field. Respects room lock.

    Returns:
        {"success": true}
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    # Check permanent room lock only (metadata doesn't use trajectory lock)
    redis_client = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    locked = redis_client.get(keys.locked())
    if locked == "1":
        return {"error": "Room is locked and cannot be modified"}, 403

    manager = RoomMetadataManager(redis_client, room_id)

    try:
        deleted = manager.delete(field)
        return {"success": True, "deleted": deleted}, 200
    except Exception as e:
        log.error(f"Error deleting metadata field '{field}' for room '{room_id}': {e}")
        return {"error": str(e)}, 500


@rooms.route("/api/rooms/<string:room_id>/locks/<string:target>", methods=["GET"])
def get_lock_status(room_id: str, target: str):
    """Get current lock status and metadata for a specific target.

    Returns who holds the lock and what they're doing with it.
    Useful for frontend to display lock status to users.

    Parameters
    ----------
    room_id : str
        Room identifier
    target : str
        Lock target (e.g., 'trajectory:meta')

    Returns
    -------
    dict
        Lock status including holder, metadata, and TTL

    Example response (unlocked):
        {"locked": false, "target": "trajectory:meta"}

    Example response (locked):
        {
            "locked": true,
            "target": "trajectory:meta",
            "holder": "alice",
            "metadata": {
                "userName": "alice",
                "timestamp": 1234567890.123,
                "msg": "Uploading trajectory data"
            },
            "ttl": 45
        }
    """
    redis_client = current_app.extensions["redis"]
    lock_key = get_lock_key(room_id, target)

    lock_holder = redis_client.get(lock_key)
    if not lock_holder:
        return {"locked": False, "target": target}

    # Get metadata
    metadata_key = f"{lock_key}:metadata"
    metadata_json = redis_client.get(metadata_key)

    metadata = {}
    if metadata_json:
        try:
            metadata = json.loads(metadata_json)
        except json.JSONDecodeError:
            log.error(f"Failed to parse lock metadata: {metadata_json}")

    # Get TTL
    ttl = redis_client.ttl(lock_key)

    return {
        "locked": True,
        "target": target,
        "holder": lock_holder,
        "metadata": metadata,
        "ttl": ttl,
    }


@rooms.route("/api/rooms/<string:room_id>/duplicate", methods=["POST"])
def duplicate_room(room_id):
    """Duplicate a room by copying all frame mappings and metadata.

    Request body:
        {
            "newRoomId": "new-room-uuid",  // Optional, auto-generated if not provided
            "description": "Copy of room1"  // Optional, description for new room
        }

    Returns:
        {
            "status": "ok",
            "roomId": "new-room-uuid",
            "frameCount": 42
        }
    """
    import uuid

    room_service = current_app.extensions["room_service"]
    data = request.get_json() or {}

    # Check source room exists
    if not room_service.room_exists(room_id):
        return {"error": "Source room not found"}, 404

    # Generate or use provided new room ID
    new_room_id = data.get("newRoomId")
    if not new_room_id:
        new_room_id = str(uuid.uuid4())

    # Check new room doesn't already exist
    try:
        room_service.validate_room_available(new_room_id)
    except ValueError:
        return {"error": "Room with that ID already exists"}, 409

    # Create room by copying from source
    description = data.get("description")
    try:
        result = room_service.create_room(
            room_id=new_room_id,
            user_name="",  # Not needed for copy operations
            description=description,
            copy_from=room_id,
        )
    except ValueError as e:
        return {"error": str(e)}, 400

    # Emit room:update event to notify clients of new room
    from zndraw.app.room_manager import emit_room_update

    emit_room_update(
        socketio,
        new_room_id,
        created=True,
        description=description,
        frameCount=result["frameCount"],
        locked=False,
        hidden=False,
        isDefault=False,
    )

    return {
        "status": "ok",
        "roomId": new_room_id,
        "frameCount": result["frameCount"],
    }, 200


@rooms.route("/api/rooms/<string:room_id>/renormalize", methods=["POST"])
def renormalize_frame_indices(room_id):
    """Renormalize frame indices to contiguous integers starting from 0.

    This is an O(N) operation that should be called infrequently when
    precision drift becomes an issue (e.g., scores too close together).

    Returns
    -------
    dict
        JSON response with status and number of frames renormalized
    """
    redis_client = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    indices_key = keys.trajectory_indices()

    # Check if room exists
    room_exists = redis_client.exists(indices_key)
    if not room_exists:
        return {"error": "Room not found"}, 404

    # Renormalize the indices
    manager = FrameIndexManager(redis_client, indices_key)
    count = manager.renormalize()

    log.info(f"Renormalized {count} frame indices for room '{room_id}'")

    # Emit bookmarks update since logical positions haven't changed
    # but we should still notify clients
    emit_bookmarks_invalidate(room_id)

    return {"status": "ok", "framesRenormalized": count}, 200


@rooms.route("/api/rooms/<string:room_id>/schema/<string:category>", methods=["GET"])
def get_room_schema(room_id: str, category: str):
    """Get the schema for a specific room with worker and queue statistics.

    Returns schema along with metadata about each extension:
    - provider: "celery" for server-side extensions, or count of registered workers
    - queueLength: number of queued tasks for this extension
    - idleWorkers: number of idle workers available
    - progressingWorkers: number of workers currently processing tasks
    """
    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    # Map category strings to the corresponding imported objects
    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
        "analysis": analysis,
    }

    if category not in category_map:
        return {"error": f"Unknown schema category '{category}'"}

    redis_client = current_app.extensions["redis"]
    # Changed from dict to list to support duplicate names with different public flags
    schema_list = []

    # Add server-provided extensions (Celery-based - always public)
    for name, cls in category_map[category].items():
        schema_list.append({
            "name": name,
            "schema": cls.model_json_schema(),
            "provider": "celery",
            "public": True,  # Celery extensions are always public
            "queueLength": 0,
            "idleWorkers": 0,
            "progressingWorkers": 0,
        })

    # Add client-provided extensions from Redis (room-scoped)
    from .redis_keys import ExtensionKeys
    from .worker_stats import WorkerStats

    schema_key = ExtensionKeys.schema_key(room_id, category)
    redis_schema = redis_client.hgetall(schema_key)

    for name, sch_str in redis_schema.items():
        sch = json.loads(sch_str)

        # Get worker statistics for this extension
        keys = ExtensionKeys.for_extension(room_id, category, name)
        stats = WorkerStats.fetch(redis_client, keys)

        # Check if this exact combination (name + public=False) already exists
        existing = any(
            ext["name"] == name and ext["public"] is False
            for ext in schema_list
        )

        if existing:
            log.warning(
                f"{category.capitalize()} extension '{name}' (room-scoped) "
                "is already registered."
            )
        else:
            schema_list.append({
                "name": name,
                "schema": sch,
                "provider": stats.total_workers,  # Number of workers for client extensions
                "public": False,  # Room-scoped extensions
                **stats.to_dict(),
            })

    # Add global (public) extensions from Redis
    global_schema_key = ExtensionKeys.global_schema_key(category)
    global_redis_schema = redis_client.hgetall(global_schema_key)

    for name, sch_str in global_redis_schema.items():
        sch = json.loads(sch_str)

        # Get worker statistics for this global extension
        keys = ExtensionKeys.for_global_extension(category, name)
        stats = WorkerStats.fetch(redis_client, keys)

        # Check if this exact combination (name + public=True) already exists
        existing = any(
            ext["name"] == name and ext["public"] is True
            for ext in schema_list
        )

        if existing:
            log.warning(
                f"{category.capitalize()} global extension '{name}' "
                "is already registered."
            )
        else:
            schema_list.append({
                "name": name,
                "schema": sch,
                "provider": stats.total_workers,  # Number of workers for global extensions
                "public": True,  # Global extensions
                **stats.to_dict(),
            })

    return schema_list


@rooms.route("/api/rooms/<string:room_id>/step", methods=["GET"])
@require_auth
def get_step(room_id: str):
    """Get current step/frame for a room.

    Parameters
    ----------
    room_id : str
        Room identifier

    Returns
    -------
    dict
        Current step and total frame count

    Example response:
        {
            "step": 42,
            "totalFrames": 100
        }
    """
    redis_client = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]

    # Get current frame and total count
    current_frame = room_service.get_current_frame(room_id)
    total_frames = room_service.get_frame_count(room_id)

    # Clamp current_frame to valid range to prevent 404s when frames are deleted
    # or not yet uploaded
    if total_frames > 0 and current_frame is not None:
        clamped_frame = max(0, min(current_frame, total_frames - 1))
        if clamped_frame != current_frame:
            # Update Redis to persist the clamped value
            log.info(
                f"Clamping frame {current_frame} to {clamped_frame} for room {room_id} "
                f"(total frames: {total_frames})"
            )
            keys = RoomKeys(room_id)
            redis_client.set(keys.current_frame(), clamped_frame)
            current_frame = clamped_frame

    return {
        "step": current_frame,
        "totalFrames": total_frames,
    }


@rooms.route("/api/rooms/<string:room_id>/step", methods=["PUT"])
@requires_lock(target="step")
def update_step(room_id: str, session_id: str, user_id: str):
    """Set current step/frame for a room.

    Requires holding the 'step' lock for the room. Used for both:
    - Atomic updates: acquire lock → PUT step → release lock
    - Continuous updates: acquire lock → PUT step × N → release lock

    Parameters
    ----------
    room_id : str
        Room identifier
    session_id : str
        Session ID (injected by @requires_lock decorator)
    user_id : str
        User ID (injected by @requires_lock decorator)

    Request Body
    ------------
    {
        "step": 42  // Frame index (non-negative integer)
    }

    Returns
    -------
    dict
        Success status and updated step

    Raises
    ------
    400
        Invalid step value (negative or non-integer)
    401
        Invalid JWT or session
    403
        Session/user mismatch
    423
        Lock not held for target="step"

    Example response:
        {
            "success": true,
            "step": 42
        }

    Side Effects
    ------------
    Emits 'frame_update' Socket.IO event to all clients in the room
    (excluding the session that made the update).
    """
    redis_client = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    data = request.get_json()
    step = data.get("step")

    # Validate step
    # Handle both integer and float strings (Plotly may send "42.5")
    try:
        step = int(float(step))
        if step < 0:
            return {"error": "Step must be non-negative"}, 400
    except (ValueError, TypeError):
        return {"error": "Invalid step value"}, 400

    # Validate upper bound (only if room has frames)
    frame_count = room_service.get_frame_count(room_id)
    if frame_count > 0 and step >= frame_count:
        return {
            "error": f"Step {step} out of range. Room has {frame_count} frames (valid range: 0-{frame_count - 1})"
        }, 400

    # Update Redis
    keys = RoomKeys(room_id)
    redis_client.set(keys.current_frame(), step)

    # Get socket ID for this session to skip emitting to the sender
    sid = redis_client.get(SessionKeys.session_to_sid(session_id))
    if sid:
        sid = sid.decode() if isinstance(sid, bytes) else sid

    # Emit frame update for live updates (skip the session that made the change)
    socketio.emit(
        SocketEvents.FRAME_UPDATE,
        {"frame": step},
        to=f"room:{room_id}",
        skip_sid=sid,
    )
    return {"success": True, "step": step}


# ==================== Progress Tracking REST API ====================


@rooms.route("/api/rooms/<string:room_id>/progress", methods=["POST"])
@require_auth
def create_progress(room_id):
    """Start tracking progress for a long-running operation.

    Request Body
    ------------
    {
        "progressId": "unique-progress-id",
        "description": "Operation description"
    }

    Returns
    -------
    Response
        JSON response with success status
    """
    data = request.get_json()
    if not data:
        return {"error": "Request body required"}, 400

    progress_id = data.get("progressId")
    description = data.get("description")

    if not progress_id or not description:
        return {"error": "Missing required fields: progressId, description"}, 400

    r = current_app.extensions["redis"]

    try:
        # Store progress in Redis
        room_keys = RoomKeys(room_id)
        progress_data = json.dumps({"roomId": room_id, "description": description, "progress": None})
        r.hset(room_keys.progress(), progress_id, progress_data)

        # Broadcast to room
        socketio.emit(
            "progress:started",
            {"progressId": progress_id, "roomId": room_id, "description": description},
            to=f"room:{room_id}",
        )

        return {"success": True, "progressId": progress_id}, 201
    except Exception as e:
        log.error(f"Failed to start progress tracking: {e}")
        return {"error": f"Failed to start progress tracking: {str(e)}"}, 500


@rooms.route("/api/rooms/<string:room_id>/progress/<string:progress_id>", methods=["PUT"])
@require_auth
def update_progress(room_id, progress_id):
    """Update progress tracking for an ongoing operation.

    Request Body
    ------------
    {
        "description": "Updated description (optional)",
        "progress": 50.0  (optional, 0-100)
    }

    Returns
    -------
    Response
        JSON response with success status
    """
    data = request.get_json()
    if not data:
        return {"error": "Request body required"}, 400

    description = data.get("description")
    progress = data.get("progress")

    if description is None and progress is None:
        return {"error": "At least one of description or progress is required"}, 400

    r = current_app.extensions["redis"]

    try:
        # Get existing progress entry
        room_keys = RoomKeys(room_id)
        existing_data = r.hget(room_keys.progress(), progress_id)

        if not existing_data:
            return {"error": "Progress tracker not found"}, 404

        # Update progress data
        progress_dict = json.loads(existing_data)
        if description is not None:
            progress_dict["description"] = description
        if progress is not None:
            progress_dict["progress"] = progress

        r.hset(room_keys.progress(), progress_id, json.dumps(progress_dict))

        # Broadcast update to room
        update_payload = {"progressId": progress_id, "roomId": room_id}
        if description is not None:
            update_payload["description"] = description
        if progress is not None:
            update_payload["progress"] = progress

        socketio.emit(SocketEvents.PROGRESS_UPDATED, update_payload, to=f"room:{room_id}")

        return {"success": True}, 200
    except Exception as e:
        log.error(f"Failed to update progress: {e}")
        return {"error": f"Failed to update progress: {str(e)}"}, 500


@rooms.route("/api/rooms/<string:room_id>/progress/<string:progress_id>", methods=["DELETE"])
@require_auth
def complete_progress(room_id, progress_id):
    """Complete and remove progress tracking for an operation.

    Returns
    -------
    Response
        JSON response with success status
    """
    r = current_app.extensions["redis"]

    try:
        # Remove progress from Redis
        room_keys = RoomKeys(room_id)
        removed_count = r.hdel(room_keys.progress(), progress_id)

        if removed_count == 0:
            return {"error": "Progress tracker not found"}, 404

        # Broadcast completion to room
        socketio.emit(SocketEvents.PROGRESS_COMPLETED, {"progressId": progress_id, "roomId": room_id}, to=f"room:{room_id}")

        return {"success": True}, 200
    except Exception as e:
        log.error(f"Failed to complete progress: {e}")
        return {"error": f"Failed to complete progress: {str(e)}"}, 500
