"""Room management routes.

Handles room creation, listing, updates, metadata, locking, duplication, and schema management.
"""

import json
import logging
import re
import uuid

from flask import Blueprint, current_app, request

from zndraw.app.constants import SocketEvents
from zndraw.auth import get_current_user, require_admin, require_auth
from zndraw.extensions.analysis import analysis
from zndraw.extensions.modifiers import modifiers
from zndraw.extensions.selections import selections
from zndraw.server import socketio
from zndraw.utils.time import utc_now_iso

from .frame_index_manager import FrameIndexManager
from .redis_keys import GlobalIndexKeys, RoomKeys, SessionKeys
from .room_manager import emit_room_update
from .route_utils import (
    check_lock,
    emit_bookmarks_invalidate,
    get_lock_key,
    get_storage,
)

log = logging.getLogger(__name__)

rooms = Blueprint("rooms", __name__)


@rooms.route("/api/rooms", methods=["GET"])
def list_rooms():
    """List rooms with metadata based on user role.

    For admins, returns all rooms.
    For non-admins (guests and users), returns only rooms they've visited.

    Query Parameters:
        search: Optional regex pattern to search in metadata values

    Returns:
        [{
            "id": "room1",
            "description": "My room",
            "frameCount": 42,
            "locked": false,
            "metadataLocked": false,
            "isDefault": false,
            "metadata": {"relative_file_path": "...", ...}
        }]
    """
    from zndraw.app.room_data_fetcher import BatchedRoomDataFetcher
    from zndraw.auth import AuthError, decode_jwt_token, extract_token_from_request

    redis_client = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    admin_service = current_app.extensions["admin_service"]
    client_service = current_app.extensions["client_service"]
    search_pattern = request.args.get("search")

    # Get current user and role (if authenticated)
    is_admin = False
    user_name = None
    try:
        token = extract_token_from_request()
        if token:
            payload = decode_jwt_token(token)
            user_name = payload.get("sub")
            if user_name:
                is_admin = admin_service.is_admin(user_name)
    except AuthError:
        # Not authenticated - treat as non-admin
        pass

    # Get all room IDs
    all_room_ids = redis_client.smembers(GlobalIndexKeys.rooms_index())

    # Filter room IDs based on user role
    if is_admin:
        # Admins see all rooms
        room_ids = all_room_ids
    elif user_name:
        # Non-admin users see only rooms they've visited
        visited_rooms = client_service.get_visited_rooms(user_name)
        room_ids = all_room_ids & visited_rooms
    else:
        # Unauthenticated users see no rooms
        room_ids = set()

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
                "isDefault": default_room == room_id,
                "metadata": data["metadata"],
            }
        )

    return rooms, 200


@rooms.route("/api/rooms", methods=["POST"])
@require_auth
def create_room():
    """Create a new room explicitly.

    This is the preferred way to create rooms. Use this endpoint instead of
    relying on implicit room creation during join.

    Headers
    -------
    Authorization: Bearer <jwt-token> (required)

    Request
    -------
    {
        "roomId": "my-room",           // Required
        "description": "...",          // Optional
        "copyFrom": "source-room",     // Optional - copy data from existing room
        "template": "empty"            // Optional - apply template (empty, water, etc.)
    }

    Fallback Logic (handled by RoomService.create_room_with_defaults)
    ------------------------------------------------------------------
    1. If copyFrom specified → copy from that room
    2. Else if template specified → use that template (includes "none" for 0 frames)
    3. Else if admin set default_room → copy from default room
    4. Else → use "empty" template (1 empty frame)

    Response
    --------
    {
        "status": "ok",
        "roomId": "my-room",
        "frameCount": 0,
        "created": true
    }
    """
    data = request.get_json() or {}
    room_id = data.get("roomId")

    if not room_id:
        return {"error": "roomId is required"}, 400

    user_name = get_current_user()
    room_service = current_app.extensions["room_service"]

    if room_service.room_exists(room_id):
        return {"error": f"Room '{room_id}' already exists"}, 409

    description = data.get("description")
    copy_from = data.get("copyFrom")
    template = data.get("template")

    try:
        storage = get_storage(room_id)
        result = room_service.create_room_with_defaults(
            room_id=room_id,
            user_name=user_name,
            storage=storage,
            description=description,
            copy_from=copy_from,
            template=template,
        )
    except ValueError as e:
        return {"error": str(e)}, 400

    # Broadcast room creation
    emit_room_update(
        socketio,
        room_id,
        created=True,
        description=description,
        frameCount=result["frameCount"],
        locked=False,
        isDefault=False,
    )

    log.debug(f"User {user_name} created room {room_id} (source: {result['source']})")

    return {
        "status": "ok",
        "roomId": room_id,
        "frameCount": result["frameCount"],
        "created": True,
    }, 201


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
            "metadata": {"relative_file_path": "...", ...}
        }
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    redis_client = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    keys = RoomKeys(room_id)

    # Check if room exists using global index
    room_exists = redis_client.sismember(GlobalIndexKeys.rooms_index(), room_id)

    if not room_exists:
        return {"error": "Room not found"}, 404

    # Get frame count using service
    frame_count = room_service.get_frame_count(room_id)

    # Get metadata
    description = redis_client.get(keys.description())
    locked = redis_client.get(keys.locked()) == "1"

    # Get file metadata
    metadata_manager = RoomMetadataManager(redis_client, room_id)
    file_metadata = metadata_manager.get_all()

    return {
        "id": room_id,
        "description": description if description else None,
        "frameCount": frame_count,
        "locked": locked,
        "metadata": file_metadata,
    }, 200


@rooms.route("/api/rooms/<string:room_id>", methods=["PATCH"])
def update_room(room_id):
    """Update room metadata (description, locked).

    Request body:
        {
            "description": "My custom description",  // Optional
            "locked": true                           // Optional
        }

    Security:
    - In admin mode, only admins can lock/unlock rooms
    - Non-admins cannot unlock admin-locked rooms
    - Admin locks are tracked via locked_by field
    """
    from zndraw.auth import get_current_user

    redis_client = current_app.extensions["redis"]
    admin_service = current_app.extensions["admin_service"]
    config = current_app.extensions["config"]
    data = request.get_json() or {}
    keys = RoomKeys(room_id)

    # Get current user (if authenticated)
    try:
        current_user = get_current_user()
        is_admin = admin_service.is_admin(current_user)
    except Exception:
        # Not authenticated - treat as non-admin
        current_user = None
        is_admin = False

    # Check if room exists using global index
    room_exists = redis_client.sismember(GlobalIndexKeys.rooms_index(), room_id)

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

    # Update locked status with admin enforcement
    if "locked" in data:
        is_locking = bool(data["locked"])
        current_locked_by = redis_client.get(keys.locked_by())

        # In admin mode, only admins can lock/unlock
        if config.admin_username and not is_admin:
            return {"error": "Only admins can lock/unlock rooms in admin mode"}, 403

        # Prevent non-admins from unlocking admin-locked rooms
        if not is_locking and current_locked_by and not is_admin:
            return {
                "error": f"This room was locked by admin '{current_locked_by}' and cannot be unlocked by non-admins"
            }, 403

        # Update lock status
        redis_client.set(keys.locked(), "1" if is_locking else "0")
        changes["locked"] = is_locking

        # Track who locked it (only for admin locks)
        if is_locking and is_admin:
            redis_client.set(keys.locked_by(), current_user)
        elif not is_locking:
            # Remove locked_by when unlocking
            redis_client.delete(keys.locked_by())

    # Emit socket event for real-time updates
    if changes:
        from zndraw.app.room_manager import emit_room_update

        emit_room_update(socketio, room_id, **changes)
        log.debug(f"Emitted room:update for room '{room_id}': {changes}")

    log.debug(f"Updated room '{room_id}' metadata: {data}")
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
@require_admin
def set_default_room():
    """Set the default room (admin only).

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
        log.debug("Unset default room")

        # Update previous default room
        if previous_default:
            from zndraw.app.room_manager import emit_room_update

            emit_room_update(socketio, previous_default, isDefault=False)
    else:
        # Verify room exists using global index
        room_exists = redis_client.sismember(GlobalIndexKeys.rooms_index(), room_id)

        if not room_exists:
            return {"error": "Room not found"}, 404

        # Set default room
        redis_client.set("default_room", room_id)
        log.debug(f"Set default room to '{room_id}'")

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
        isDefault=False,
    )

    return {
        "status": "ok",
        "roomId": new_room_id,
        "frameCount": result["frameCount"],
    }, 200


@rooms.route("/api/rooms/<string:room_id>/renormalize", methods=["POST"])
@check_lock(target="trajectory:meta", forbid=["room:locked"])
def renormalize_frame_indices(room_id: str, session_id: str, user_id: str):
    """Renormalize frame indices to contiguous integers starting from 0.

    This is an O(N) operation that should be called infrequently when
    precision drift becomes an issue (e.g., scores too close together).

    Proceeds unless another session holds trajectory:meta lock or room is locked.

    Headers
    -------
    X-Session-ID : str
        Session ID from /join

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

    log.debug(f"Renormalized {count} frame indices for room '{room_id}'")

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
    # Map category strings to the corresponding imported objects
    # Note: "settings" category is now handled by dedicated /settings/ endpoints
    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "analysis": analysis,
    }

    if category not in category_map:
        return {"error": f"Unknown schema category '{category}'"}

    redis_client = current_app.extensions["redis"]
    # Changed from dict to list to support duplicate names with different public flags
    schema_list = []

    # Add server-provided extensions (Celery-based - always public)
    for name, cls in category_map[category].items():
        schema_list.append(
            {
                "name": name,
                "schema": cls.model_json_schema(),
                "provider": "celery",
                "public": True,  # Celery extensions are always public
                "queueLength": 0,
                "idleWorkers": 0,
                "progressingWorkers": 0,
            }
        )

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
            ext["name"] == name and ext["public"] is False for ext in schema_list
        )

        if existing:
            log.warning(
                f"{category.capitalize()} extension '{name}' (room-scoped) "
                "is already registered."
            )
        else:
            schema_list.append(
                {
                    "name": name,
                    "schema": sch,
                    "provider": stats.total_workers,  # Number of workers for client extensions
                    "public": False,  # Room-scoped extensions
                    **stats.to_dict(),
                }
            )

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
            ext["name"] == name and ext["public"] is True for ext in schema_list
        )

        if existing:
            log.warning(
                f"{category.capitalize()} global extension '{name}' "
                "is already registered."
            )
        else:
            schema_list.append(
                {
                    "name": name,
                    "schema": sch,
                    "provider": stats.total_workers,  # Number of workers for global extensions
                    "public": True,  # Global extensions
                    **stats.to_dict(),
                }
            )

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
            log.debug(
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
@check_lock(target="step")
def update_step(room_id: str, session_id: str, user_id: str):
    """Set current step/frame for a room.

    Proceeds unless another session holds the 'step' lock. Used for both:
    - Simple updates: PUT step (no lock needed unless coordinating)
    - Continuous updates with lock: acquire lock → PUT step × N → release lock

    Parameters
    ----------
    room_id : str
        Room identifier
    session_id : str
        Session ID (injected by @check_lock decorator)
    user_id : str
        User ID (injected by @check_lock decorator)

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
        Step lock held by another session

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
        progress_data = json.dumps(
            {"roomId": room_id, "description": description, "progress": None}
        )
        r.hset(room_keys.progress(), progress_id, progress_data)

        # Broadcast to room
        socketio.emit(
            "progress:start",
            {"progressId": progress_id, "roomId": room_id, "description": description},
            to=f"room:{room_id}",
        )

        return {"success": True, "progressId": progress_id}, 201
    except Exception as e:
        log.error(f"Failed to start progress tracking: {e}")
        return {"error": f"Failed to start progress tracking: {str(e)}"}, 500


@rooms.route(
    "/api/rooms/<string:room_id>/progress/<string:progress_id>", methods=["PUT"]
)
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

        socketio.emit(
            SocketEvents.PROGRESS_UPDATED, update_payload, to=f"room:{room_id}"
        )

        return {"success": True}, 200
    except Exception as e:
        log.error(f"Failed to update progress: {e}")
        return {"error": f"Failed to update progress: {str(e)}"}, 500


@rooms.route(
    "/api/rooms/<string:room_id>/progress/<string:progress_id>", methods=["DELETE"]
)
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
        socketio.emit(
            SocketEvents.PROGRESS_COMPLETED,
            {"progressId": progress_id, "roomId": room_id},
            to=f"room:{room_id}",
        )

        return {"success": True}, 200
    except Exception as e:
        log.error(f"Failed to complete progress: {e}")
        return {"error": f"Failed to complete progress: {str(e)}"}, 500
