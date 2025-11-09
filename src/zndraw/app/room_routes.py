"""Room management routes.

Handles room creation, listing, updates, metadata, locking, duplication, and schema management.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio

from .frame_index_manager import FrameIndexManager
from .redis_keys import RoomKeys
from .room_manager import emit_room_update
from .route_utils import (
    emit_bookmarks_invalidate,
    get_lock_key,
    get_metadata_lock_info,
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

    from zndraw.app.metadata_manager import RoomMetadataManager

    redis_client = current_app.extensions["redis"]
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

    # Build detailed room objects
    room_service = current_app.extensions["room_service"]
    rooms = []
    for room_id in sorted(room_ids):
        keys = RoomKeys(room_id)

        # Get frame count using service
        frame_count = room_service.get_frame_count(room_id)

        # Get metadata
        description = redis_client.get(keys.description())
        locked = redis_client.get(keys.locked()) == "1"
        hidden = redis_client.get(keys.hidden()) == "1"
        is_default = default_room == room_id

        # Check if metadata lock is held (trajectory:meta is the target used by vis.lock)
        metadata_lock_key = get_lock_key(room_id, "trajectory:meta")
        metadata_locked = redis_client.get(metadata_lock_key) is not None

        # Get file metadata
        metadata_manager = RoomMetadataManager(redis_client, room_id)
        file_metadata = metadata_manager.get_all()

        # Filter by search pattern if provided
        if search_pattern:
            try:
                pattern = re.compile(search_pattern, re.IGNORECASE)
                # Search in metadata values and room ID
                search_targets = list(file_metadata.values()) + [room_id]
                if description:
                    search_targets.append(description)
                if not any(pattern.search(str(v)) for v in search_targets):
                    continue
            except re.error:
                # Invalid regex, skip filtering for this room
                pass

        rooms.append(
            {
                "id": room_id,
                "description": description if description else None,
                "frameCount": frame_count,
                "locked": locked,
                "metadataLocked": metadata_locked,
                "hidden": hidden,
                "isDefault": is_default,
                "metadata": file_metadata,
            }
        )

    return rooms, 200


@rooms.route("/api/rooms/<string:room_id>", methods=["GET"])
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
    for key in redis_client.scan_iter(match=f"room:{room_id}:*", count=1):
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
        "frameCount": 0,
        "step": 0,
        "created": true,
        ...
    }
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
    session_key = f"session:{session_id}"
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
        "frameCount": 0,
        "roomId": room_id,
        "selections": None,
        "frame_selection": None,
        "created": True,
        "presenter-lock": False,
        "step": None,
        "geometries": None,
    }

    response["created"] = not room_exists

    if not room_exists:
        # Create new room using RoomService
        try:
            result = room_service.create_room(
                room_id, user_name, description, copy_from
            )
            frame_count = result["frameCount"]
        except ValueError as e:
            return {"error": str(e)}, 404

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

    # Get frame count using service
    response["frameCount"] = room_service.get_frame_count(room_id)

    keys = RoomKeys(room_id)

    selections_raw = r.hgetall(keys.selections())
    selections = {k: json.loads(v) for k, v in selections_raw.items()}
    response["selections"] = selections

    frame_selection = r.get(keys.frame_selection())
    response["frame_selection"] = (
        json.loads(frame_selection) if frame_selection else None
    )

    presenter_lock = r.get(keys.presenter_lock())
    response["presenter-lock"] = presenter_lock

    # Get current frame using service (handles validation and error cases)
    response["step"] = room_service.get_current_frame(room_id)

    bookmarks_raw = r.hgetall(keys.bookmarks())

    geometries = r.hgetall(keys.geometries())
    response["geometries"] = {k: json.loads(v) for k, v in geometries.items()}

    # Add geometry defaults from Pydantic models (single source of truth)
    from zndraw.geometries import geometries as geometry_models

    geometry_defaults = {}
    for name, model in geometry_models.items():
        try:
            # Try to instantiate with no arguments to get defaults
            geometry_defaults[name] = model().model_dump()
        except Exception:
            # Skip geometries that require arguments (e.g., Camera requires curve references)
            # These geometries don't have meaningful defaults without context
            pass
    response["geometryDefaults"] = geometry_defaults

    # Convert bookmark keys from strings (Redis) to integers
    if bookmarks_raw:
        response["bookmarks"] = {int(k): v for k, v in bookmarks_raw.items()}
    else:
        response["bookmarks"] = None

    # Fetch all settings for the user using service
    response["settings"] = settings_service.get_all(room_id, user_name)

    # Fetch selection groups
    groups_raw = r.hgetall(keys.selection_groups())
    selection_groups = {}
    for group_name, group_data in groups_raw.items():
        selection_groups[group_name] = json.loads(group_data)
    response["selectionGroups"] = selection_groups

    # Fetch active selection group
    active_group = r.get(keys.active_selection_group())
    response["activeSelectionGroup"] = active_group if active_group else None

    # Check if metadata lock is held
    response["metadataLocked"] = get_metadata_lock_info(room_id)

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
    for key in redis_client.scan_iter(match=f"room:{room_id}:*", count=1):
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
        room_exists = False
        for key in redis_client.scan_iter(match=f"room:{room_id}:*", count=1):
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
    schema = {}

    # Add server-provided extensions (Celery-based)
    for name, cls in category_map[category].items():
        schema[name] = {
            "schema": cls.model_json_schema(),
            "provider": "celery",
            "queueLength": 0,
            "idleWorkers": 0,
            "progressingWorkers": 0,
        }

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

        if name in schema:
            if schema[name]["schema"] != sch:
                log.warning(
                    f"{category.capitalize()} extension '{name}' schema "
                    "in Redis differs from server schema."
                )
        else:
            schema[name] = {
                "schema": sch,
                "provider": stats.total_workers,  # Number of workers for client extensions
                **stats.to_dict(),
            }

    # Add global (public) extensions from Redis
    global_schema_key = ExtensionKeys.global_schema_key(category)
    global_redis_schema = redis_client.hgetall(global_schema_key)

    for name, sch_str in global_redis_schema.items():
        sch = json.loads(sch_str)

        # Get worker statistics for this global extension
        keys = ExtensionKeys.for_global_extension(category, name)
        stats = WorkerStats.fetch(redis_client, keys)

        if name in schema:
            if schema[name]["schema"] != sch:
                log.warning(
                    f"{category.capitalize()} global extension '{name}' schema "
                    "differs from existing schema (room or server)."
                )
        else:
            schema[name] = {
                "schema": sch,
                "provider": stats.total_workers,  # Number of workers for global extensions
                **stats.to_dict(),
            }

    return schema
