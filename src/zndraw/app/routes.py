import json
import logging
import traceback
from pathlib import Path

import msgpack
import zarr
from flask import Response, current_app, request, send_file, send_from_directory
from flask_socketio import disconnect
from zarr.storage import MemoryStore

from zndraw.server import socketio
from zndraw.storage import ZarrStorageSequence, decode_data, encode_data

from . import main
from .constants import SocketEvents
from .frame_index_manager import FrameIndexManager
from .job_manager import JobManager
from .models import LockMetadata
from .queue_manager import emit_queue_update
from .redis_keys import ExtensionKeys
from .worker_stats import WorkerStats

# --- Logging Setup ---
log = logging.getLogger(__name__)

STORAGE: dict[str, MemoryStore] = {}


# TODO: move to utils
def get_lock_key(room: str, target: str) -> str:
    """Constructs a standardized Redis key for a lock."""
    return f"room:{room}:lock:{target}"


def get_zarr_store_path(room_id: str) -> str:
    """Returns the path to the Zarr store for a given room."""
    storage_path = current_app.config.get("STORAGE_PATH", "./zndraw-data.zarr")
    # Remove .zarr extension if present to append room_id
    base_path = storage_path.rstrip("/").removesuffix(".zarr")
    return f"{base_path}/{room_id}.zarr"


def get_storage(room_id: str) -> ZarrStorageSequence:
    # store_path = get_zarr_store_path(room_id)
    if room_id not in STORAGE:
        STORAGE[room_id] = MemoryStore()
    root = zarr.group(STORAGE[room_id])
    storage = ZarrStorageSequence(root)
    return storage


def check_room_locked(
    room_id: str, client_id: str | None = None
) -> tuple[dict[str, str], int] | None:
    """
    Check if a room is locked. Returns error response tuple if locked, None if not locked.

    This checks BOTH:
    1. Permanent room lock (room:{room_id}:locked)
    2. Temporary trajectory lock (trajectory:meta) - used by vis.lock context manager

    Args:
        room_id: The room ID to check
        client_id: Optional client ID - if provided and this client holds the trajectory lock,
                   the check will pass

    Returns:
        Error tuple if locked and client doesn't hold lock, None otherwise
    """
    redis_client = current_app.extensions["redis"]

    # Check permanent room lock
    locked = redis_client.get(f"room:{room_id}:locked")
    if locked == "1":
        return {"error": "Room is locked and cannot be modified"}, 403

    # Check trajectory lock (from vis.lock context manager)
    # This protects ALL trajectory operations: append, delete, upload, etc.
    trajectory_lock_key = get_lock_key(room_id, "trajectory:meta")
    lock_holder = redis_client.get(trajectory_lock_key)

    if lock_holder:
        log.debug(
            f"Lock check: room={room_id}, lock_holder={lock_holder}, client_id={client_id}"
        )
        # If a client_id is provided and it holds the lock, allow the operation
        if client_id and lock_holder == client_id:
            return None
        # Otherwise, the room is locked by another operation
        log.warning(
            f"Lock rejected: lock_holder={lock_holder} != client_id={client_id}"
        )
        return {"error": "Room is temporarily locked by another operation"}, 423

    return None


def shift_bookmarks_on_delete(room_id: str, deleted_indices: list[int]):
    """
    Shift bookmark indices after frames are deleted.
    Bookmarks at deleted indices are removed, bookmarks after are shifted down.

    Parameters
    ----------
    room_id : str
        The room ID
    deleted_indices : list[int]
        Sorted list of deleted frame indices
    """
    if not deleted_indices:
        return

    r = current_app.extensions["redis"]
    bookmarks_key = f"room:{room_id}:bookmarks"
    bookmarks_raw = r.hgetall(bookmarks_key)

    if not bookmarks_raw:
        return

    # Convert to int keys and sort deletion indices
    bookmarks = {int(k): v for k, v in bookmarks_raw.items()}
    deleted_set = set(deleted_indices)

    # Build new bookmarks dict with shifted indices
    new_bookmarks = {}
    for idx, label in bookmarks.items():
        # Remove bookmarks at deleted indices
        if idx in deleted_set:
            continue

        # Calculate how many indices below this one were deleted
        shift = sum(1 for del_idx in deleted_indices if del_idx < idx)
        new_idx = idx - shift
        new_bookmarks[new_idx] = label

    # Update Redis with shifted bookmarks
    r.delete(bookmarks_key)
    if new_bookmarks:
        r.hset(bookmarks_key, mapping={str(k): v for k, v in new_bookmarks.items()})


def shift_bookmarks_on_insert(room_id: str, insert_position: int):
    """
    Shift bookmark indices after a frame is inserted.
    All bookmarks at or after the insert position shift up by 1.

    Parameters
    ----------
    room_id : str
        The room ID
    insert_position : int
        The index where the frame was inserted
    """
    r = current_app.extensions["redis"]
    bookmarks_key = f"room:{room_id}:bookmarks"
    bookmarks_raw = r.hgetall(bookmarks_key)

    if not bookmarks_raw:
        return

    # Convert to int keys
    bookmarks = {int(k): v for k, v in bookmarks_raw.items()}

    # Build new bookmarks dict with shifted indices
    new_bookmarks = {}
    for idx, label in bookmarks.items():
        if idx >= insert_position:
            new_bookmarks[idx + 1] = label
        else:
            new_bookmarks[idx] = label

    # Update Redis with shifted bookmarks
    r.delete(bookmarks_key)
    if new_bookmarks:
        r.hset(bookmarks_key, mapping={str(k): v for k, v in new_bookmarks.items()})


def remove_bookmark_at_index(room_id: str, index: int):
    """
    Remove bookmark at a specific index (used when frame is replaced).

    Parameters
    ----------
    room_id : str
        The room ID
    index : int
        The frame index
    """
    r = current_app.extensions["redis"]
    bookmarks_key = f"room:{room_id}:bookmarks"
    r.hdel(bookmarks_key, str(index))


def emit_bookmarks_invalidate(room_id: str):
    """
    Emit bookmarks invalidate event to all clients after frame mapping changes.
    Clients will refetch bookmarks from the server when they receive this event.
    """
    socketio.emit(
        SocketEvents.INVALIDATE_BOOKMARK,
        {"operation": "frame_mapping_changed"},
        to=f"room:{room_id}",
    )


def emit_frames_invalidate(
    room_id: str,
    operation: str,
    affected_index: int | None = None,
    affected_from: int | None = None,
):
    """
    Emit frames:invalidate event to tell clients to clear their frame cache.

    Args:
        room_id: The room identifier
        operation: Type of operation - 'delete', 'insert', or 'replace'
        affected_index: For 'replace' - the specific frame index that was replaced
        affected_from: For 'delete'/'insert' - all frames from this index onward are affected
    """
    data = {
        "roomId": room_id,
        "operation": operation,
    }

    if affected_index is not None:
        data["affectedIndex"] = affected_index
    if affected_from is not None:
        data["affectedFrom"] = affected_from

    log.info(f"Emitting frames:invalidate for room '{room_id}': {data}")

    socketio.emit(
        "frames:invalidate",
        data,
        to=f"room:{room_id}",
    )


def emit_len_frames_update(room_id: str):
    """
    Emit len_frames event to update clients about the frame count after mutations.
    """
    room_service = current_app.extensions["room_service"]
    frame_count = room_service.get_frame_count(room_id)

    # Broadcast frame count update via room:update event
    from zndraw.app.room_manager import emit_room_update

    emit_room_update(socketio, room_id, frameCount=frame_count)


@main.route("/health")
def health_check():
    """Health check endpoint for server status verification."""
    return {"status": "ok"}, 200


@main.route("/api/version")
def get_version():
    """Get the ZnDraw server version."""
    import zndraw

    return {"version": zndraw.__version__}, 200


@main.route("/api/login", methods=["POST"])
def login():
    """Authenticate user and issue JWT token.

    Request
    -------
    {
        "userName": "John Doe"  // Required
    }

    Response
    --------
    {
        "status": "ok",
        "token": "eyJhbGc...",     // JWT token
        "clientId": "uuid-string"  // Server-generated client ID
    }
    """
    import datetime
    import uuid

    from zndraw.auth import create_jwt_token

    data = request.get_json() or {}
    user_name = data.get("userName")

    if not user_name or not user_name.strip():
        return {"error": "userName is required"}, 400

    # Generate server-side client ID
    client_id = str(uuid.uuid4())

    # Create JWT token
    token = create_jwt_token(client_id, user_name)

    # Store client metadata in Redis
    r = current_app.extensions["redis"]
    client_key = f"client:{client_id}"
    current_time = datetime.datetime.utcnow().isoformat()
    r.hset(
        client_key,
        mapping={
            "userName": user_name,
            "createdAt": current_time,
            "lastLogin": current_time,
        },
    )

    log.info(f"User '{user_name}' logged in with client ID: {client_id}")

    return {
        "status": "ok",
        "token": token,
        "clientId": client_id,
    }


@main.route("/assets/<path:filename>")
def serve_static_assets(filename: str):
    static_folder = Path(__file__).parent.parent / "static" / "assets"
    return send_from_directory(static_folder, filename)


@main.route("/", defaults={"path": ""})
@main.route("/<path:path>")
def serve_react_router_paths(path: str):
    """Catch-all route to serve index.html for React Router paths.

    This allows React Router to handle client-side routing for paths like /rooms, /rooms/:id, etc.
    API routes are registered with higher priority and won't match this pattern.
    """
    static_folder = Path(__file__).parent.parent / "static"
    return send_from_directory(static_folder, "index.html")


@main.route("/api/disconnect/<string:client_sid>", methods=["POST"])
def disconnect_sid(client_sid: str):
    """Disconnects the client from the room.

    Args:
        client_sid: Can be either a socket sid OR a client_id.
                   We try both lookups to support both cases.
    """
    try:
        r = current_app.extensions["redis"]

        # First, try to interpret client_sid as a client_id and get the socket sid
        socket_sid = r.hget(f"client:{client_sid}", "currentSid")

        if socket_sid:
            disconnect(socket_sid, namespace="/")
            return {"success": True}
        else:
            return {"success": False, "error": "Client not found or not connected"}
    except Exception as e:
        log.error(f"Error disconnecting client {client_sid}: {e}")
        return {"success": False, "error": str(e)}


@main.route("/internal/emit", methods=["POST"])
def internal_emit():
    """Internal endpoint to emit Socket.IO events. Secured via a shared secret."""
    data = request.get_json()

    event = data.get("event")
    sid = data.get("sid")
    payload = data.get("data", {})

    if not event or not sid:
        return {"error": "Event and sid are required"}, 400

    socketio.emit(event, payload, to=sid)
    return {"success": True}


@main.route("/api/rooms/<string:room_id>/frames", methods=["GET"])
def get_frames(room_id):
    """
    Serves multiple frames' data from the room's Zarr store using either
    indices or slice parameters. This implementation is optimized and uses
    the FrameIndexManager for all index operations.
    """
    r = current_app.extensions["redis"]
    try:
        indices_key = f"room:{room_id}:trajectory:indices"
        manager = FrameIndexManager(r, indices_key)

        frame_count = manager.get_count()

        if frame_count == 0:
            return {
                "error": f"Index out of range for data with {frame_count} frames in room '{room_id}'",
                "type": "IndexError",
            }, 404

        max_frame_idx = frame_count - 1
        physical_keys = []

        if "indices" in request.args:
            indices_str = request.args.get("indices")
            if not indices_str:
                frame_indices = []
            else:
                try:
                    frame_indices = [int(idx.strip()) for idx in indices_str.split(",")]
                except ValueError:
                    return {"error": "Indices must be comma-separated integers"}, 400

            for frame_id in frame_indices:
                if not (0 <= frame_id <= max_frame_idx):
                    return Response(
                        json.dumps(
                            {
                                "error": f"Invalid frame index {frame_id}, valid range: 0-{max_frame_idx}",
                                "type": "IndexError",
                            }
                        ),
                        status=404,
                        content_type="application/json",
                    )

            # Use the new manager method for non-contiguous indices
            physical_keys = manager.get_by_indices(frame_indices)

        else:  # Slice logic
            start_param = request.args.get("start")
            stop_param = request.args.get("stop")
            step_param = request.args.get("step")

            start = int(start_param) if start_param is not None else None
            stop = int(stop_param) if stop_param is not None else None
            step = int(step_param) if step_param is not None else 1

            if step is not None and step == 0:
                return {"error": "step cannot be zero"}, 400

            slice_obj = slice(start, stop, step)
            start, stop, step_val = slice_obj.indices(frame_count)

            if step_val == 1:
                # Use the optimized range method for contiguous slices
                physical_keys = manager.get_range(start, stop - 1)
            else:
                # For stepped slices, generate indices and use the by_indices method
                frame_indices = list(range(start, stop, step_val))
                physical_keys = manager.get_by_indices(frame_indices)

        # If no keys were found (e.g., empty indices list), return empty response
        if not physical_keys:
            return Response(msgpack.packb([]), content_type="application/octet-stream")

        keys_str = request.args.get("keys")
        requested_keys = [k.strip() for k in keys_str.split(",")] if keys_str else None

        frames_data = []
        try:
            for mapping_entry in physical_keys:
                mapping_entry_str = (
                    mapping_entry.decode("utf-8")
                    if isinstance(mapping_entry, bytes)
                    else mapping_entry
                )

                if ":" in mapping_entry_str:
                    source_room_id, physical_index_str = mapping_entry_str.split(":", 1)
                    physical_index = int(physical_index_str)
                    source_storage = get_storage(source_room_id)
                else:
                    physical_index = int(mapping_entry_str)
                    source_storage = get_storage(room_id)

                frame_data = source_storage.get(physical_index, keys=requested_keys)
                frames_data.append(encode_data(frame_data))
        except (KeyError, IndexError) as e:
            # Handle both Zarr key errors and physical index errors
            return Response(
                json.dumps(
                    {
                        "error": f"Error accessing physical storage: {e}",
                        "type": type(e).__name__,
                    }
                ),
                status=404,
                content_type="application/json",
            )

        packed_data = msgpack.packb(frames_data)
        return Response(packed_data, content_type="application/octet-stream")

    except Exception as e:
        log.error(f"Server error in get_frames: {e}\n{traceback.format_exc()}")
        return Response(
            json.dumps({"error": f"Server error: {e}", "type": type(e).__name__}),
            status=500,
            content_type="application/json",
        )


@main.route(
    "/api/rooms/<string:room_id>/frames/<int:frame_id>/metadata", methods=["GET"]
)
def get_frame_metadata(room_id: str, frame_id: int):
    """Get available keys and their shapes for a specific frame.

    Returns metadata about what data is available for the given frame,
    including the list of valid keys and their shapes/dtypes.
    """
    r = current_app.extensions["redis"]

    try:
        storage = get_storage(room_id)

        # Get logical-to-physical mapping from Redis
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {
                "error": f"No frames found in room '{room_id}'",
                "type": "IndexError",
            }, 404

        max_frame = len(frame_mapping) - 1

        # Validate frame_id
        if frame_id < 0 or frame_id > max_frame:
            return {
                "error": f"Invalid frame index {frame_id}, valid range: 0-{max_frame}",
                "type": "IndexError",
            }, 404

        # Get the physical index for this frame
        mapping_entry = frame_mapping[frame_id]

        if ":" not in mapping_entry:
            return {
                "error": f"Invalid frame mapping format for frame {frame_id}",
                "type": "ValueError",
            }, 500

        source_room_id, physical_index_str = mapping_entry.split(":", 1)
        physical_index = int(physical_index_str)
        source_storage = get_storage(source_room_id)

        # Get the Zarr group
        zarr_group = source_storage.group

        # Get valid keys for this frame
        if "__valid_keys__" not in zarr_group:
            return {
                "error": "No metadata available for this frame",
                "type": "MetadataError",
            }, 404

        try:
            valid_keys_array = zarr_group["__valid_keys__"]
            valid_keys_json = valid_keys_array[physical_index].item()  # type: ignore
            valid_keys = json.loads(valid_keys_json)
        except IndexError:
            return {
                "error": f"Frame index {frame_id} is out of bounds",
                "type": "IndexError",
            }, 404

        # Build metadata for each key
        keys_metadata = {}

        for key in valid_keys:
            if key in zarr_group:
                item = zarr_group[key]

                if isinstance(item, zarr.Array):
                    # Get shape and dtype for this array
                    is_json = item.attrs.get("format") == "json"

                    if is_json:
                        # For JSON-encoded arrays, we can't determine shape without loading
                        keys_metadata[key] = {"type": "json", "dtype": "json"}
                    else:
                        # Regular numpy array
                        shape = item.shape[1:]  # Remove the frame dimension

                        # Check if there's a mask for variable-sized arrays
                        mask_key = f"__mask__{key}__"
                        if mask_key in zarr_group:
                            mask_array = zarr_group[mask_key]
                            actual_size = int(mask_array[physical_index])  # type: ignore
                            shape = (actual_size,) + shape[1:]

                        keys_metadata[key] = {
                            "type": "array",
                            "shape": list(shape),
                            "dtype": str(item.dtype),
                        }
                elif isinstance(item, zarr.Group):
                    # For nested groups, indicate it's a group
                    keys_metadata[key] = {"type": "group"}

        return {
            "frameId": frame_id,
            # "physicalIndex": physical_index,
            "sourceRoom": source_room_id if ":" in mapping_entry else room_id,
            "keys": valid_keys,
            "metadata": keys_metadata,
        }, 200

    except Exception as e:
        error_data = {
            "error": f"Server error: {e}",
            "type": type(e).__name__,
        }
        log.error(f"Error getting frame metadata: {e}\n{traceback.format_exc()}")
        return Response(
            json.dumps(error_data), status=500, content_type="application/json"
        )


@main.route("/api/rooms/<string:room_id>/frames", methods=["DELETE"])
def delete_frames_batch(room_id):
    """Deletes frames using either a single frame_id, indices, or slice parameters from query params."""
    # Get client_id from query params (if provided)
    client_id = request.args.get("client_id")

    # Check if room is locked (passes client_id for lock holder check)
    lock_error = check_room_locked(room_id, client_id)
    if lock_error:
        return lock_error

    r = current_app.extensions["redis"]

    try:
        indices_key = f"room:{room_id}:trajectory:indices"
        index_manager = FrameIndexManager(r, indices_key)
        frame_mapping = index_manager.get_all()

        if not frame_mapping:
            return {"error": "No frames found in room"}, 404

        max_frame = len(frame_mapping) - 1

        # Determine frame indices based on request parameters
        if "frame_id" in request.args:
            # Single frame deletion
            frame_id_str = request.args.get("frame_id")
            try:
                frame_id = int(frame_id_str)
            except ValueError:
                return {"error": "frame_id must be an integer"}, 400

            if frame_id < 0 or frame_id > max_frame:
                error_data = {
                    "error": f"Invalid frame index {frame_id}, valid range: 0-{max_frame}",
                    "type": "IndexError",
                }
                return Response(
                    json.dumps(error_data),
                    status=404,
                    content_type="application/json",
                )
            frame_indices = [frame_id]

        elif "indices" in request.args:
            # Indices parameter was provided as comma-separated values
            indices_str = request.args.get("indices")

            # Handle empty string case
            if not indices_str:
                frame_indices = []
            else:
                # Split by comma and convert to integers
                try:
                    frame_indices = [int(idx.strip()) for idx in indices_str.split(",")]
                except ValueError:
                    return {"error": "Indices must be comma-separated integers"}, 400

            # Validate frame indices
            for frame_id in frame_indices:
                if (
                    not isinstance(frame_id, int)
                    or frame_id < 0
                    or frame_id > max_frame
                ):
                    error_data = {
                        "error": f"Invalid frame index {frame_id}, valid range: 0-{max_frame}",
                        "type": "IndexError",
                    }
                    return Response(
                        json.dumps(error_data),
                        status=404,
                        content_type="application/json",
                    )

        else:
            # Use slice parameters
            start_param = request.args.get("start")
            stop_param = request.args.get("stop")
            step_param = request.args.get("step")

            start = int(start_param) if start_param is not None else None
            stop = int(stop_param) if stop_param is not None else None
            step = int(step_param) if step_param is not None else None

            if step == 0:
                return {"error": "step cannot be zero"}, 400

            # Use slice.indices to properly handle negative indices
            slice_obj = slice(start, stop, step)
            start, stop, step = slice_obj.indices(len(frame_mapping))

            # Generate frame indices from slice
            frame_indices = list(range(start, stop, step))

        # Check for template frames before deleting
        for frame_id in frame_indices:
            mapping_entry = frame_mapping[frame_id]
            if isinstance(mapping_entry, bytes):
                mapping_entry = mapping_entry.decode()

            # Check if this frame belongs to a different room (e.g., a template)
            if ":" in mapping_entry:
                source_room_id = mapping_entry.split(":", 1)[0]
                if source_room_id != room_id:
                    return {
                        "error": f"Cannot delete template frame at index {frame_id}. This frame belongs to '{source_room_id}'",
                        "type": "PermissionError",
                    }, 403

        # Delete frames using gap method - no need to rebuild the entire sorted set
        members_to_delete = [frame_mapping[idx] for idx in frame_indices]

        if members_to_delete:
            # Batch delete all members at once
            r.zrem(indices_key, *members_to_delete)

        log.info(
            f"Deleted {len(frame_indices)} frames from room '{room_id}'. Physical data preserved."
        )

        # Shift bookmarks after deletion
        shift_bookmarks_on_delete(room_id, frame_indices)
        # Emit bookmarks update to all clients
        emit_bookmarks_invalidate(room_id)
        # Invalidate all frames from the first deleted position onward
        if frame_indices:
            emit_frames_invalidate(
                room_id, operation="delete", affected_from=min(frame_indices)
            )
        # Emit len_frames update to notify clients of frame count change
        emit_len_frames_update(room_id)

        return {"success": True, "deleted_count": len(frame_indices)}
    except Exception as e:
        error_data = {
            "error": f"Server error: {e}",
            "type": type(e).__name__,
            "success": False,
        }
        log.error(f"Server error: {e}\n{traceback.format_exc()}")
        return Response(
            json.dumps(error_data), status=500, content_type="application/json"
        )


@main.route("/api/rooms/<string:room_id>/frames", methods=["POST"])
def append_frame(room_id):
    """Handles frame operations (append, extend, replace, insert) based on action query parameter."""
    # Get client_id from query params (if provided)
    client_id = request.args.get("client_id")

    # Check if room is locked (passes client_id for lock holder check)
    lock_error = check_room_locked(room_id, client_id)
    if lock_error:
        return lock_error

    r = current_app.extensions["redis"]

    # Get action from query parameters
    action = request.args.get("action", "append")

    if action not in {"append", "replace", "insert", "extend"}:
        return {"error": "Invalid action specified"}, 400

    try:
        # Unpack the msgpack data
        serialized_data = msgpack.unpackb(request.data, strict_map_key=False)

        storage = get_storage(room_id)

        indices_key = f"room:{room_id}:trajectory:indices"
        index_manager = FrameIndexManager(r, indices_key)

        if action == "replace":
            # Get frame_id from query parameters
            target_frame_id = request.args.get("frame_id")
            if target_frame_id is None:
                return {"error": "frame_id is required for replace operations"}, 400

            try:
                target_frame_id = int(target_frame_id)
            except ValueError:
                return {"error": "frame_id must be an integer"}, 400

            frame_mapping = r.zrange(indices_key, 0, -1)
            frame_mapping_with_scores = r.zrange(indices_key, 0, -1, withscores=True)

            # Validate target_frame_id
            if not (0 <= target_frame_id < len(frame_mapping)):
                return {
                    "error": f"Invalid or missing frame_id for replace. Valid range: 0-{len(frame_mapping) - 1}",
                    "type": "IndexError",
                }, 404

            new_physical_index = len(storage)
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            old_mapping_entry = frame_mapping[target_frame_id]
            # Get the original score to preserve gap-based indexing
            _, old_score = frame_mapping_with_scores[target_frame_id]

            pipeline = r.pipeline()
            pipeline.zrem(indices_key, old_mapping_entry)
            pipeline.zadd(indices_key, {f"{room_id}:{new_physical_index}": old_score})
            pipeline.execute()

            # Remove bookmark at the logical index (if it exists)
            remove_bookmark_at_index(room_id, target_frame_id)

            # Emit bookmarks update to reflect removal
            emit_bookmarks_invalidate(room_id)
            # Invalidate only the replaced frame
            emit_frames_invalidate(
                room_id, operation="replace", affected_index=target_frame_id
            )

            log.info(
                f"Replaced frame {target_frame_id} (old: {old_mapping_entry}, new: {room_id}:{new_physical_index}) in room '{room_id}'"
            )
            return {"success": True, "replaced_frame": target_frame_id}

        elif action == "extend":
            # Extend operation: add multiple frames in one go
            if not isinstance(serialized_data, list):
                return {
                    "error": "For extend action, data must be a list of frame dictionaries"
                }, 400

            # 1. Determine starting logical and physical positions
            start_logical_pos = index_manager.get_count()
            start_physical_pos = len(storage)
            num_frames = len(serialized_data)

            # 2. Decode all frames and extend the physical storage
            decoded_frames = [decode_data(frame) for frame in serialized_data]
            storage.extend(decoded_frames)

            # 3. Append all new frames using gap method
            for i in range(num_frames):
                index_manager.append(f"{room_id}:{start_physical_pos + i}")

            # 4. Prepare response data
            new_indices = list(range(start_logical_pos, start_logical_pos + num_frames))

            # Emit len_frames update to notify clients of frame count change
            emit_len_frames_update(room_id)

            log.info(
                f"Extended trajectory with {num_frames} frames (physical: {start_physical_pos}-{start_physical_pos + num_frames - 1}) to room '{room_id}'"
            )
            return {"success": True, "new_indices": new_indices}

        elif action == "insert":
            # Get insert_position from query parameters
            insert_position_str = request.args.get("insert_position")
            if insert_position_str is None:
                return {
                    "error": "insert_position is required for insert operations"
                }, 400

            try:
                insert_position = int(insert_position_str)
            except ValueError:
                return {"error": "insert_position must be an integer"}, 400

            current_length = index_manager.get_count()

            if not (0 <= insert_position <= current_length):
                return {
                    "error": f"Insert position {insert_position} out of range [0, {current_length}]"
                }, 400

            # 1. Determine the new physical index
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data to the physical store
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Insert using gap method (O(log N) instead of O(N))
            index_manager.insert(insert_position, f"{room_id}:{new_physical_index}")

            # Shift bookmarks after insertion
            shift_bookmarks_on_insert(room_id, insert_position)
            # Emit bookmarks update (logical indices shifted by the insert)
            emit_bookmarks_invalidate(room_id)
            # Invalidate all frames from insert position onward (they all shift up)
            emit_frames_invalidate(
                room_id, operation="insert", affected_from=insert_position
            )
            # Emit len_frames update to notify clients of frame count change
            emit_len_frames_update(room_id)

            log.info(
                f"Inserted frame at position {insert_position} (physical: {new_physical_index}) in room '{room_id}'"
            )
            return {"success": True, "inserted_position": insert_position}

        elif action == "append":
            # Append operation: add a new frame to the end of the logical sequence
            # 1. Determine logical and physical positions
            logical_position = index_manager.get_count()
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Add the new physical index to the logical sequence using gap method
            index_manager.append(f"{room_id}:{new_physical_index}")

            # Emit len_frames update to notify clients of frame count change
            emit_len_frames_update(room_id)

            log.debug(
                f"Appended frame {logical_position} (physical: {room_id}:{new_physical_index}) to room '{room_id}'"
            )
            return {"success": True, "new_index": logical_position}

        else:
            # Default case for any unknown actions
            return {"error": f"The requested action '{action}' is not supported."}, 400
    except Exception as e:
        log.error(f"Failed to write to Zarr store: {e}\n{traceback.format_exc()}")
        return {"error": "Failed to write to data store"}, 500


@main.route("/api/rooms/<string:room_id>/frames/bulk", methods=["PATCH"])
def bulk_replace_frames(room_id):
    """
    Bulk replace frames using either a slice (start/stop) or a list of indices.

    This endpoint is optimized for replacing multiple frames in a single operation,
    significantly reducing the number of HTTP requests and Redis operations.

    Query parameters:
    - For slice operations: start, stop (optional)
    - For list operations: indices (comma-separated)

    Body: msgpack-encoded list of frame dictionaries
    """
    r = current_app.extensions["redis"]

    try:
        # Unpack the msgpack data - should be a list of frames
        serialized_data = msgpack.unpackb(request.data, strict_map_key=False)

        if not isinstance(serialized_data, list):
            return {"error": "Body must contain a list of frame dictionaries"}, 400

        storage = get_storage(room_id)
        indices_key = f"room:{room_id}:trajectory:indices"
        index_manager = FrameIndexManager(r, indices_key)
        frame_mapping = index_manager.get_all()
        frame_mapping_with_scores = index_manager.get_all(withscores=True)

        # Determine which indices to replace
        if "indices" in request.args:
            # List operations require existing frames
            if not frame_mapping:
                return {"error": "No frames found in room"}, 404
            # List of specific indices
            indices_str = request.args.get("indices")
            if not indices_str:
                return {"error": "indices parameter cannot be empty"}, 400

            try:
                target_indices = [int(idx.strip()) for idx in indices_str.split(",")]
            except ValueError:
                return {"error": "Indices must be comma-separated integers"}, 400

            # Validate indices
            for idx in target_indices:
                if idx < 0 or idx >= len(frame_mapping):
                    return {
                        "error": f"Invalid index {idx}, valid range: 0-{len(frame_mapping) - 1}",
                        "type": "IndexError",
                    }, 404

            # Check length matches
            if len(serialized_data) != len(target_indices):
                return {
                    "error": f"Number of frames ({len(serialized_data)}) must match number of indices ({len(target_indices)})"
                }, 400

            # Replace each frame
            start_physical_pos = len(storage)
            decoded_frames = [decode_data(frame) for frame in serialized_data]
            storage.extend(decoded_frames)

            pipeline = r.pipeline()
            for i, logical_idx in enumerate(target_indices):
                old_mapping_entry = frame_mapping[logical_idx]
                # Get the original score to preserve gap-based indexing
                old_member, old_score = frame_mapping_with_scores[logical_idx]
                new_physical_idx = start_physical_pos + i

                # Remove old mapping and add new one with same score
                pipeline.zrem(indices_key, old_mapping_entry)
                pipeline.zadd(indices_key, {f"{room_id}:{new_physical_idx}": old_score})

            pipeline.execute()

            emit_bookmarks_invalidate(room_id)
            emit_frames_invalidate(
                room_id, operation="replace", affected_from=min(target_indices)
            )

            log.info(f"Bulk replaced {len(target_indices)} frames in room '{room_id}'")
            return {"success": True, "replaced_count": len(target_indices)}

        else:
            # Slice operation
            start_param = request.args.get("start")
            stop_param = request.args.get("stop")

            start = int(start_param) if start_param is not None else None
            stop = int(stop_param) if stop_param is not None else None

            # Use slice.indices to handle negative indices and defaults
            slice_obj = slice(start, stop, 1)  # step is always 1 for simple slices
            start, stop, _ = slice_obj.indices(len(frame_mapping))

            # Calculate how many frames are being replaced
            old_count = stop - start
            new_count = len(serialized_data)

            # Delete the old frames and insert the new ones
            # This is more complex because the length can change

            # 1. Decode all new frames and add to physical storage
            start_physical_pos = len(storage)
            decoded_frames = [decode_data(frame) for frame in serialized_data]
            storage.extend(decoded_frames)

            # 2. Build the new frame mapping
            # Keep frames before the slice
            new_mapping_list = list(frame_mapping[:start])

            # Add new frames
            for i in range(new_count):
                new_mapping_list.append(f"{room_id}:{start_physical_pos + i}")

            # Add frames after the slice (if any)
            new_mapping_list.extend(frame_mapping[stop:])

            # 3. Rebuild the entire Redis sorted set with the new mapping
            pipeline = r.pipeline()
            pipeline.delete(indices_key)

            new_mapping = {
                physical_idx_str: logical_pos
                for logical_pos, physical_idx_str in enumerate(new_mapping_list)
            }
            if new_mapping:
                pipeline.zadd(indices_key, new_mapping)

            pipeline.execute()

            emit_bookmarks_invalidate(room_id)
            emit_len_frames_update(room_id)
            emit_frames_invalidate(
                room_id, operation="bulk_replace", affected_from=start
            )

            log.info(
                f"Bulk replaced slice [{start}:{stop}] ({old_count} frames) with {new_count} frames in room '{room_id}'"
            )
            return {
                "success": True,
                "replaced_range": [start, stop],
                "old_count": old_count,
                "new_count": new_count,
                "new_length": len(new_mapping_list),
            }

    except Exception as e:
        log.error(f"Failed to bulk replace frames: {e}\n{traceback.format_exc()}")
        return {"error": "Failed to bulk replace frames"}, 500


@main.route("/api/rooms", methods=["GET"])
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
        # Get frame count using service
        frame_count = room_service.get_frame_count(room_id)

        # Get metadata
        description = redis_client.get(f"room:{room_id}:description")
        locked = redis_client.get(f"room:{room_id}:locked") == "1"
        hidden = redis_client.get(f"room:{room_id}:hidden") == "1"
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


@main.route("/api/rooms/<string:room_id>", methods=["GET"])
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
    description = redis_client.get(f"room:{room_id}:description")
    locked = redis_client.get(f"room:{room_id}:locked") == "1"
    hidden = redis_client.get(f"room:{room_id}:hidden") == "1"

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


@main.route("/api/rooms/<string:room_id>", methods=["PATCH"])
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
            redis_client.delete(f"room:{room_id}:description")
            changes["description"] = None
        else:
            redis_client.set(f"room:{room_id}:description", data["description"])
            changes["description"] = data["description"]

    # Update locked status
    if "locked" in data:
        redis_client.set(f"room:{room_id}:locked", "1" if data["locked"] else "0")
        changes["locked"] = bool(data["locked"])

    # Update hidden status
    if "hidden" in data:
        redis_client.set(f"room:{room_id}:hidden", "1" if data["hidden"] else "0")
        changes["hidden"] = bool(data["hidden"])

    # Emit socket event for real-time updates
    if changes:
        from zndraw.app.room_manager import emit_room_update

        emit_room_update(socketio, room_id, **changes)
        log.debug(f"Emitted room:update for room '{room_id}': {changes}")

    log.info(f"Updated room '{room_id}' metadata: {data}")
    return {"status": "ok"}, 200


@main.route("/api/rooms/default", methods=["GET"])
def get_default_room():
    """Get the default room ID.

    Returns:
        {"roomId": "room1"} or {"roomId": null}
    """
    redis_client = current_app.extensions["redis"]
    default_room = redis_client.get("default_room")
    return {"roomId": default_room if default_room else None}, 200


@main.route("/api/rooms/default", methods=["PUT"])
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


@main.route("/api/rooms/<string:room_id>/metadata", methods=["GET"])
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


@main.route("/api/rooms/<string:room_id>/metadata", methods=["POST"])
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
    locked = redis_client.get(f"room:{room_id}:locked")
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


@main.route("/api/rooms/<string:room_id>/metadata/<string:field>", methods=["DELETE"])
def delete_room_metadata_field(room_id: str, field: str):
    """Delete specific metadata field. Respects room lock.

    Returns:
        {"success": true}
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    # Check permanent room lock only (metadata doesn't use trajectory lock)
    redis_client = current_app.extensions["redis"]
    locked = redis_client.get(f"room:{room_id}:locked")
    if locked == "1":
        return {"error": "Room is locked and cannot be modified"}, 403

    manager = RoomMetadataManager(redis_client, room_id)

    try:
        deleted = manager.delete(field)
        return {"success": True, "deleted": deleted}, 200
    except Exception as e:
        log.error(f"Error deleting metadata field '{field}' for room '{room_id}': {e}")
        return {"error": str(e)}, 500


@main.route("/api/rooms/<string:room_id>/locks/<string:target>", methods=["GET"])
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
            "holder": "client_123",
            "metadata": {
                "clientId": "client_123",
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


@main.route("/api/rooms/<string:room_id>/duplicate", methods=["POST"])
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
    redis_client = current_app.extensions["redis"]
    data = request.get_json() or {}

    # Check source room exists
    source_indices_key = f"room:{room_id}:trajectory:indices"
    if not redis_client.exists(source_indices_key):
        return {"error": "Source room not found"}, 404

    # Generate or use provided new room ID
    new_room_id = data.get("newRoomId")
    if not new_room_id:
        import uuid

        new_room_id = str(uuid.uuid4())

    # Check new room doesn't already exist
    if redis_client.exists(f"room:{new_room_id}:trajectory:indices"):
        return {"error": "Room with that ID already exists"}, 409

    # 1. Copy trajectory indices (sorted set) - this shares frame data
    source_indices = redis_client.zrange(source_indices_key, 0, -1, withscores=True)
    if source_indices:
        redis_client.zadd(
            f"room:{new_room_id}:trajectory:indices",
            {member: score for member, score in source_indices},
        )

    # 2. Copy geometries hash
    geometries = redis_client.hgetall(f"room:{room_id}:geometries")
    if geometries:
        redis_client.hset(f"room:{new_room_id}:geometries", mapping=geometries)

    # 3. Copy bookmarks hash (physical keys remain valid)
    bookmarks = redis_client.hgetall(f"room:{room_id}:bookmarks")
    if bookmarks:
        redis_client.hset(f"room:{new_room_id}:bookmarks", mapping=bookmarks)

    # 4. Initialize new room metadata
    redis_client.set(f"room:{new_room_id}:current_frame", 0)
    redis_client.set(f"room:{new_room_id}:locked", 0)
    redis_client.set(f"room:{new_room_id}:hidden", 0)

    description = data.get("description")
    if description:
        redis_client.set(f"room:{new_room_id}:description", description)

    # 5. Initialize default geometries for new room
    from zndraw.geometries import Bond, Cell, Curve, Floor, Sphere

    if not geometries:  # Only if source had no geometries
        # Create particles with explicit arrays.positions/colors for initial /join
        particles_data = Sphere(
            position="arrays.positions", color="arrays.colors"
        ).model_dump()
        redis_client.hset(
            f"room:{new_room_id}:geometries",
            "particles",
            json.dumps({"type": Sphere.__name__, "data": particles_data}),
        )
        redis_client.hset(
            f"room:{new_room_id}:geometries",
            "bonds",
            json.dumps({"type": Bond.__name__, "data": Bond().model_dump()}),
        )
        redis_client.hset(
            f"room:{new_room_id}:geometries",
            "curve",
            json.dumps({"type": Curve.__name__, "data": Curve().model_dump()}),
        )
        redis_client.hset(
            f"room:{new_room_id}:geometries",
            "cell",
            json.dumps({"type": Cell.__name__, "data": Cell().model_dump()}),
        )
        redis_client.hset(
            f"room:{new_room_id}:geometries",
            "floor",
            json.dumps({"type": Floor.__name__, "data": Floor().model_dump()}),
        )

    log.info(
        f"Duplicated room '{room_id}' to '{new_room_id}' with {len(source_indices)} frames"
    )

    # Emit room:update event to notify clients of new room
    from zndraw.app.room_manager import emit_room_update

    emit_room_update(
        socketio,
        new_room_id,
        created=True,
        description=description,
        frameCount=len(source_indices),
        locked=False,
        hidden=False,
        isDefault=False,
    )

    return {
        "status": "ok",
        "roomId": new_room_id,
        "frameCount": len(source_indices),
    }, 200


@main.route("/api/rooms/<string:room_id>/renormalize", methods=["POST"])
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
    indices_key = f"room:{room_id}:trajectory:indices"

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


@main.route("/api/shutdown", methods=["POST"])
def exit_app():
    """Endpoint to gracefully shut down the server. Secured via a shared secret."""
    socketio.stop()
    return {"success": True}  # this might never be seen


@main.route("/api/rooms/<string:room_id>/schema/<string:category>", methods=["GET"])
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

    # Add client-provided extensions from Redis
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

    return schema


@main.route(
    "/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/submit",
    methods=["POST"],
)
def log_room_extension(room_id: str, category: str, extension: str):
    """Submits a user extension action to create a job."""
    from zndraw.auth import AuthError, get_current_client

    # Authenticate and get user from JWT token
    try:
        client = get_current_client()
        user_id = client["userName"]
    except AuthError as e:
        return {"error": e.message}, e.status_code

    json_data = request.json
    if json_data is None:
        json_data = {}

    data = json_data.pop("data", None)
    if data is None:
        # If no "data" key, treat the entire JSON body as data
        data = json_data

    log.info(
        f"Logging extension for room {room_id}: category={category}, extension={extension}, data={json.dumps(data)}"
    )

    # store in redis
    redis_client = current_app.extensions["redis"]

    # Check if this is a server-side (Celery) extension
    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "analysis": analysis,
        "settings": settings,
    }

    is_celery_extension = (
        category in category_map and extension in category_map[category]
    )

    # Check if extension exists (either as server-side or client-registered)
    if not is_celery_extension:
        # Check if any client has registered this extension
        schema_key = ExtensionKeys.schema_key(room_id, category)
        extension_schema = redis_client.hget(schema_key, extension)

        if extension_schema is None:
            return {"error": f"No workers available for extension {extension}"}, 400

    # Store the entire extension data as a JSON string
    redis_client.hset(
        f"room:{room_id}:user:{user_id}:{category}", extension, json.dumps(data)
    )

    # Handle settings differently - no job queue, just update and notify
    if category == "settings":
        # Emit socket event to notify all clients in room
        socketio.emit(
            SocketEvents.INVALIDATE,
            {
                "userId": user_id,
                "category": category,
                "extension": extension,
                "roomId": room_id,
            },
            to=f"room:{room_id}",  # Broadcast to all users in room
        )
        log.info(
            f"Updated settings for room {room_id}, user {user_id}: {extension} = {json.dumps(data)}"
        )
        return {"status": "success", "message": "Settings updated"}, 200

    # Create job
    provider = "celery" if is_celery_extension else "client"
    job_id = JobManager.create_job(
        redis_client, room_id, category, extension, data, user_id, provider
    )

    queue_position = 0

    if is_celery_extension:
        # Queue job for Celery workers to poll via /jobs/next endpoint
        # Celery workers will check provider="celery" and pick up these jobs
        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add to queue with provider info
        redis_client.rpush(
            keys.queue,
            json.dumps(
                {
                    "user_id": user_id,
                    "data": data,
                    "room": room_id,
                    "jobId": job_id,
                    "provider": "celery",
                }
            ),
        )
        log.info(
            f"Queued Celery task for user {user_id}, category {category}, extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1

        # Trigger a celery worker task to pick up the job
        from zndraw.app.tasks import celery_job_worker

        server_url = current_app.config.get("SERVER_URL", "http://localhost:5000")
        _ = celery_job_worker.delay(room_id, server_url)

        # Notify all clients in room about queue update
        emit_queue_update(redis_client, room_id, category, extension, socketio)
    else:
        # Queue job for client workers to poll via /jobs/next endpoint
        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add to queue
        redis_client.rpush(
            keys.queue,
            json.dumps(
                {"user_id": user_id, "data": data, "room": room_id, "jobId": job_id}
            ),
        )
        log.info(
            f"Queued task for user {user_id}, category {category}, extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1  # Zero-indexed position

        # Notify all clients in room about queue update
        emit_queue_update(redis_client, room_id, category, extension, socketio)

    log.info(
        f"Emitting invalidate for user {user_id}, category {category}, extension {extension}, room {room_id} to user:{user_id}"
    )
    socketio.emit(
        SocketEvents.INVALIDATE,
        {
            "userId": user_id,
            "category": category,
            "extension": extension,
            "roomId": room_id,
        },
        to=f"user:{user_id}",
    )
    return {"status": "success", "queuePosition": queue_position, "jobId": job_id}, 200


@main.route(
    "/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/data",
    methods=["GET"],
)
def get_extension_data(room_id: str, category: str, extension: str):
    """Gets the cached data for a user's extension."""
    from zndraw.auth import AuthError, get_current_client

    # Authenticate and get user from JWT token
    try:
        client = get_current_client()
        user_id = client["userName"]
    except AuthError as e:
        return {"error": e.message}, e.status_code

    print(
        f"get_extension_data called with userId={user_id}, category={category}, extension={extension} for room {room_id}"
    )

    redis_client = current_app.extensions["redis"]
    extension_data = redis_client.hget(
        f"room:{room_id}:user:{user_id}:{category}", extension
    )
    if extension_data is None:
        return {"data": None}, 200
    extension_data = json.loads(extension_data)
    return {"data": extension_data}, 200


@main.route(
    "/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/workers",
    methods=["GET"],
)
def get_extension_workers(room_id: str, category: str, extension: str):
    """Get worker details for a specific extension.

    Returns:
        {
            "idleWorkers": ["worker_id_1", "worker_id_2"],
            "progressingWorkers": ["worker_id_3"],
            "queueLength": 5
        }
    """
    redis_client = current_app.extensions["redis"]
    keys = ExtensionKeys.for_extension(room_id, category, extension)

    idle_workers = list(redis_client.smembers(keys.idle_workers))
    progressing_workers = list(redis_client.smembers(keys.progressing_workers))
    queue_length = redis_client.llen(keys.queue)

    return {
        "idleWorkers": idle_workers,
        "progressingWorkers": progressing_workers,
        "queueLength": queue_length,
        "totalWorkers": len(idle_workers) + len(progressing_workers),
    }, 200


# Job management endpoints
@main.route("/api/rooms/<string:room_id>/jobs", methods=["GET"])
def list_jobs(room_id: str):
    """List active jobs for a room."""
    redis_client = current_app.extensions["redis"]
    jobs = JobManager.list_all_jobs(redis_client, room_id)
    return jobs, 200


@main.route(
    "/api/rooms/<string:room_id>/jobs/<string:job_id>", methods=["GET", "DELETE"]
)
def get_job(room_id: str, job_id: str):
    """Get details for a specific job."""
    redis_client = current_app.extensions["redis"]
    if request.method == "DELETE":
        job = JobManager.get_job(redis_client, job_id)
        if not job:
            return {"error": "Job not found"}, 404
        if job.get("status") == "running":
            return {"error": "Cannot delete a running job"}, 400

        category = job.get("category")
        extension = job.get("extension")

        JobManager.delete_job(redis_client, job_id)

        # Emit queue update instead of job:deleted
        if category and extension:
            emit_queue_update(redis_client, room_id, category, extension, socketio)

        return {"status": "success"}, 200

    job = JobManager.get_job(redis_client, job_id)
    if not job:
        return {"error": "Job not found"}, 404
    return job, 200


@main.route("/api/rooms/<string:room_id>/jobs/<string:job_id>/status", methods=["PUT"])
def update_job_status(room_id: str, job_id: str):
    """Update a job's status (complete or fail) and transition worker back to idle.

    Called by both Celery tasks and client workers when they finish.
    For client workers: transitions worker from progressing → idle, then checks queue.
    For Celery workers: just marks job complete/failed (no worker state to manage).

    Request body:
        {
            "status": "completed" | "failed",
            "result": <result data>,  // for completed status
            "error": <error message>,  // for failed status
            "workerId": <worker_id>  // optional, only for client workers
        }
    """
    data = request.get_json() or {}
    status = data.get("status")
    worker_id = data.get("workerId")  # Optional: only for client workers

    if status not in ["completed", "failed"]:
        return {"error": "Status must be 'completed' or 'failed'"}, 400

    redis_client = current_app.extensions["redis"]

    # Get job details to know category/extension
    job = JobManager.get_job(redis_client, job_id)
    if not job:
        log.error(f"Job {job_id} not found in Redis")
        return {"error": "Job not found"}, 404

    log.info(
        f"Update job status called: job_id={job_id}, status={status}, worker_id={worker_id}, room_id={room_id}"
    )
    log.info(
        f"Job data: category={job.get('category')}, extension={job.get('extension')}"
    )

    # Validate job is in running state
    if job.get("status") != "running":
        log.error(f"Job {job_id} is not running (status: {job.get('status')})")
        return {"error": "Job is not running"}, 400

    # Validate worker ID matches the job's assigned worker
    if worker_id and job.get("worker_id") != worker_id:
        log.error(
            f"Worker ID mismatch: job assigned to {job.get('worker_id')}, but {worker_id} tried to update it"
        )
        return {"error": "Worker ID does not match job's worker ID"}, 400

    # Update job status based on requested status
    if status == "completed":
        result = data.get("result")
        success = JobManager.complete_job(redis_client, job_id, result)
        if not success:
            log.error(f"Failed to mark job {job_id} as completed")
            return {"error": "Failed to complete job"}, 400
        log.info(f"Job {job_id} completed in room {room_id}")
        worker_success = True
    else:  # status == "failed"
        error = data.get("error", "Unknown error")
        success = JobManager.fail_job(redis_client, job_id, error)
        if not success:
            return {"error": "Failed to mark job as failed"}, 400
        log.error(f"Job {job_id} failed in room {room_id}: {error}")
        worker_success = False

    # If this is a client worker (not Celery), handle worker state transition
    if worker_id and not worker_id.startswith("celery:"):
        log.info(f"Transitioning worker {worker_id} to idle")
        _transition_worker_to_idle(
            redis_client, socketio, worker_id, job, room_id, success=worker_success
        )
    else:
        log.info(f"Skipping worker transition (worker_id={worker_id})")

    return {"status": "success"}, 200


@main.route("/api/rooms/<string:room_id>/extensions/register", methods=["POST"])
def register_extension(room_id: str):
    data = request.get_json()
    try:
        name = data["name"]
        category = data["category"]
        schema = data["schema"]
        client_id = data["clientId"]
    except KeyError as e:
        return {"error": f"Missing required field: {e}"}, 400

    redis_client = current_app.extensions["redis"]

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

    if category in category_map and name in category_map[category]:
        log.warning(
            f"Blocked attempt to register extension '{name}' in category '{category}' "
            f"- name conflicts with server-side extension (security violation)"
        )
        return {
            "error": f"Cannot register extension '{name}': name is reserved for server-side extensions"
        }

    log.info(
        f"Registering extension for room {room_id}: name={name}, category={category}"
    )

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
            }, 400
        redis_client.sadd(keys.idle_workers, client_id)
        redis_client.sadd(user_extensions_key, name)  # Keep this for disconnect cleanup

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
        # if room_id:  # Null check for type safety
        #     dispatch_next_task(redis_client, socketio, client_id, room_id, category)

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
        # if room_id:  # Null check for type safety
        #     dispatch_next_task(redis_client, socketio, client_id, room_id, category)

    return {"status": "success"}


@main.route("/api/workers/<string:worker_id>", methods=["GET"])
def get_worker_state(worker_id: str):
    """Get the current state of a worker.

    Returns:
        {
            "idle": true/false,
            "currentJob": job_id or null
        }
    """
    redis_client = current_app.extensions["redis"]

    # Check both modifiers and selections categories to find if worker is progressing
    current_job_id = None
    is_idle = True

    for category in ["modifiers", "selections", "analysis"]:
        # Get all extensions this worker might be registered for
        # We need to check all extensions to find if the worker is in any progressing set
        # First, get the worker's registered extensions
        user_extensions_key = ExtensionKeys.user_extensions_key(
            "*", category, worker_id
        )

        # Since we don't know the room, we need to scan for worker state
        # Check if worker has a current job by scanning all rooms
        # For simplicity, we'll check the most common pattern

        # Alternative approach: check if there's any job assigned to this worker
        # by looking for jobs where worker_id matches
        for key in redis_client.scan_iter(match="job:*"):
            job_data = redis_client.hgetall(key)
            if (
                job_data.get("worker_id") == worker_id
                and job_data.get("status") == "running"
            ):
                current_job_id = job_data.get("id")
                is_idle = False
                break

        if current_job_id:
            break

    return {"idle": is_idle, "currentJob": current_job_id}, 200


@main.route("/api/rooms/<string:room_id>/jobs/next", methods=["POST"])
def get_next_job(room_id: str):
    """Poll for the next available job for a worker.

    Workers call this endpoint to get assigned their next job from the queue.
    The worker should be registered for at least one extension category.

    Request body:
        {
            "workerId": "worker_session_id"
        }

    Returns:
        Job object directly (jobId, category, extension, data, etc.) or null if no jobs available
    """
    try:
        data = request.get_json() or {}
        worker_id = data.get("workerId")

        if not worker_id:
            return {"error": "workerId is required"}, 400

        redis_client = current_app.extensions["redis"]

        # Check if worker already has a running job
        for key in redis_client.scan_iter(match="job:*"):
            job_data = redis_client.hgetall(key)
            if (
                job_data.get("worker_id") == worker_id
                and job_data.get("status") == "running"
            ):
                return {"error": "Worker is not idle"}, 400

        # If worker is "celery-worker", check for celery jobs across all extensions
        is_celery_worker = worker_id.startswith("celery")

        # Check both modifiers and selections categories for queued jobs
        for category in ["modifiers", "selections", "analysis"]:
            if is_celery_worker:
                # For celery workers, check all extensions for celery jobs
                try:
                    from zndraw.extensions.analysis import analysis
                    from zndraw.extensions.modifiers import modifiers
                    from zndraw.extensions.selections import selections
                except ImportError as e:
                    log.error(f"Failed to import extensions: {e}")
                    continue

                category_map = {
                    "modifiers": modifiers,
                    "selections": selections,
                    "analysis": analysis,
                }

                if category in category_map:
                    for extension in category_map[category].keys():
                        keys = ExtensionKeys.for_extension(room_id, category, extension)
                        queue_length = redis_client.llen(keys.queue)

                        if queue_length > 0:
                            # Peek at the first job to check if it's a celery job
                            task_data = redis_client.lindex(keys.queue, 0)
                            if task_data:
                                try:
                                    task_info = json.loads(task_data)
                                except json.JSONDecodeError as e:
                                    log.error(
                                        f"Invalid JSON in queue {keys.queue}: {e}, data: {task_data}"
                                    )
                                    # Skip this malformed job and continue
                                    redis_client.lpop(keys.queue)
                                    continue

                                if task_info.get("provider") == "celery":
                                    # This is a celery job, pop it
                                    redis_client.lpop(keys.queue)
                                    job_id = task_info.get("jobId")

                                    # Get full job details
                                    job = JobManager.get_job(redis_client, job_id)
                                    if job:
                                        # Mark job as started by this worker
                                        JobManager.start_job(
                                            redis_client, job_id, worker_id
                                        )

                                        # Emit queue update
                                        emit_queue_update(
                                            redis_client,
                                            room_id,
                                            category,
                                            extension,
                                            socketio,
                                        )

                                        # Get updated job details (now with running status)
                                        job = JobManager.get_job(redis_client, job_id)

                                        # Rename 'id' to 'jobId' for consistency with client expectations
                                        job["jobId"] = job.pop("id")

                                        log.info(
                                            f"Assigned celery job {job_id} to worker {worker_id} from queue"
                                        )
                                        return job, 200
            else:
                # Get all extensions this worker is registered for in this category
                user_extensions_key = ExtensionKeys.user_extensions_key(
                    room_id, category, worker_id
                )
                registered_extensions = redis_client.smembers(user_extensions_key)

                # Check each registered extension for queued jobs
                for extension in registered_extensions:
                    keys = ExtensionKeys.for_extension(room_id, category, extension)
                    queue_length = redis_client.llen(keys.queue)

                    if queue_length > 0:
                        # Get the first job from the queue
                        task_data = redis_client.lpop(keys.queue)
                        if task_data:
                            try:
                                task_info = json.loads(task_data)
                            except json.JSONDecodeError as e:
                                log.error(
                                    f"Invalid JSON in queue {keys.queue}: {e}, data: {task_data}"
                                )
                                # Continue to next extension, the malformed job was already popped
                                continue

                            job_id = task_info.get("jobId")

                            # Get full job details
                            job = JobManager.get_job(redis_client, job_id)
                            if job:
                                # Move worker from idle to progressing
                                redis_client.smove(
                                    keys.idle_workers,
                                    keys.progressing_workers,
                                    worker_id,
                                )

                                # Mark job as started by this worker
                                JobManager.start_job(redis_client, job_id, worker_id)

                                # Emit queue update (this covers job:started notification)
                                emit_queue_update(
                                    redis_client, room_id, category, extension, socketio
                                )

                                # Get updated job details (now with running status)
                                job = JobManager.get_job(redis_client, job_id)

                                # Rename 'id' to 'jobId' for consistency with client expectations
                                job["jobId"] = job.pop("id")

                                log.info(
                                    f"Assigned job {job_id} to worker {worker_id} from queue"
                                )
                                return job, 200

        # No jobs available
        return {"error": "No jobs available"}, 400

    except Exception as e:
        log.error(
            f"Unexpected error in get_next_job for worker {worker_id if 'worker_id' in locals() else 'unknown'}: {e}",
            exc_info=True,
        )
        return {"error": f"Internal server error: {str(e)}"}, 500


@main.route("/api/rooms/<string:room_id>/chat/messages", methods=["GET"])
def get_chat_messages(room_id: str):
    """
    Get paginated chat messages for a room.

    Query Parameters:
        - limit (int): Number of messages (default: 30, max: 100)
        - before (int): Get messages before this timestamp
        - after (int): Get messages after this timestamp

    Returns:
    {
        "messages": [Message],
        "metadata": {
            "hasMore": bool,
            "totalCount": int,
            "oldestTimestamp": int | null,
            "newestTimestamp": int | null
        }
    }
    """
    r = current_app.extensions["redis"]

    # Parse and validate query parameters
    try:
        limit = int(request.args.get("limit", 30))
        limit = max(1, min(limit, 100))  # Clamp between 1 and 100
    except (ValueError, TypeError):
        return {"error": "Invalid limit parameter"}, 400

    before = request.args.get("before")
    after = request.args.get("after")

    try:
        before = int(before) if before else None
        after = int(after) if after else None
    except (ValueError, TypeError):
        return {"error": "Invalid before/after parameter"}, 400

    index_key = f"room:{room_id}:chat:index"
    data_key = f"room:{room_id}:chat:data"

    # Get total count
    total_count = r.zcard(index_key)

    # Determine range query based on before/after
    if after is not None:
        # Get messages after timestamp (ascending order, then reverse)
        message_ids = r.zrangebyscore(
            index_key, f"({after}", "+inf", start=0, num=limit
        )
        # Reverse to maintain newest-first order
        message_ids = list(reversed(message_ids))
    elif before is not None:
        # Get messages before timestamp (descending order)
        message_ids = r.zrevrangebyscore(
            index_key, f"({before}", "-inf", start=0, num=limit
        )
    else:
        # Get latest messages (descending order)
        message_ids = r.zrevrangebyscore(index_key, "+inf", "-inf", start=0, num=limit)

    # Fetch message data
    messages = []
    if message_ids:
        with r.pipeline() as pipe:
            for msg_id in message_ids:
                pipe.hget(data_key, msg_id)
            message_data_list = pipe.execute()

        messages = [json.loads(msg_data) for msg_data in message_data_list if msg_data]

    # Calculate metadata
    has_more = False
    oldest_timestamp = None
    newest_timestamp = None

    if messages:
        oldest_timestamp = messages[-1]["createdAt"]
        newest_timestamp = messages[0]["createdAt"]

        # Check if there are more messages
        if after is not None:
            # Check if there are more recent messages
            count_after = r.zcount(index_key, f"({newest_timestamp}", "+inf")
            has_more = count_after > 0
        else:
            # Check if there are older messages
            count_before = r.zcount(index_key, "-inf", f"({oldest_timestamp}")
            has_more = count_before > 0

    return {
        "messages": messages,
        "metadata": {
            "hasMore": has_more,
            "totalCount": total_count,
            "oldestTimestamp": oldest_timestamp,
            "newestTimestamp": newest_timestamp,
        },
    }, 200


@main.route("/api/rooms/<string:room_id>/join", methods=["POST"])
def join_room(room_id):
    """Join a room (requires JWT authentication).

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
        "clientId": "uuid-string",
        "frameCount": 0,
        "step": 0,
        "created": true,
        ...
    }
    """
    from zndraw.auth import AuthError, get_current_client

    data = request.get_json() or {}
    if ":" in room_id:
        return {"error": "Room ID cannot contain ':' character"}, 400

    # Authenticate request (JWT required)
    try:
        client = get_current_client()
        client_id = client["clientId"]
        user_name = client["userName"]
    except AuthError as e:
        return {"error": e.message}, e.status_code
    description = data.get("description")
    copy_from = data.get("copyFrom")
    allow_create = data.get("allowCreate", True)
    r = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    client_service = current_app.extensions["client_service"]
    settings_service = current_app.extensions["settings_service"]

    # Check if room already exists
    room_exists = room_service.room_exists(room_id)

    # If allowCreate is False and room doesn't exist, return 404
    if not allow_create and not room_exists:
        return {
            "status": "not_found",
            "message": f"Room '{room_id}' does not exist yet. It may still be loading.",
        }, 404

    # Update client metadata (userName is not room-specific, so handled separately)
    client_key = f"client:{client_id}"
    r.hset(client_key, "userName", user_name)

    # Update client room membership atomically
    client_service.update_client_and_room_membership(client_id, room_id)

    log.info(f"Client {client_id} ({user_name}) joined room: {room_id}")

    response = {
        "status": "ok",
        "clientId": client_id,
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
            result = room_service.create_room(room_id, user_name, description, copy_from)
            frame_count = result["frameCount"]
        except ValueError as e:
            return {"error": str(e)}, 404

        # Initialize default settings for new user in room
        settings_service.initialize_defaults(room_id, user_name)

        # Broadcast room creation to all connected clients

        # Check if metadata lock is held (trajectory:meta is the target used by vis.lock)
        metadata_lock_key = get_lock_key(room_id, "trajectory:meta")
        metadata_locked = None
        if r.exists(metadata_lock_key):
            # Lock exists - get metadata
            metadata_raw = r.get(f"{metadata_lock_key}:metadata")
            if metadata_raw:
                metadata = json.loads(metadata_raw)
                lock_metadata = LockMetadata(
                    msg=metadata.get("msg"),
                    userName=metadata.get("userName"),
                    timestamp=metadata.get("timestamp"),
                )
                metadata_locked = lock_metadata.model_dump()
            else:
                # Lock exists but no metadata - get basic info from client
                lock_holder_client_id = r.get(metadata_lock_key)
                user_name = (
                    r.hget(f"client:{lock_holder_client_id}", "userName")
                    if lock_holder_client_id
                    else None
                )
                lock_metadata = LockMetadata(
                    msg=None, userName=user_name or "unknown", timestamp=None
                )
                metadata_locked = lock_metadata.model_dump()

        from zndraw.app.room_manager import emit_room_update

        emit_room_update(
            socketio,
            room_id,
            created=True,
            description=description,
            frameCount=frame_count,
            locked=False,
            hidden=False,
            isDefault=False,
            metadataLocked=metadata_locked,
        )

    # Get frame count using service
    response["frameCount"] = room_service.get_frame_count(room_id)

    selections_raw = r.hgetall(f"room:{room_id}:selections")
    selections = {k: json.loads(v) for k, v in selections_raw.items()}
    response["selections"] = selections

    frame_selection = r.get(f"room:{room_id}:frame_selection:default")
    response["frame_selection"] = (
        json.loads(frame_selection) if frame_selection else None
    )

    presenter_lock = r.get(f"room:{room_id}:presenter_lock")
    response["presenter-lock"] = presenter_lock

    # Get current frame using service (handles validation and error cases)
    response["step"] = room_service.get_current_frame(room_id)

    bookmarks_key = f"room:{room_id}:bookmarks"
    bookmarks_raw = r.hgetall(bookmarks_key)

    geometries = r.hgetall(f"room:{room_id}:geometries")
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
    groups_raw = r.hgetall(f"room:{room_id}:selection_groups")
    selection_groups = {}
    for group_name, group_data in groups_raw.items():
        selection_groups[group_name] = json.loads(group_data)
    response["selectionGroups"] = selection_groups

    # Fetch active selection group
    active_group = r.get(f"room:{room_id}:active_selection_group")
    response["activeSelectionGroup"] = active_group if active_group else None

    # Check if metadata lock is held (trajectory:meta is the target used by vis.lock)
    metadata_lock_key = get_lock_key(room_id, "trajectory:meta")
    if r.exists(metadata_lock_key):
        # Lock exists - get metadata
        metadata_raw = r.get(f"{metadata_lock_key}:metadata")
        if metadata_raw:
            metadata = json.loads(metadata_raw)
            lock_metadata = LockMetadata(
                msg=metadata.get("msg"),
                userName=metadata.get("userName"),
                timestamp=metadata.get("timestamp"),
            )
            response["metadataLocked"] = lock_metadata.model_dump()
        else:
            # Lock exists but no metadata (shouldn't happen, but handle gracefully)
            response["metadataLocked"] = None
    else:
        response["metadataLocked"] = None

    return response


@main.route("/api/rooms/<string:room_id>/geometries", methods=["POST"])
def create_geometry(room_id: str):
    """Create or update a geometry in the room.

    Request body:
        {
            "key": "geometry_name",
            "type": "Sphere" | "Arrow" | "Bond" | "Curve",
            "data": {...}  // geometry-specific data
        }
    """
    from zndraw.auth import AuthError, get_current_client

    # Authenticate and get client ID from JWT token
    try:
        client = get_current_client()
        client_id = client["clientId"]
    except AuthError as e:
        return {"error": e.message}, e.status_code

    data = request.get_json() or {}
    key = data.get("key")
    geometry_type = data.get("type")
    geometry_data = data.get("data")

    if not key or not geometry_type or geometry_data is None:
        return {
            "error": "'key', 'type', and 'data' are required",
            "type": "ValueError",
        }, 400

    from zndraw.geometries import geometries

    if geometry_type not in geometries:
        return {
            "error": f"Unknown geometry type '{geometry_type}'",
            "type": "ValueError",
        }, 400

    r = current_app.extensions["redis"]
    existing_geometry_json = r.hget(f"room:{room_id}:geometries", key)

    if existing_geometry_json:
        existing_geometry = json.loads(existing_geometry_json)
        # Make sure we're updating the same type of geometry
        if existing_geometry.get("type") == geometry_type:
            # Merge new data into existing data
            merged_data = existing_geometry.get("data", {}).copy()
            merged_data.update(geometry_data)
            geometry_data = merged_data
        else:
            log.warning(f"Geometry type mismatch for key '{key}'. Overwriting.")

    # Validate and apply defaults through Pydantic model
    try:
        geometry_class = geometries[geometry_type]
        validated_geometry = geometry_class(**geometry_data)
        value_to_store = json.dumps(
            {"type": geometry_type, "data": validated_geometry.model_dump()}
        )
    except Exception as e:
        return {
            "error": f"Invalid geometry data: {str(e)}",
            "type": "ValidationError",
        }, 400

    r.hset(f"room:{room_id}:geometries", key, value_to_store)
    socketio.emit(
        SocketEvents.INVALIDATE_GEOMETRY,
        {
            "key": key,
            "operation": "set",
        },
        to=f"room:{room_id}",
    )

    return {"status": "success"}, 200


@main.route("/api/rooms/<string:room_id>/geometries/<string:key>", methods=["GET"])
def get_geometry(room_id: str, key: str):
    """Get a specific geometry by key.

    Returns:
        {
            "key": "geometry_name",
            "geometry": {
                "type": "Sphere",
                "data": {...}
            }
        }
    """
    r = current_app.extensions["redis"]
    geometry_data = r.hget(f"room:{room_id}:geometries", key)
    if not geometry_data:
        return {
            "error": f"Geometry with key '{key}' not found",
            "type": "KeyError",
        }, 404
    geometry = json.loads(geometry_data)
    return {"key": key, "geometry": geometry}, 200


@main.route("/api/rooms/<string:room_id>/geometries/<string:key>", methods=["DELETE"])
def delete_geometry(room_id: str, key: str):
    from zndraw.auth import AuthError, get_current_client

    # Authenticate and get client ID from JWT token
    try:
        client = get_current_client()
        client_id = client["clientId"]
    except AuthError as e:
        return {"error": e.message}, e.status_code

    r = current_app.extensions["redis"]
    response = r.hdel(f"room:{room_id}:geometries", key)
    if response == 0:
        return {
            "error": f"Geometry with key '{key}' does not exist",
            "type": "KeyError",
        }, 404

    socketio.emit(
        SocketEvents.INVALIDATE_GEOMETRY,
        {
            "key": key,
            "operation": "delete",
        },
        to=f"room:{room_id}",
    )
    return {"status": "success"}, 200


@main.route("/api/rooms/<string:room_id>/geometries", methods=["GET"])
def list_geometries(room_id: str):
    """List all geometry keys in the room.

    Returns:
        {"geometries": ["key1", "key2", ...]}
    """
    r = current_app.extensions["redis"]
    all_keys = r.hkeys(f"room:{room_id}:geometries")
    return {"geometries": list(all_keys)}, 200


@main.route("/api/rooms/<string:room_id>/geometries/schemas", methods=["GET"])
def list_geometry_schemas(room_id: str):
    """Return JSON schemas for all geometry types for form generation."""
    from zndraw.geometries import geometries

    schemas = {name: model.model_json_schema() for name, model in geometries.items()}
    return {"schemas": schemas}, 200


# -------------#
### FIGURES ###
# -------------#


@main.route("/api/rooms/<string:room_id>/figures", methods=["POST"])
def create_figure(room_id: str):
    data = request.get_json() or {}
    key = data.get("key")
    figure = data.get("figure")

    if not key or not figure:
        return {
            "error": "Both 'key' and 'figure' are required",
            "type": "ValueError",
        }, 400

    # store in hash
    r = current_app.extensions["redis"]
    r.hset(f"room:{room_id}:figures", key, json.dumps(figure))
    socketio.emit(
        SocketEvents.INVALIDATE_FIGURE,
        {
            "key": key,
            "operation": "set",
        },
        to=f"room:{room_id}",
    )
    return {"status": "success"}, 200


@main.route("/api/rooms/<string:room_id>/figures/<string:key>", methods=["GET"])
def get_figure(room_id: str, key: str):
    r = current_app.extensions["redis"]
    figure_data = r.hget(f"room:{room_id}:figures", key)
    if not figure_data:
        return {"error": f"Figure with key '{key}' not found", "type": "KeyError"}, 404
    figure = json.loads(figure_data)
    return {"key": key, "figure": figure}, 200


@main.route("/api/rooms/<string:room_id>/figures/<string:key>", methods=["DELETE"])
def delete_figure(room_id: str, key: str):
    r = current_app.extensions["redis"]
    response = r.hdel(f"room:{room_id}:figures", key)
    if response == 0:
        return {
            "error": f"Figure with key '{key}' does not exist",
            "type": "KeyError",
        }, 404
    socketio.emit(
        SocketEvents.INVALIDATE_FIGURE,
        {
            "key": key,
            "operation": "delete",
        },
        to=f"room:{room_id}",
    )
    return {"status": "success"}, 200


@main.route("/api/rooms/<string:room_id>/figures", methods=["GET"])
def list_figures(room_id: str):
    r = current_app.extensions["redis"]
    all_keys = r.hkeys(f"room:{room_id}:figures")
    return {"figures": list(all_keys)}, 200


# ============================================================================
# Bookmarks API Routes
# ============================================================================


@main.route("/api/rooms/<string:room_id>/bookmarks", methods=["GET"])
def get_all_bookmarks(room_id: str):
    """Get all bookmarks for a room.

    Returns:
        {"bookmarks": {1: "First Frame", 5: "Middle Frame"}}
    """
    r = current_app.extensions["redis"]
    bookmarks_raw = r.hgetall(f"room:{room_id}:bookmarks")
    # Convert byte keys to integers
    bookmarks = {int(k): v for k, v in bookmarks_raw.items()}
    return {"bookmarks": bookmarks}, 200


@main.route("/api/rooms/<string:room_id>/bookmarks/<int:index>", methods=["GET"])
def get_bookmark(room_id: str, index: int):
    """Get a specific bookmark by frame index."""
    r = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]

    # Check if frame index is valid
    frame_count = room_service.get_frame_count(room_id)
    if index < 0 or index >= frame_count:
        return {
            "error": f"Bookmark index {index} out of range (0-{frame_count - 1})",
            "type": "IndexError",
        }, 404

    label = r.hget(f"room:{room_id}:bookmarks", str(index))
    if label is None:
        return {
            "error": f"Bookmark at index {index} does not exist",
            "type": "KeyError",
        }, 404

    return {"index": index, "label": label}, 200


@main.route("/api/rooms/<string:room_id>/bookmarks/<int:index>", methods=["PUT"])
def set_bookmark(room_id: str, index: int):
    """Set or update a bookmark at a specific frame index.

    Body: {"label": "Frame Label"}
    """
    r = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    data = request.get_json() or {}
    label = data.get("label")

    if not label or not isinstance(label, str):
        return {
            "error": "Bookmark label must be a non-empty string",
            "type": "ValueError",
        }, 400

    # Check if frame index is valid
    frame_count = room_service.get_frame_count(room_id)
    if index < 0 or index >= frame_count:
        return {
            "error": f"Bookmark index {index} out of range (0-{frame_count - 1})",
            "type": "IndexError",
        }, 400

    # Set the bookmark
    r.hset(f"room:{room_id}:bookmarks", str(index), label)

    # Emit invalidate event
    socketio.emit(
        SocketEvents.INVALIDATE_BOOKMARK,
        {"index": index, "operation": "set"},
        to=f"room:{room_id}",
    )

    return {"status": "success"}, 200


@main.route("/api/rooms/<string:room_id>/bookmarks/<int:index>", methods=["DELETE"])
def delete_bookmark(room_id: str, index: int):
    """Delete a bookmark at a specific frame index."""
    r = current_app.extensions["redis"]

    response = r.hdel(f"room:{room_id}:bookmarks", str(index))
    if response == 0:
        return {
            "error": f"Bookmark at index {index} does not exist",
            "type": "KeyError",
        }, 404

    # Emit invalidate event
    socketio.emit(
        SocketEvents.INVALIDATE_BOOKMARK,
        {"index": index, "operation": "delete"},
        to=f"room:{room_id}",
    )

    return {"status": "success"}, 200


# ============================================================================
# Selection API Routes
# ============================================================================


@main.route("/api/rooms/<string:room_id>/selections", methods=["GET"])
def get_all_selections(room_id: str):
    """Get all current selections and groups.

    Returns:
        {
            "selections": {"particles": [1,2,3], "forces": [2,3]},
            "groups": {"group1": {"particles": [1,3], "forces": [1,3]}},
            "activeGroup": "group1" | null
        }
    """
    r = current_app.extensions["redis"]

    # Get current selections
    selections_raw = r.hgetall(f"room:{room_id}:selections")
    selections = {k: json.loads(v) for k, v in selections_raw.items()}

    # Get selection groups
    groups_raw = r.hgetall(f"room:{room_id}:selection_groups")
    groups = {k: json.loads(v) for k, v in groups_raw.items()}

    # Get active group
    active_group = r.get(f"room:{room_id}:active_selection_group")

    return {
        "selections": selections,
        "groups": groups,
        "activeGroup": active_group,
    }, 200


@main.route("/api/rooms/<string:room_id>/selections/<string:geometry>", methods=["GET"])
def get_selection(room_id: str, geometry: str):
    """Get selection for a specific geometry."""
    r = current_app.extensions["redis"]
    selection = r.hget(f"room:{room_id}:selections", geometry)

    if selection is None:
        return {"selection": []}, 200

    return {"selection": json.loads(selection)}, 200


@main.route("/api/rooms/<string:room_id>/selections/<string:geometry>", methods=["PUT"])
def update_selection(room_id: str, geometry: str):
    """Update selection for a specific geometry.

    Body: {"indices": [1, 2, 3]}
    """
    r = current_app.extensions["redis"]
    data = request.get_json()

    indices = data.get("indices", [])
    if not isinstance(indices, list):
        return {"error": "indices must be a list"}, 400

    if any(not isinstance(idx, int) or idx < 0 for idx in indices):
        return {"error": "All indices must be non-negative integers"}, 400

    # Store selection
    r.hset(f"room:{room_id}:selections", geometry, json.dumps(indices))

    # Clear active group (manual edit breaks group association)
    r.delete(f"room:{room_id}:active_selection_group")

    # Emit invalidation
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"geometry": geometry},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200


@main.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>",
    methods=["GET"],
)
def get_selection_group(room_id: str, group_name: str):
    """Get a specific selection group."""
    r = current_app.extensions["redis"]
    group = r.hget(f"room:{room_id}:selection_groups", group_name)

    if group is None:
        return {"error": "Group not found"}, 404

    return {"group": json.loads(group)}, 200


@main.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>",
    methods=["PUT"],
)
def create_update_selection_group(room_id: str, group_name: str):
    """Create or update a selection group.

    Body: {"particles": [1, 3], "forces": [1, 3]}
    """
    r = current_app.extensions["redis"]
    data = request.get_json()

    # Validate data is a dict of geometry -> indices
    if not isinstance(data, dict):
        return {"error": "Group must be a dictionary"}, 400

    for geometry, indices in data.items():
        if not isinstance(indices, list):
            return {"error": f"Indices for '{geometry}' must be a list"}, 400
        if any(not isinstance(idx, int) or idx < 0 for idx in indices):
            return {"error": f"Invalid indices for '{geometry}'"}, 400

    # Store group
    r.hset(f"room:{room_id}:selection_groups", group_name, json.dumps(data))

    # Emit invalidation (groups list changed)
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION_GROUPS,
        {"operation": "group_saved", "group": group_name},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200


@main.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>",
    methods=["DELETE"],
)
def delete_selection_group(room_id: str, group_name: str):
    """Delete a selection group."""
    r = current_app.extensions["redis"]

    # Check if group exists
    if not r.hexists(f"room:{room_id}:selection_groups", group_name):
        return {"error": "Group not found"}, 404

    # Delete group
    r.hdel(f"room:{room_id}:selection_groups", group_name)

    # Clear active group if it was this one
    active_group = r.get(f"room:{room_id}:active_selection_group")
    if active_group == group_name:
        r.delete(f"room:{room_id}:active_selection_group")

    # Emit invalidation
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION_GROUPS,
        {"operation": "group_deleted", "group": group_name},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200


@main.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>/load",
    methods=["POST"],
)
def load_selection_group(room_id: str, group_name: str):
    """Load a selection group (apply it to current selections).

    This sets the active group and updates all selections to match the group.
    """
    r = current_app.extensions["redis"]

    # Get group
    group_data = r.hget(f"room:{room_id}:selection_groups", group_name)
    if group_data is None:
        return {"error": "Group not found"}, 404

    group = json.loads(group_data)

    # Apply group to current selections
    for geometry, indices in group.items():
        r.hset(f"room:{room_id}:selections", geometry, json.dumps(indices))

    # Set as active group
    r.set(f"room:{room_id}:active_selection_group", group_name)

    # Emit invalidation (all selections changed)
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"operation": "group_loaded", "group": group_name},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200


# ============================================================================
# End Selection API Routes
# ============================================================================


# ============================================================================
# Screenshot API Routes
# ============================================================================

from zndraw.screenshot_manager import ScreenshotManager

# 10MB max file size for screenshots
MAX_SCREENSHOT_SIZE = 10 * 1024 * 1024


def _get_screenshot_manager(room_id: str) -> ScreenshotManager:
    """Helper to create ScreenshotManager instance."""
    storage_path = current_app.config.get("STORAGE_PATH", "./zndraw-data.zarr")
    return ScreenshotManager(room_id, storage_path)


@main.route("/api/rooms/<string:room_id>/screenshots/upload", methods=["POST"])
def upload_screenshot(room_id: str):
    """Upload screenshot from frontend.

    Request:
        multipart/form-data with:
        - file: image file
        - format: png/jpeg/webp
        - width: optional image width
        - height: optional image height

    Response:
        JSON with screenshot metadata
    """
    if "file" not in request.files:
        return {"error": "No file provided"}, 400

    file = request.files["file"]
    if file.filename == "":
        return {"error": "Empty filename"}, 400

    format = request.form.get("format", "png")
    width = request.form.get("width", type=int)
    height = request.form.get("height", type=int)

    try:
        image_data = file.read()

        # Validate file size
        if len(image_data) > MAX_SCREENSHOT_SIZE:
            return {
                "error": f"File too large. Maximum size is {MAX_SCREENSHOT_SIZE // 1024 // 1024}MB"
            }, 400

        manager = _get_screenshot_manager(room_id)
        screenshot = manager.save(image_data, format, width, height)

        socketio.emit(
            "screenshot:created",
            {"id": screenshot.id},
            to=f"room:{room_id}",
        )

        return {
            "id": screenshot.id,
            "format": screenshot.format,
            "size": screenshot.size,
            "url": f"/api/rooms/{room_id}/screenshots/{screenshot.id}",
        }, 201
    except ValueError as e:
        return {"error": str(e)}, 400
    except Exception as e:
        return {"error": f"Failed to save screenshot: {str(e)}"}, 500


@main.route("/api/rooms/<string:room_id>/screenshots", methods=["GET"])
def list_screenshots(room_id: str):
    """List all screenshots for a room.

    Query params:
        limit: max results (default 20)
        offset: skip N results (default 0)

    Response:
        JSON with screenshots array and total count
    """
    limit = request.args.get("limit", 20, type=int)
    offset = request.args.get("offset", 0, type=int)

    if limit < 1 or limit > 100:
        return {"error": "Limit must be between 1 and 100"}, 400
    if offset < 0:
        return {"error": "Offset must be non-negative"}, 400

    try:
        manager = _get_screenshot_manager(room_id)
        screenshots = manager.list(limit, offset)
        total = manager.count()

        return {
            "screenshots": [
                {
                    "id": s.id,
                    "format": s.format,
                    "size": s.size,
                    "width": s.width,
                    "height": s.height,
                    "url": f"/api/rooms/{room_id}/screenshots/{s.id}",
                }
                for s in screenshots
            ],
            "total": total,
            "limit": limit,
            "offset": offset,
        }, 200
    except Exception as e:
        return {"error": f"Failed to list screenshots: {str(e)}"}, 500


@main.route(
    "/api/rooms/<string:room_id>/screenshots/<int:screenshot_id>", methods=["GET"]
)
def get_screenshot(room_id: str, screenshot_id: int):
    """Download a specific screenshot.

    Returns the image file with appropriate Content-Type header.
    """
    try:
        manager = _get_screenshot_manager(room_id)
        result = manager.get(screenshot_id)

        if not result:
            return {"error": "Screenshot not found"}, 404

        filepath, metadata = result
        return send_from_directory(
            filepath.parent,
            filepath.name,
            mimetype=f"image/{metadata.format}",
        )
    except Exception as e:
        return {"error": f"Failed to get screenshot: {str(e)}"}, 500


@main.route(
    "/api/rooms/<string:room_id>/screenshots/<int:screenshot_id>/metadata",
    methods=["GET"],
)
def get_screenshot_metadata(room_id: str, screenshot_id: int):
    """Get metadata for a specific screenshot."""
    try:
        manager = _get_screenshot_manager(room_id)
        result = manager.get(screenshot_id)

        if not result:
            return {"error": "Screenshot not found"}, 404

        _, metadata = result
        return {
            "id": metadata.id,
            "format": metadata.format,
            "size": metadata.size,
            "width": metadata.width,
            "height": metadata.height,
            "url": f"/api/rooms/{room_id}/screenshots/{metadata.id}",
        }, 200
    except Exception as e:
        return {"error": f"Failed to get screenshot metadata: {str(e)}"}, 500


@main.route(
    "/api/rooms/<string:room_id>/screenshots/<int:screenshot_id>", methods=["DELETE"]
)
def delete_screenshot(room_id: str, screenshot_id: int):
    """Delete a screenshot."""
    try:
        manager = _get_screenshot_manager(room_id)

        if manager.delete(screenshot_id):
            socketio.emit(
                "screenshot:deleted",
                {"id": screenshot_id},
                to=f"room:{room_id}",
            )
            return {"success": True}, 200

        return {"error": "Screenshot not found"}, 404
    except Exception as e:
        return {"error": f"Failed to delete screenshot: {str(e)}"}, 500


# ============================================================================
# End Screenshot API Routes
# ============================================================================


def _transition_worker_to_idle(
    redis_client,
    socketio_instance,
    worker_id: str,
    job: dict,
    room_id: str,
    success: bool = True,
) -> None:
    """Transition worker from progressing to idle and dispatch next task.

    Args:
        redis_client: Redis client
        socketio_instance: SocketIO instance
        worker_id: Worker session ID
        job: Job data dict containing category/extension
        room_id: Room identifier
        success: True if job completed, False if failed (affects log message)
    """
    category = job.get("category")
    extension = job.get("extension")

    log.info(
        f"_transition_worker_to_idle called: worker_id={worker_id}, category={category}, extension={extension}"
    )

    if not category or not extension:
        log.warning(
            f"Missing category or extension in job data: category={category}, extension={extension}"
        )
        return

    # Move worker from progressing back to idle for the extension they just completed
    keys = ExtensionKeys.for_extension(room_id, category, extension)
    moved = redis_client.smove(keys.progressing_workers, keys.idle_workers, worker_id)

    log.info(
        f"Worker transition: moved={moved}, progressing_key={keys.progressing_workers}, idle_key={keys.idle_workers}"
    )

    if moved:
        status_msg = (
            "finished and is now idle" if success else "marked idle after failure"
        )
        log.info(f"Worker {worker_id} {status_msg}")

        # Get all extensions this worker is registered for in this category
        user_extensions_key = ExtensionKeys.user_extensions_key(
            room_id, category, worker_id
        )
        registered_extensions = redis_client.smembers(user_extensions_key)

        log.info(f"Worker {worker_id} registered extensions: {registered_extensions}")

        # Add worker to idle set for ALL registered extensions (not just the one they completed)
        # This ensures dispatch_next_task can find them for any extension
        for ext_name in registered_extensions:
            if ext_name != extension:  # Already moved above for completed extension
                ext_keys = ExtensionKeys.for_extension(room_id, category, ext_name)
                # Remove from progressing (if present) and add to idle
                redis_client.srem(ext_keys.progressing_workers, worker_id)
                redis_client.sadd(ext_keys.idle_workers, worker_id)
                log.info(
                    f"Added worker {worker_id} to idle set for extension {ext_name}"
                )

        # Emit queue update for this extension
        emit_queue_update(redis_client, room_id, category, extension, socketio_instance)

        # Workers will poll for the next task via /jobs/next endpoint
        log.info(f"Worker {worker_id} is now idle and can poll for next task")
    else:
        log.warning(
            f"Failed to move worker {worker_id} from progressing to idle (may already be idle or not in progressing)"
        )


@main.route("/api/rooms/<string:room_id>/download", methods=["GET"])
def download_frames(room_id: str):
    """Download frames in ExtendedXYZ format using streaming.

    Query Parameters
    ----------------
    indices : str, optional
        Comma-separated frame indices (e.g., '0,5,10'). If not provided, downloads all frames.
    selection : str, optional
        Comma-separated particle indices to filter.
    filename : str, optional
        Custom filename for download.

    Returns
    -------
    Response
        Streaming XYZ file download or error.
    """
    from io import StringIO

    import ase.io

    from zndraw.utils import atoms_from_dict

    redis_client = current_app.extensions["redis"]

    # Get parameters
    indices_param = request.args.get("indices")
    selection_param = request.args.get("selection")
    custom_filename = request.args.get("filename")

    # Parse selection if provided
    selection = None
    if selection_param:
        try:
            selection = [int(i.strip()) for i in selection_param.split(",")]
        except ValueError:
            return {"error": "Invalid selection format"}, 400

    # Get frame indices manager with correct Redis key
    indices_key = f"room:{room_id}:trajectory:indices"
    fim = FrameIndexManager(redis_client, indices_key)

    # Check if room has frames
    frame_count = len(fim)
    if frame_count == 0:
        return {"error": "No frames available in room"}, 400

    # Determine which frames to export
    frame_indices = []

    if indices_param:
        # Use explicit indices parameter
        try:
            frame_indices = [int(i.strip()) for i in indices_param.split(",")]
        except ValueError:
            return {"error": "Invalid indices format"}, 400
    else:
        # No indices provided = download all frames
        frame_indices = list(range(frame_count))

    # Validate frame indices
    for idx in frame_indices:
        if idx < 0 or idx >= frame_count:
            return {
                "error": f"Frame index {idx} out of range (0-{frame_count - 1})"
            }, 400

    # Generate filename
    if not custom_filename:
        if len(frame_indices) == 1:
            custom_filename = f"{room_id}_frame_{frame_indices[0]}.xyz"
        else:
            custom_filename = f"{room_id}_{len(frame_indices)}_frames.xyz"

    def generate_xyz_chunks():
        """Generator that yields XYZ format chunks for streaming."""
        buffer = StringIO()

        # Get physical keys for all frame indices at once (efficient)
        physical_keys = fim.get_by_indices(frame_indices)

        for idx, mapping_entry in zip(frame_indices, physical_keys):
            try:
                # Decode bytes to string if necessary
                mapping_entry_str = (
                    mapping_entry.decode("utf-8")
                    if isinstance(mapping_entry, bytes)
                    else mapping_entry
                )

                # Parse the mapping entry to get source room and physical index
                if ":" in mapping_entry_str:
                    source_room_id, physical_index_str = mapping_entry_str.split(":", 1)
                    physical_index = int(physical_index_str)
                    source_storage = get_storage(source_room_id)
                else:
                    physical_index = int(mapping_entry_str)
                    source_storage = get_storage(room_id)

                # Access storage using physical index
                atoms_dict = source_storage[physical_index]

                # Convert to ASE Atoms object using the proper utility function
                atoms = atoms_from_dict(atoms_dict)

                # Apply selection filter if provided
                if selection:
                    atoms = atoms[selection]

                # Write single frame to buffer in ExtendedXYZ format
                ase.io.write(buffer, atoms, format="extxyz", write_info=True)

                # Get string content and encode to bytes
                chunk = buffer.getvalue().encode("utf-8")

                # Yield the chunk
                yield chunk

                # Clear buffer for next frame
                buffer.seek(0)
                buffer.truncate(0)

            except Exception as e:
                log.error(f"Error streaming frame {idx}: {e}")
                # Continue with next frame instead of failing entire download
                continue

    # Return streaming response
    return Response(
        generate_xyz_chunks(),
        mimetype="chemical/x-xyz",
        headers={
            "Content-Disposition": f'attachment; filename="{custom_filename}"',
            "Content-Type": "chemical/x-xyz",
        },
    )
