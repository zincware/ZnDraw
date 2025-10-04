import json
import logging
import traceback

import msgpack
import zarr
from flask import Response, current_app, request
from flask_socketio import disconnect
from zarr.storage import MemoryStore

from zndraw.server import socketio
from zndraw.storage import ZarrStorageSequence, decode_data, encode_data

from . import main
from .constants import SocketEvents
from .job_manager import JobManager
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


def emit_bookmarks_update(room_id: str):
    """
    Emit bookmarks update to all clients after frame mapping changes.
    Converts physical bookmarks to logical indices based on current mapping.
    """
    r = current_app.extensions["redis"]
    indices_key = f"room:{room_id}:trajectory:indices"
    bookmarks_key = f"room:{room_id}:bookmarks"

    physical_bookmarks = r.hgetall(bookmarks_key)

    if physical_bookmarks:
        # Get current frame mapping
        frame_mapping = r.zrange(indices_key, 0, -1)

        # Build reverse mapping: physical_key -> logical_index
        physical_to_logical = {
            physical_key: idx for idx, physical_key in enumerate(frame_mapping)
        }

        # Convert bookmarks from physical to logical indices
        logical_bookmarks = {}
        for physical_key, label in physical_bookmarks.items():
            if physical_key in physical_to_logical:
                logical_idx = physical_to_logical[physical_key]
                logical_bookmarks[str(logical_idx)] = label

        # Emit bookmarks update to all clients
        socketio.emit(
            "bookmarks:update",
            {"bookmarks": logical_bookmarks},
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
    """Serves multiple frames' data from the room's Zarr store using either indices or slice parameters."""
    r = current_app.extensions["redis"]
    try:
        storage = get_storage(room_id)

        # Get logical-to-physical mapping from Redis
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {
                "error": f"Index out of range for data with 0 frames in room '{room_id}'",
                "type": "IndexError",
            }, 404

        max_frame = len(frame_mapping) - 1

        # Determine frame indices based on request parameters
        if "indices" in request.args:
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
            # Default to slice behavior for any remaining cases (including empty payload)
            # This handles slice parameters and slice(None, None, None) which sends empty payload
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

        # Get keys parameter if specified (comma-separated)
        keys_str = request.args.get("keys")
        if keys_str:
            requested_keys = [k.strip() for k in keys_str.split(",")]
        else:
            requested_keys = None

        # TODO: requested keys and KeyError handling
        # TODO: instead of iterate load at once
        try:
            frames_data = []
            for frame_id in frame_indices:
                # Parse the frame mapping entry (format: "room_id:physical_index")
                mapping_entry = frame_mapping[frame_id]

                if ":" in mapping_entry:
                    source_room_id, physical_index_str = mapping_entry.split(":", 1)
                    physical_index = int(physical_index_str)
                    source_storage = get_storage(source_room_id)
                else:
                    # Fallback for any legacy data without colon separator
                    physical_index = int(mapping_entry)
                    source_storage = storage

                frame_data = source_storage.get(physical_index, keys=requested_keys)
                frames_data.append(encode_data(frame_data))
        except KeyError as e:
            error_data = {
                "error": f"Key(s) not found: {e}",
                "type": "KeyError",
            }
            return Response(
                json.dumps(error_data), status=404, content_type="application/json"
            )
        except IndexError as e:
            error_data = {
                "error": f"Frame index out of range: {e}",
                "type": "IndexError",
            }
            return Response(
                json.dumps(error_data), status=404, content_type="application/json"
            )

        packed_data = msgpack.packb(frames_data)
        return Response(packed_data, content_type="application/octet-stream")
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

        if ":" in mapping_entry:
            source_room_id, physical_index_str = mapping_entry.split(":", 1)
            physical_index = int(physical_index_str)
            source_storage = get_storage(source_room_id)
        else:
            # Fallback for any legacy data without colon separator
            physical_index = int(mapping_entry)
            source_storage = storage

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
    r = current_app.extensions["redis"]

    try:
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

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

        # Optimized batch delete: calculate final state once, then apply in single transaction
        # Convert to set for O(1) lookup
        indices_to_delete = set(frame_indices)

        # Build the list of frames that should remain after deletion
        remaining_frames = []
        for idx, physical_idx_str in enumerate(frame_mapping):
            if idx not in indices_to_delete:
                remaining_frames.append(physical_idx_str)

        # Use a pipeline to delete and rebuild the mapping in a single atomic transaction
        pipeline = r.pipeline()
        pipeline.delete(indices_key)

        if remaining_frames:
            # Build new mapping with sequential logical indices
            new_mapping = {
                physical_idx_str: logical_pos
                for logical_pos, physical_idx_str in enumerate(remaining_frames)
            }
            pipeline.zadd(indices_key, new_mapping)

        pipeline.execute()

        log.info(
            f"Deleted {len(frame_indices)} frames from room '{room_id}'. Physical data preserved."
        )

        # Emit bookmarks update to all clients
        emit_bookmarks_update(room_id)
        # Invalidate all frames from the first deleted position onward
        if frame_indices:
            emit_frames_invalidate(
                room_id, operation="delete", affected_from=min(frame_indices)
            )

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

            pipeline = r.pipeline()
            pipeline.zrem(indices_key, old_mapping_entry)
            pipeline.zadd(
                indices_key, {f"{room_id}:{new_physical_index}": target_frame_id}
            )
            pipeline.execute()

            # Remove bookmark for old physical frame (if it exists)
            bookmarks_key = f"room:{room_id}:bookmarks"
            r.hdel(bookmarks_key, old_mapping_entry)

            # Emit bookmarks update to reflect removal
            emit_bookmarks_update(room_id)
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
            start_logical_pos = r.zcard(indices_key)
            start_physical_pos = len(storage)
            num_frames = len(serialized_data)

            # 2. Decode all frames and extend the physical storage
            decoded_frames = [decode_data(frame) for frame in serialized_data]
            storage.extend(decoded_frames)

            # 3. Create the new logical-to-physical mapping for all new frames with room_id prefix
            new_mapping = {
                f"{room_id}:{start_physical_pos + i}": start_logical_pos + i
                for i in range(num_frames)
            }
            if new_mapping:
                r.zadd(indices_key, new_mapping)

            # 4. Prepare response data
            new_indices = list(range(start_logical_pos, start_logical_pos + num_frames))

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

            current_length = r.zcard(indices_key)

            if not (0 <= insert_position <= current_length):
                return {
                    "error": f"Insert position {insert_position} out of range [0, {current_length}]"
                }, 400

            # 1. Determine the new physical index (always at the end of the physical store)
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data to the physical store
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Atomically shift existing logical indices and insert the new one
            pipeline = r.pipeline()
            # Get members whose scores (logical positions) need to be incremented
            members_to_shift = r.zrangebyscore(indices_key, insert_position, "+inf")
            if members_to_shift:
                for member in members_to_shift:
                    pipeline.zincrby(indices_key, 1, member)

            # Add the new frame at the correct logical position with room_id prefix
            pipeline.zadd(
                indices_key, {f"{room_id}:{new_physical_index}": insert_position}
            )
            pipeline.execute()

            # Emit bookmarks update (logical indices shifted by the insert)
            emit_bookmarks_update(room_id)
            # Invalidate all frames from insert position onward (they all shift up)
            emit_frames_invalidate(
                room_id, operation="insert", affected_from=insert_position
            )

            log.info(
                f"Inserted frame at position {insert_position} (physical: {new_physical_index}) in room '{room_id}'"
            )
            return {"success": True, "inserted_position": insert_position}

        elif action == "append":
            # Append operation: add a new frame to the end of the logical sequence
            # 1. Determine the new logical and physical positions
            logical_position = r.zcard(indices_key)
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Add the new physical index to the logical sequence with room_id prefix
            r.zadd(indices_key, {f"{room_id}:{new_physical_index}": logical_position})

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
        frame_mapping = r.zrange(indices_key, 0, -1)

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
                new_physical_idx = start_physical_pos + i

                # Remove old mapping and add new one
                pipeline.zrem(indices_key, old_mapping_entry)
                pipeline.zadd(
                    indices_key, {f"{room_id}:{new_physical_idx}": logical_idx}
                )

            pipeline.execute()

            emit_bookmarks_update(room_id)
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

            emit_bookmarks_update(room_id)
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


def ensure_empty_template_exists():
    """Ensure the 'empty' template exists in Redis."""
    redis_client = current_app.extensions["redis"]
    template_key = "template:empty"

    # Check if template already exists
    if not redis_client.exists(template_key):
        # Create the empty template
        template_data = {
            "id": "empty",
            "name": "Empty Room Template",
            "description": "Empty room template",
        }
        redis_client.hset(template_key, mapping=template_data)
        log.info("Created 'empty' template")


@main.route("/api/rooms", methods=["GET"])
def list_rooms():
    """List all active rooms."""
    redis_client = current_app.extensions["redis"]

    # Scan for all room keys to find unique room IDs
    room_ids = set()
    for key in redis_client.scan_iter(match="room:*"):
        # Extract room ID from keys like "room:{room_id}:..."
        parts = key.split(":")
        if len(parts) >= 2:
            room_ids.add(parts[1])

    # Return list of room objects with id and template fields
    rooms = []
    for room_id in sorted(room_ids):
        template_id = redis_client.get(f"room:{room_id}:template") or "empty"
        rooms.append({"id": room_id, "template": template_id})

    return rooms, 200


@main.route("/api/rooms/<string:room_id>", methods=["GET"])
def get_room(room_id):
    """Get details for a specific room."""
    redis_client = current_app.extensions["redis"]

    # Check if room exists
    room_exists = False
    for key in redis_client.scan_iter(match=f"room:{room_id}:*", count=1):
        room_exists = True
        break

    if not room_exists:
        return {"error": "Room not found"}, 404

    template_id = redis_client.get(f"room:{room_id}:template") or "empty"

    return {"id": room_id, "template": template_id}, 200


@main.route("/api/templates", methods=["GET"])
def list_templates():
    """List all available room templates."""
    redis_client = current_app.extensions["redis"]

    # Ensure empty template exists
    ensure_empty_template_exists()

    # Get all templates from Redis
    templates = []
    for key in redis_client.scan_iter(match="template:*"):
        template_data = redis_client.hgetall(key)
        if template_data:
            templates.append(template_data)

    # Sort by id to ensure "empty" comes first
    templates.sort(key=lambda t: (t["id"] != "empty", t["id"]))

    return templates, 200


@main.route("/api/templates/default", methods=["GET"])
def get_default_template():
    """Get the default template."""
    redis_client = current_app.extensions["redis"]

    # Ensure empty template exists
    ensure_empty_template_exists()

    # Get default template ID from Redis, defaulting to "empty"
    default_template_id = redis_client.get("default_template") or "empty"

    # Get the template data
    template_key = f"template:{default_template_id}"
    template_data = redis_client.hgetall(template_key)

    if not template_data:
        return {"error": "Default template not found"}, 404

    return template_data, 200


@main.route("/api/templates/default", methods=["PUT"])
def set_default_template():
    """Set the default template."""
    redis_client = current_app.extensions["redis"]
    data = request.get_json()

    template_id = data.get("template_id")
    if not template_id:
        return {"error": "template_id is required"}, 400

    # Verify template exists
    template_key = f"template:{template_id}"
    if not redis_client.exists(template_key):
        return {"error": f"Template '{template_id}' not found"}, 404

    # Set as default
    redis_client.set("default_template", template_id)

    log.info(f"Set default template to '{template_id}'")
    return {"status": "ok"}, 200


@main.route("/api/rooms/<string:room_id>/promote", methods=["POST"])
def promote_room_to_template(room_id):
    """Promote a room to a template and make it read-only."""
    redis_client = current_app.extensions["redis"]
    data = request.get_json()

    name = data.get("name")
    description = data.get("description")

    if not name or not description:
        return {"error": "name and description are required"}, 400

    # Sanitize room_id: replace ':' with '_' to avoid conflicts with frame index format
    if ":" in room_id:
        original_room_id = room_id
        room_id = room_id.replace(":", "_")
        log.warning(
            f"Template/Room ID '{original_room_id}' contains ':' which is reserved. "
            f"Using sanitized ID '{room_id}' instead."
        )

    # Check if room exists by looking for any room-related keys
    room_exists = False
    for key in redis_client.scan_iter(match=f"room:{room_id}:*", count=1):
        room_exists = True
        break

    if not room_exists:
        return {"error": "Room not found"}, 404

    # Create template from room
    template_key = f"template:{room_id}"
    template_data = {"id": room_id, "name": name, "description": description}
    redis_client.hset(template_key, mapping=template_data)

    # Lock the room permanently by setting a lock without TTL
    # This prevents any client from acquiring the trajectory lock
    lock_key = get_lock_key(room_id, "trajectory:meta")
    redis_client.set(lock_key, "template-lock")

    log.info(
        f"Promoted room '{room_id}' to template '{name}' and locked it permanently"
    )
    return {"status": "ok"}, 200


@main.route("/api/shutdown", methods=["POST"])
def exit_app():
    """Endpoint to gracefully shut down the server. Secured via a shared secret."""
    socketio.stop()
    return {"success": True}


@main.route("/api/rooms/<string:room_id>/schema/<string:category>", methods=["GET"])
def get_room_schema(room_id: str, category: str):
    """Get the schema for a specific room with worker and queue statistics.

    Returns schema along with metadata about each extension:
    - provider: "celery" for server-side extensions, or count of registered workers
    - queueLength: number of queued tasks for this extension
    - idleWorkers: number of idle workers available
    - progressingWorkers: number of workers currently processing tasks
    """
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    # Map category strings to the corresponding imported objects
    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
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
    json_data = request.json
    if json_data is None:
        json_data = {}

    # Try to get userId from query params first, then from JSON body
    user_id = request.args.get("userId") or json_data.pop("userId", None)
    if user_id is None:
        return {"error": "User ID is required"}, 400

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
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
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
    user_id = request.args.get("userId")
    print(
        f"get_extension_data called with userId={user_id}, category={category}, extension={extension} for room {room_id}"
    )

    if not user_id:
        return {"error": "User ID is required"}, 400

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

    for category in ["modifiers", "selections"]:
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
    for category in ["modifiers", "selections"]:
        if is_celery_worker:
            # For celery workers, check all extensions for celery jobs
            from zndraw.extensions.modifiers import modifiers
            from zndraw.extensions.selections import selections

            category_map = {
                "modifiers": modifiers,
                "selections": selections,
            }

            if category in category_map:
                for extension in category_map[category].keys():
                    keys = ExtensionKeys.for_extension(room_id, category, extension)
                    queue_length = redis_client.llen(keys.queue)

                    if queue_length > 0:
                        # Peek at the first job to check if it's a celery job
                        task_data = redis_client.lindex(keys.queue, 0)
                        if task_data:
                            task_info = json.loads(task_data)
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
                        task_info = json.loads(task_data)
                        job_id = task_info.get("jobId")

                        # Get full job details
                        job = JobManager.get_job(redis_client, job_id)
                        if job:
                            # Move worker from idle to progressing
                            redis_client.smove(
                                keys.idle_workers, keys.progressing_workers, worker_id
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
    data = request.get_json() or {}
    if ":" in room_id:
        return {"error": "Room ID cannot contain ':' character"}, 400

    user_name = data.get("userId", "Anonymous")
    client_id = data.get("clientId")
    r = current_app.extensions["redis"]

    # Generate clientId if not provided
    if not client_id:
        import uuid

        client_id = str(uuid.uuid4())
        log.info(f"Generated new clientId: {client_id}")

    # Get or create client data
    client_key = f"client:{client_id}"

    # Update client data (but NOT currentSid - that happens on socket connect)
    r.hset(client_key, "userName", user_name)
    r.hset(client_key, "currentRoom", room_id)

    # Check if room already exists
    room_exists = r.exists(f"room:{room_id}:template")

    # Add client to room's client set (they're joining the room, even if not connected yet)
    r.sadd(f"room:{room_id}:clients", client_id)

    # Generate temporary join token for socket authentication
    import uuid

    join_token = str(uuid.uuid4())
    token_key = f"join_token:{join_token}"
    # Store token with clientId, room, userName (expires in 60 seconds)
    token_data = json.dumps(
        {"clientId": client_id, "roomId": room_id, "userName": user_name}
    )
    r.setex(token_key, 60, token_data)

    log.info(
        f"Client {client_id} ({user_name}) prepared to join room: {room_id}, token: {join_token}"
    )

    response = {
        "status": "ok",
        "clientId": client_id,
        "joinToken": join_token,
        "frameCount": 0,
        "roomId": room_id,
        "template": "empty",
        "selection": None,
        "frame_selection": None,
        "created": True,
        "presenter-lock": False,
        "step": None,
    }

    response["created"] = not room_exists

    if not room_exists:
        if "template" not in data:
            # Use default template
            from .routes import ensure_empty_template_exists

            ensure_empty_template_exists()
            template_id = r.get("default_template") or "empty"
        elif data["template"] is None:
            # Explicit None means use "empty" template
            template_id = "empty"
        else:
            # Use specified template
            template_id = data["template"]
            # Verify template exists, fall back to "empty" if it doesn't
            template_key = f"template:{template_id}"
            if not r.exists(template_key):
                log.warning(
                    f"Template '{template_id}' not found for room '{room_id}', falling back to 'empty'"
                )
                template_id = "empty"

        # Store the template used for this room
        r.set(f"room:{room_id}:template", template_id)
        log.info(f"New room '{room_id}' created from template '{template_id}'")

        # Initialize room with template frames if not empty
        if template_id != "empty":
            template_indices_key = f"room:{template_id}:trajectory:indices"
            template_frames = r.zrange(template_indices_key, 0, -1, withscores=True)

            if template_frames:
                room_indices_key = f"room:{room_id}:trajectory:indices"
                # Copy template frames with the new room_id:index format
                for member, score in template_frames:
                    r.zadd(room_indices_key, {member: score})
                log.info(
                    f"Copied {len(template_frames)} frames from template '{template_id}' to room '{room_id}'"
                )
        response["template"] = template_id

    indices_key = f"room:{room_id}:trajectory:indices"
    response["frameCount"] = r.zcard(indices_key)

    # Initialize current_frame to 0 for newly created rooms
    if not room_exists:
        r.set(f"room:{room_id}:current_frame", 0)
        log.info(f"Initialized current_frame to 0 for new room '{room_id}'")

        # Initialize default settings for new room
        from zndraw.settings import RoomConfig

        default_config = RoomConfig()
        config_dict = default_config.model_dump()

        # Store each settings category
        for category_name, category_data in config_dict.items():
            r.hset(
                f"room:{room_id}:user:{user_name}:settings",
                category_name,
                json.dumps(category_data),
            )
        log.info(f"Initialized default settings for new room '{room_id}'")

    selection = r.get(f"room:{room_id}:selection:default")
    response["selection"] = json.loads(selection) if selection else None

    frame_selection = r.get(f"room:{room_id}:frame_selection:default")
    response["frame_selection"] = (
        json.loads(frame_selection) if frame_selection else None
    )

    presenter_lock = r.get(f"room:{room_id}:presenter_lock")
    response["presenter-lock"] = presenter_lock

    step = r.get(f"room:{room_id}:current_frame")
    response["step"] = int(step) if step is not None else 0

    bookmarks_key = f"room:{room_id}:bookmarks"
    physical_bookmarks = r.hgetall(bookmarks_key)

    if physical_bookmarks:
        frame_mapping = r.zrange(indices_key, 0, -1)

        # Build reverse mapping: physical_key -> logical_index
        physical_to_logical = {
            physical_key: idx for idx, physical_key in enumerate(frame_mapping)
        }

        # Convert bookmarks from physical to logical indices
        logical_bookmarks = {}
        for physical_key, label in physical_bookmarks.items():
            if physical_key in physical_to_logical:
                logical_idx = physical_to_logical[physical_key]
                logical_bookmarks[logical_idx] = label
            # If physical_key not in mapping, skip it (orphaned bookmark)

        response["bookmarks"] = logical_bookmarks
    else:
        response["bookmarks"] = None

    # Fetch all settings for the user
    settings_data = {}
    settings_hash = r.hgetall(f"room:{room_id}:user:{user_name}:settings")
    for setting_name, setting_value in settings_hash.items():
        settings_data[setting_name] = json.loads(setting_value)
    response["settings"] = settings_data

    return response


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
