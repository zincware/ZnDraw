"""Frame and trajectory operation routes.

Handles all frame-related operations including getting, creating, updating,
deleting frames and downloading trajectories.
"""

import json
import logging
import traceback

import msgpack
import zarr
from flask import Blueprint, Response, current_app, request, send_file

from zndraw.storage import decode_data, encode_data

from .frame_index_manager import FrameIndexManager
from .route_utils import (
    check_room_locked,
    emit_bookmarks_invalidate,
    emit_frames_invalidate,
    emit_len_frames_update,
    get_storage,
    remove_bookmark_at_index,
    shift_bookmarks_on_delete,
    shift_bookmarks_on_insert,
)

log = logging.getLogger(__name__)

frames = Blueprint("frames", __name__)


@frames.route("/api/rooms/<string:room_id>/frames", methods=["GET"])
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


@frames.route(
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


@frames.route("/api/rooms/<string:room_id>/frames", methods=["DELETE"])
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


@frames.route("/api/rooms/<string:room_id>/frames", methods=["POST"])
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


@frames.route("/api/rooms/<string:room_id>/frames/bulk", methods=["PATCH"])
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


@frames.route("/api/rooms/<string:room_id>/download", methods=["GET"])
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


# =============================================================================
# Extension Analytics Endpoints
# =============================================================================
