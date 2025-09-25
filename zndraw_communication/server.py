from flask import Flask, Response, request
from flask_socketio import SocketIO, join_room, rooms, leave_room
import redis
import numpy as np
import zarr
import typing as t
import uuid
import logging
import msgpack
import shutil
import json

# --- Logging Setup ---
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)


app = Flask(__name__)
socketio = SocketIO(app, cors_allowed_origins="*")
r = redis.Redis(decode_responses=True)


# --- Helper Functions ---
def get_project_room_from_session(sid: str) -> t.Optional[str]:
    """Finds the project room a client has joined."""
    for room in rooms(sid=sid):
        if room != sid:
            return room
    return None

def get_lock_key(room: str, target: str) -> str:
    """Constructs a standardized Redis key for a lock."""
    return f"room:{room}:lock:{target}"

def get_zarr_store_path(room_id: str) -> str:
    """Returns the path to the Zarr store for a given room."""
    # In a real app, this path might come from a config file.
    return f"data/{room_id}.zarr"


# --- HTTP Data Endpoints ---
@app.route("/frame/<string:room_id>/<int:frame_id>")
def get_frame(room_id, frame_id):
    """Serves a single frame's data from the room's Zarr store."""
    # Get keys parameter from query string
    keys_param = request.args.get('keys')
    requested_keys = keys_param.split(',') if keys_param else None

    store_path = get_zarr_store_path(room_id)
    try:
        root = zarr.group(store_path)

        # Get logical-to-physical mapping from Redis
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return Response("No frames found in room", status=404)

        if frame_id >= len(frame_mapping):
            error_data = {"error": f"Frame {frame_id} not found, max frame: {len(frame_mapping)-1}", "type": "IndexError"}
            return Response(json.dumps(error_data), status=404, content_type='application/json')

        # Get the physical index for this logical frame
        physical_index = int(frame_mapping[frame_id])

        # Build response dict with arrays for this frame
        frame_data = {}

        # Get all keys from zarr store (arrays) and metadata
        available_keys = set(root.keys())
        metadata_keys = set()

        # If metadata exists, extract the keys it contains
        if '_metadata' in root:
            metadata_dataset = root['_metadata']
            if physical_index < metadata_dataset.shape[0]:
                metadata_array = metadata_dataset[physical_index]
                try:
                    # Reconstruct the JSON string from the metadata array
                    # metadata_array is a numpy array with Unicode strings
                    json_str = metadata_array.item() if metadata_array.size > 0 else '{}'
                    metadata_dict = json.loads(json_str)
                    metadata_keys = set(metadata_dict.keys())
                except Exception as e:
                    log.debug(f"Failed to parse metadata: {e}")
                    pass

        # Determine which keys to process
        if requested_keys:
            keys_to_process = requested_keys
            # Validate that all requested keys exist
            all_available_keys = available_keys | metadata_keys
            missing_keys = set(requested_keys) - all_available_keys
            if missing_keys:
                error_data = {"error": f"Key(s) not found: {', '.join(sorted(missing_keys))}", "type": "KeyError"}
                return Response(json.dumps(error_data), status=404, content_type='application/json')
        else:
            keys_to_process = list(available_keys) + list(metadata_keys)

        # Process regular zarr arrays
        for key in keys_to_process:
            if key in root and key != '_metadata':
                dataset = root[key]
                if physical_index < dataset.shape[0]:
                    frame_data[key] = dataset[physical_index]

        # Process metadata keys if metadata exists
        if '_metadata' in root and metadata_keys:
            metadata_dataset = root['_metadata']
            if physical_index < metadata_dataset.shape[0]:
                metadata_array = metadata_dataset[physical_index]
                try:
                    json_str = metadata_array.item() if metadata_array.size > 0 else '{}'
                    metadata_dict = json.loads(json_str)

                    # Filter metadata to only include requested keys that exist in metadata
                    filtered_metadata = {}
                    for key in keys_to_process:
                        if key in metadata_dict:
                            filtered_metadata[key] = metadata_dict[key]

                    # If we have filtered metadata to include, add it back to frame_data as _metadata
                    if filtered_metadata:
                        # Convert metadata back to the format expected by client
                        filtered_json_str = json.dumps(filtered_metadata)
                        filtered_json_array = np.array([filtered_json_str], dtype='U')
                        frame_data['_metadata'] = filtered_json_array
                except:
                    pass

        # Serialize using msgpack with bytes and shape info
        serialized_data = {}
        for key, array in frame_data.items():
            serialized_data[key] = {
                'data': array.tobytes(),
                'shape': array.shape,
                'dtype': str(array.dtype)
            }

        packed_data = msgpack.packb(serialized_data)
        return Response(packed_data, content_type='application/octet-stream')
    except (IOError, KeyError, Exception) as e:
        log.error(f"Error retrieving frame {frame_id} from room '{room_id}': {e}")
        return Response(f"Room '{room_id}' not found or is invalid.", status=404)

@app.route("/frames/<string:room_id>", methods=["POST"])
def get_frames(room_id):
    """Serves multiple frames' data from the room's Zarr store using either indices or slice parameters."""
    try:
        # Parse the request data
        request_data = request.get_json()
        if request_data is None:
            return {"error": "Request body required"}, 400

        store_path = get_zarr_store_path(room_id)
        root = zarr.group(store_path)

        # Get logical-to-physical mapping from Redis
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {"error": "No frames found in room"}, 404

        max_frame = len(frame_mapping) - 1

        # Determine frame indices based on request parameters
        if 'indices' in request_data:
            # Direct list of indices
            frame_indices = request_data['indices']
            if not isinstance(frame_indices, list):
                return {"error": "Indices must be a list"}, 400

            # Validate frame indices
            for frame_id in frame_indices:
                if not isinstance(frame_id, int) or frame_id < 0 or frame_id > max_frame:
                    error_data = {"error": f"Invalid frame index {frame_id}, valid range: 0-{max_frame}", "type": "IndexError"}
                    return Response(json.dumps(error_data), status=404, content_type='application/json')

        else:
            # Default to slice behavior for any remaining cases (including empty payload)
            # This handles slice parameters and slice(None, None, None) which sends empty payload
            start = request_data.get('start', 0)
            stop = request_data.get('stop', len(frame_mapping))
            step = request_data.get('step', 1)

            # Validate slice parameters
            if not all(isinstance(x, int) for x in [start, stop, step]):
                return {"error": "start, stop, and step must be integers"}, 400

            if step == 0:
                return {"error": "step cannot be zero"}, 400

            # Generate frame indices from slice
            try:
                frame_indices = list(range(start, stop, step))
                # Filter out invalid indices
                frame_indices = [i for i in frame_indices if 0 <= i <= max_frame]
            except ValueError as e:
                return {"error": f"Invalid slice parameters: {e}"}, 400

        # Get keys parameter if specified
        requested_keys = request_data.get('keys')

        # Validate requested keys before processing any frames
        if requested_keys:
            # For frames endpoint, we need comprehensive validation similar to single frame endpoint
            available_keys = set(root.keys())
            metadata_keys = set()

            # If metadata exists, get the keys from the first frame to check availability
            if '_metadata' in root:
                metadata_dataset = root['_metadata']
                if len(frame_mapping) > 0 and metadata_dataset.shape[0] > 0:
                    try:
                        physical_index = int(frame_mapping[0])
                        if physical_index < metadata_dataset.shape[0]:
                            metadata_array = metadata_dataset[physical_index]
                            json_str = metadata_array.item() if metadata_array.size > 0 else '{}'
                            metadata_dict = json.loads(json_str)
                            metadata_keys = set(metadata_dict.keys())
                    except:
                        pass

            # Validate that all requested keys exist
            all_available_keys = available_keys | metadata_keys
            missing_keys = set(requested_keys) - all_available_keys
            if missing_keys:
                error_data = {"error": f"Key(s) not found: {', '.join(sorted(missing_keys))}", "type": "KeyError"}
                return Response(json.dumps(error_data), status=404, content_type='application/json')

        # Build response with all requested frames (empty list if no indices)
        frames_data = []
        if not frame_indices:
            # Return empty result for valid but empty requests
            packed_data = msgpack.packb(frames_data)
            return Response(packed_data, content_type='application/octet-stream')
        for frame_id in frame_indices:
            # Get the physical index for this logical frame
            physical_index = int(frame_mapping[frame_id])

            # Build frame data dict with arrays for this frame
            frame_data = {}

            # Determine which keys to process
            keys_to_process = requested_keys if requested_keys else root.keys()

            for key in keys_to_process:
                if key in root:
                    dataset = root[key]
                    if physical_index < dataset.shape[0]:
                        frame_data[key] = dataset[physical_index]

            # Serialize frame using msgpack with bytes and shape info
            serialized_frame = {}
            for key, array in frame_data.items():
                serialized_frame[key] = {
                    'data': array.tobytes(),
                    'shape': array.shape,
                    'dtype': str(array.dtype)
                }

            frames_data.append(serialized_frame)

        packed_data = msgpack.packb(frames_data)
        return Response(packed_data, content_type='application/octet-stream')
    except Exception as e:
        log.error(f"Error retrieving frames from room '{room_id}': {e}")
        return {"error": "Internal server error"}, 500

@app.route("/rooms/<string:room_id>/frames", methods=["POST"])
def append_frame(room_id):
    """Appends a new frame. Authorized via a short-lived Bearer token."""
    auth_header = request.headers.get("Authorization")
    if not auth_header or not auth_header.startswith("Bearer "):
        return {"error": "Authorization token is missing or invalid"}, 401

    token = auth_header.split(" ")[1]
    token_key = f"room:{room_id}:upload_token:{token}"

    # Get token metadata
    token_data = r.hgetall(token_key)
    if not token_data:
        return {"error": "Token is invalid or has expired"}, 403

    sid_from_token = token_data.get("sid")
    action = token_data.get("action", "append")
    target_frame_id = int(token_data.get("frame_id", -1)) if token_data.get("frame_id") != "-1" else None

    r.delete(token_key) # Invalidate the token after first use

    lock_key = get_lock_key(room_id, "trajectory:meta")
    if r.get(lock_key) != sid_from_token:
        return {"error": "Client does not hold the trajectory lock"}, 403

    try:
        # Unpack the msgpack data
        serialized_data = msgpack.unpackb(request.data, strict_map_key=False)

        store_path = get_zarr_store_path(room_id)
        root = zarr.group(store_path)

        indices_key = f"room:{room_id}:trajectory:indices"

        if action == "replace":
            # Replace operation: update the physical data but keep the same logical mapping
            frame_mapping = r.zrange(indices_key, 0, -1)
            if target_frame_id >= len(frame_mapping):
                return {"error": f"Frame {target_frame_id} does not exist"}, 404

            # Get the physical index for this logical frame
            physical_index = int(frame_mapping[target_frame_id])

            # Process each array in the frame data
            for key, array_info in serialized_data.items():
                # Reconstruct numpy array from bytes, shape, and dtype
                data_bytes = array_info['data']
                shape = tuple(array_info['shape'])
                dtype = array_info['dtype']
                array = np.frombuffer(data_bytes, dtype=dtype).reshape(shape)

                # Update the existing data at the physical index
                if key not in root:
                    return {"error": f"Array '{key}' does not exist for replacement"}, 400

                dataset = root[key]
                if physical_index >= dataset.shape[0]:
                    return {"error": f"Physical index {physical_index} out of range for array '{key}'"}, 400

                dataset[physical_index] = array

            socketio.emit("trajectory_updated", {"action": "replace", "frame_id": target_frame_id}, to=room_id, skip_sid=sid_from_token)

            log.info(f"Replaced frame {target_frame_id} (physical: {physical_index}) in room '{room_id}' with keys: {list(serialized_data.keys())}")
            return {"success": True, "replaced_frame": target_frame_id}

        elif action == "extend":
            # Extend operation: add multiple frames in one go
            # Find next available physical indices
            used_physical_indices = [int(x) for x in r.zrange(indices_key, 0, -1)]
            next_physical_index = max(used_physical_indices) + 1 if used_physical_indices else 0

            # serialized_data should be a list of dictionaries for multiple frames
            if not isinstance(serialized_data, list):
                return {"error": "For extend action, data must be a list of frame dictionaries"}, 400

            num_frames = len(serialized_data)
            new_indices = []

            # Process each frame in the list
            for frame_idx, frame_data in enumerate(serialized_data):
                current_physical_index = next_physical_index + frame_idx

                # Process each array in the frame data
                for key, array_info in frame_data.items():
                    # Reconstruct numpy array from bytes, shape, and dtype
                    data_bytes = array_info['data']
                    shape = tuple(array_info['shape'])
                    dtype = array_info['dtype']
                    array = np.frombuffer(data_bytes, dtype=dtype).reshape(shape)

                    # Create array if it doesn't exist
                    if key not in root:
                        # Create with expandable first dimension for frames
                        initial_shape = (current_physical_index + num_frames - frame_idx,) + array.shape
                        chunks = (1,) + array.shape
                        dataset = root.create_array(
                            name=key,
                            shape=initial_shape,
                            chunks=chunks,
                            dtype=array.dtype
                        )
                        dataset[current_physical_index] = array
                        log.info(f"Created new array '{key}' with shape {initial_shape}")
                    else:
                        dataset = root[key]
                        # Resize array to accommodate all new frames if needed
                        required_size = current_physical_index + 1
                        if required_size > dataset.shape[0]:
                            new_shape = (required_size,) + dataset.shape[1:]
                            dataset.resize(new_shape)
                        dataset[current_physical_index] = array

                # Add the physical index to the logical sequence
                logical_position = len(used_physical_indices) + frame_idx
                new_indices.append(logical_position)
                r.zadd(indices_key, {str(current_physical_index): logical_position})

            socketio.emit("trajectory_updated", {"action": "extend", "new_indices": new_indices}, to=room_id, skip_sid=sid_from_token)

            log.info(f"Extended trajectory with {num_frames} frames (physical: {next_physical_index}-{next_physical_index + num_frames - 1}) to room '{room_id}'")
            return {"success": True, "new_indices": new_indices}

        elif action == "insert":
            # Insert operation: add new data and shift existing logical indices
            insert_position = int(token_data.get("insert_position", 0))

            # Get current frame mapping and validate insert position
            used_physical_indices = [int(x) for x in r.zrange(indices_key, 0, -1)]
            current_length = len(used_physical_indices)

            if insert_position < 0 or insert_position > current_length:
                return {"error": f"Insert position {insert_position} out of range [0, {current_length}]"}, 400

            # Find next available physical index for the new frame
            next_physical_index = max(used_physical_indices) + 1 if used_physical_indices else 0

            # Process the new frame data
            for key, array_info in serialized_data.items():
                # Reconstruct numpy array from bytes, shape, and dtype
                data_bytes = array_info['data']
                shape = tuple(array_info['shape'])
                dtype = array_info['dtype']
                array = np.frombuffer(data_bytes, dtype=dtype).reshape(shape)

                # Create array if it doesn't exist
                if key not in root:
                    # Create with expandable first dimension for frames
                    initial_shape = (next_physical_index + 1,) + array.shape
                    chunks = (1,) + array.shape
                    dataset = root.create_array(
                        name=key,
                        shape=initial_shape,
                        chunks=chunks,
                        dtype=array.dtype
                    )
                    dataset[next_physical_index] = array
                    log.info(f"Created new array '{key}' with shape {initial_shape}")
                else:
                    dataset = root[key]
                    # Resize array to accommodate new frame if needed
                    if next_physical_index >= dataset.shape[0]:
                        new_shape = (next_physical_index + 1,) + dataset.shape[1:]
                        dataset.resize(new_shape)
                    dataset[next_physical_index] = array

            # Get all current mappings and rebuild them with the insertion
            all_mappings = r.zrange(indices_key, 0, -1, withscores=True)

            # Clear the current mappings
            r.delete(indices_key)

            # Rebuild mappings with the new frame inserted
            pipeline = r.pipeline()

            # Add all existing frames, shifting those at position >= insert_position
            for physical_idx_str, logical_pos in all_mappings:
                if logical_pos >= insert_position:
                    # Shift this frame one position to the right
                    pipeline.zadd(indices_key, {physical_idx_str: logical_pos + 1})
                else:
                    # Keep this frame at its current position
                    pipeline.zadd(indices_key, {physical_idx_str: logical_pos})

            # Add the new frame at the insert position
            pipeline.zadd(indices_key, {str(next_physical_index): insert_position})
            pipeline.execute()

            socketio.emit("trajectory_updated", {"action": "insert", "position": insert_position}, to=room_id, skip_sid=sid_from_token)

            log.info(f"Inserted frame at position {insert_position} (physical: {next_physical_index}) in room '{room_id}' with keys: {list(serialized_data.keys())}")
            return {"success": True, "inserted_position": insert_position}

        else:
            # Append operation: add new physical data and update logical mapping
            # Find next available physical index
            used_physical_indices = [int(x) for x in r.zrange(indices_key, 0, -1)]
            next_physical_index = max(used_physical_indices) + 1 if used_physical_indices else 0

            # Process each array in the frame data
            for key, array_info in serialized_data.items():
                # Reconstruct numpy array from bytes, shape, and dtype
                data_bytes = array_info['data']
                shape = tuple(array_info['shape'])
                dtype = array_info['dtype']
                array = np.frombuffer(data_bytes, dtype=dtype).reshape(shape)

                # Create array if it doesn't exist
                if key not in root:
                    # Create with expandable first dimension for frames
                    initial_shape = (next_physical_index + 1,) + array.shape
                    chunks = (1,) + array.shape
                    dataset = root.create_array(
                        name=key,
                        shape=initial_shape,
                        chunks=chunks,
                        dtype=array.dtype
                    )
                    dataset[next_physical_index] = array
                    log.info(f"Created new array '{key}' with shape {initial_shape}")
                else:
                    dataset = root[key]
                    # Resize array to accommodate new frame if needed
                    if next_physical_index >= dataset.shape[0]:
                        new_shape = (next_physical_index + 1,) + dataset.shape[1:]
                        dataset.resize(new_shape)
                    dataset[next_physical_index] = array

            # Add the physical index to the logical sequence
            logical_position = len(used_physical_indices)  # This will be the next logical position
            r.zadd(indices_key, {str(next_physical_index): logical_position})

            socketio.emit("trajectory_updated", {"action": "append", "new_index": logical_position}, to=room_id, skip_sid=sid_from_token)

            log.info(f"Appended frame {logical_position} (physical: {next_physical_index}) to room '{room_id}' with keys: {list(serialized_data.keys())}")
            return {"success": True, "new_index": logical_position}
    except Exception as e:
        log.error(f"Failed to write to Zarr store: {e}")
        return {"error": "Failed to write to data store"}, 500


# --- Socket.IO Control Endpoints ---
@socketio.on('disconnect')
def handle_disconnect():
    sid = request.sid
    log.info(f"Client disconnected: {sid}")
    lock_keys = r.scan_iter(f"*:lock:*")
    for key in lock_keys:
        if r.get(key) == sid:
            log.warning(f"Cleaning up orphaned lock '{key}' held by disconnected client {sid}")
            r.delete(key)

@socketio.on("join_room")
def handle_join(data):
    room = data["room"]
    sid = request.sid
    if previous_room := get_project_room_from_session(sid):
        leave_room(previous_room)
        log.info(f"Client {sid} left room: {previous_room}")
        
    join_room(room)
    log.info(f"Client {sid} joined room: {room}")

@socketio.on("lock:acquire")
def acquire_lock(data):
    sid = request.sid
    target = data.get("target")
    room = get_project_room_from_session(sid)

    if not room or not target:
        return {"success": False, "error": "Room or target missing"}

    lock_key = get_lock_key(room, target)
    if r.set(lock_key, sid, nx=True, ex=60):
        log.info(f"Lock acquired for '{target}' in room '{room}' by {sid}")
        return {"success": True}
    else:
        log.info(f"Lock for '{target}' in room '{room}' already held by {r.get(lock_key)}, denied for {sid}")
        return {"success": False}

@socketio.on("lock:release")
def release_lock(data):
    sid = request.sid
    target = data.get("target")
    room = get_project_room_from_session(sid)

    if not room or not target:
        return {"success": False, "error": "Room or target missing"}
    
    lock_key = get_lock_key(room, target)
    if r.get(lock_key) == sid:
        r.delete(lock_key)
        log.info(f"Lock released for '{target}' in room '{room}' by {sid}")
        return {"success": True}
    
    log.warning(f"Failed release: Lock for '{target}' in room '{room}' not held by {sid}")
    return {"success": False}

@socketio.on("upload:prepare")
def handle_upload_prepare(data):
    sid = request.sid
    room = get_project_room_from_session(sid)
    action = data.get("action", "append")  # Default to append for backward compatibility
    frame_id = data.get("frame_id")  # For replace operations
    insert_position = data.get("insert_position")  # For insert operations

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    lock_key = get_lock_key(room, "trajectory:meta")
    if r.get(lock_key) != sid:
        return {"success": False, "error": "Client does not hold the trajectory lock."}

    # For replace operations, validate the frame exists
    if action == "replace":
        if frame_id is None:
            return {"success": False, "error": "frame_id is required for replace operations"}

        indices_key = f"room:{room}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping or frame_id >= len(frame_mapping):
            return {"success": False, "error": f"Frame {frame_id} does not exist"}

    token = str(uuid.uuid4())
    token_key = f"room:{room}:upload_token:{token}"
    # Store additional metadata with the token
    token_data = {
        "sid": sid,
        "action": action,
        "frame_id": str(frame_id) if frame_id is not None else "-1"
    }

    # Add insert_position for insert operations
    if insert_position is not None:
        token_data["insert_position"] = str(insert_position)
    r.hset(token_key, mapping=token_data)
    r.expire(token_key, 60)

    log.info(f"Issued {action} token for room '{room}' to {sid}" + (f" (frame {frame_id})" if frame_id is not None else ""))
    return {"success": True, "token": token}

@socketio.on("frames:count")
def handle_len_frames(data):
    sid = request.sid
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

@socketio.on("frame:delete")
def handle_delete_frame(data):
    sid = request.sid
    room = get_project_room_from_session(sid)
    frame_id = data.get("frame_id")

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    if frame_id is None:
        return {"success": False, "error": "frame_id is required"}

    lock_key = get_lock_key(room, "trajectory:meta")
    if r.get(lock_key) != sid:
        return {"success": False, "error": "Client does not hold the trajectory lock."}

    try:
        indices_key = f"room:{room}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {"success": False, "error": "No frames found in room"}

        if frame_id >= len(frame_mapping):
            return {"success": False, "error": f"Frame {frame_id} not found, max frame: {len(frame_mapping)-1}"}

        # Get the physical index that we're "deleting" (just removing from mapping)
        physical_index_to_remove = int(frame_mapping[frame_id])

        # Remove the mapping entry for this logical position
        # We need to rebuild the mapping without this entry
        remaining_physical_indices = frame_mapping[:frame_id] + frame_mapping[frame_id+1:]

        # Clear and rebuild the Redis mapping
        r.delete(indices_key)
        for logical_pos, physical_idx_str in enumerate(remaining_physical_indices):
            r.zadd(indices_key, {physical_idx_str: logical_pos})

        socketio.emit("trajectory_updated", {"action": "delete", "frame_id": frame_id}, to=room, skip_sid=sid)

        log.info(f"Deleted logical frame {frame_id} (physical: {physical_index_to_remove}) from room '{room}'. Physical data preserved.")
        return {"success": True, "deleted_frame": frame_id, "physical_preserved": physical_index_to_remove}
    except Exception as e:
        log.error(f"Failed to delete frame: {e}")
        return {"success": False, "error": "Failed to delete frame"}


def main():
    """Main entry point for the ZnDraw server."""
    try:
        log.info("Starting ZnDraw Server")
        socketio.run(app, debug=True, host="0.0.0.0", port=5000)
    finally:
        r.flushall()
        # Clean up Zarr data directory on shutdown for this example
        shutil.rmtree("data", ignore_errors=True)


if __name__ == '__main__':
    main()
