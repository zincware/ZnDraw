import json
import logging
import shutil
import typing as t
import uuid

import msgpack
import redis
import zarr
from flask import Flask, Response, request
from flask_socketio import SocketIO, join_room, leave_room, rooms

from zndraw_communication.storage import ZarrStorageSequence, decode_data, encode_data

# --- Logging Setup ---
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
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
        storage = ZarrStorageSequence(root)

        # Get logical-to-physical mapping from Redis
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {"error": "No frames found in room"}, 404

        max_frame = len(frame_mapping) - 1

        # Determine frame indices based on request parameters
        if "indices" in request_data:
            # Direct list of indices
            frame_indices = request_data["indices"]
            if not isinstance(frame_indices, list):
                return {"error": "Indices must be a list"}, 400

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
            start = request_data.get("start", 0)
            stop = request_data.get("stop", len(frame_mapping))
            step = request_data.get("step", 1)

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
        requested_keys = request_data.get("keys")
        # error_data = {"error": f"Key(s) not found: {', '.join(sorted(missing_keys))}", "type": "KeyError"}
        # return Response(json.dumps(error_data), status=404, content_type='application/json')

        # TODO: requested keys and KeyError handling
        # TODO: instead of iterate load at once
        try:
            frames_data = []
            for frame_id in frame_indices:
                # Get the physical index for this logical frame
                physical_index = int(frame_mapping[frame_id])
                frame_data = storage.get(physical_index, keys=requested_keys)
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
        return Response(
            json.dumps(error_data), status=500, content_type="application/json"
        )


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
    target_frame_id = (
        int(token_data.get("frame_id", -1))
        if token_data.get("frame_id") != "-1"
        else None
    )

    r.delete(token_key)  # Invalidate the token after first use

    lock_key = get_lock_key(room_id, "trajectory:meta")
    if r.get(lock_key) != sid_from_token:
        return {"error": "Client does not hold the trajectory lock"}, 403

    try:
        # Unpack the msgpack data
        serialized_data = msgpack.unpackb(request.data, strict_map_key=False)

        store_path = get_zarr_store_path(room_id)
        root = zarr.group(store_path)
        storage = ZarrStorageSequence(root)

        indices_key = f"room:{room_id}:trajectory:indices"

        if action == "replace":
            frame_mapping = r.zrange(indices_key, 0, -1)

            # Validate target_frame_id
            if target_frame_id is None or not (
                0 <= target_frame_id < len(frame_mapping)
            ):
                return {
                    "error": f"Invalid or missing frame_id for replace. Valid range: 0-{len(frame_mapping) - 1}"
                }, 404

            new_physical_index = len(storage)
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            old_physical_index_to_unmap = frame_mapping[target_frame_id]

            pipeline = r.pipeline()
            pipeline.zrem(indices_key, old_physical_index_to_unmap)
            pipeline.zadd(indices_key, {str(new_physical_index): target_frame_id})
            pipeline.execute()

            log.info(
                f"Replaced frame {target_frame_id} (old physical: {old_physical_index_to_unmap}, new physical: {new_physical_index}) in room '{room_id}'"
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

            # 3. Create the new logical-to-physical mapping for all new frames
            new_mapping = {
                str(start_physical_pos + i): start_logical_pos + i
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
            # Insert operation: add a new physical frame and insert it into the logical sequence
            insert_position = int(token_data.get("insert_position", 0))
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

            # Add the new frame at the correct logical position
            pipeline.zadd(indices_key, {str(new_physical_index): insert_position})
            pipeline.execute()

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

            # 3. Add the new physical index to the logical sequence
            r.zadd(indices_key, {str(new_physical_index): logical_position})

            log.info(
                f"Appended frame {logical_position} (physical: {new_physical_index}) to room '{room_id}'"
            )
            return {"success": True, "new_index": logical_position}

        else:
            # Default case for any unknown actions
            return {"error": f"The requested action '{action}' is not supported."}, 400
    except Exception as e:
        log.error(f"Failed to write to Zarr store: {e}")
        return {"error": "Failed to write to data store"}, 500


# --- Socket.IO Control Endpoints ---
@socketio.on("disconnect")
def handle_disconnect():
    sid = request.sid
    log.info(f"Client disconnected: {sid}")
    lock_keys = r.scan_iter(f"*:lock:*")
    for key in lock_keys:
        if r.get(key) == sid:
            log.warning(
                f"Cleaning up orphaned lock '{key}' held by disconnected client {sid}"
            )
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
        log.info(
            f"Lock for '{target}' in room '{room}' already held by {r.get(lock_key)}, denied for {sid}"
        )
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

    log.warning(
        f"Failed release: Lock for '{target}' in room '{room}' not held by {sid}"
    )
    return {"success": False}


@socketio.on("upload:prepare")
def handle_upload_prepare(data):
    sid = request.sid
    room = get_project_room_from_session(sid)
    action = data.get("action")  # Default to append for backward compatibility
    frame_id = data.get("frame_id")  # For replace operations
    insert_position = data.get("insert_position")  # For insert operations

    if action not in {"append", "replace", "insert", "extend"}:
        return {"success": False, "error": "Invalid action specified"}

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    lock_key = get_lock_key(room, "trajectory:meta")
    if r.get(lock_key) != sid:
        return {"success": False, "error": "Client does not hold the trajectory lock."}

    # For replace operations, validate the frame exists
    if action == "replace":
        if frame_id is None:
            return {
                "success": False,
                "error": "frame_id is required for replace operations",
            }

        indices_key = f"room:{room}:trajectory:indices"
        # Use zcard for an efficient count of logical frames
        num_frames = r.zcard(indices_key)

        if frame_id >= num_frames:
            return {
                "success": False,
                "error": f"Frame {frame_id} does not exist. Max index is {num_frames - 1}.",
            }

    token = str(uuid.uuid4())
    token_key = f"room:{room}:upload_token:{token}"
    # Store additional metadata with the token
    token_data = {
        "sid": sid,
        "action": action,
        "frame_id": str(frame_id) if frame_id is not None else "-1",
    }

    # Add insert_position for insert operations
    if insert_position is not None:
        token_data["insert_position"] = str(insert_position)
    r.hset(token_key, mapping=token_data)
    r.expire(token_key, 60)

    log.info(
        f"Issued {action} token for room '{room}' to {sid}"
        + (f" (frame {frame_id})" if frame_id is not None else "")
    )
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
            return {
                "success": False,
                "error": f"Frame {frame_id} not found, max frame: {len(frame_mapping) - 1}",
            }

        # Get the physical index that we're "deleting" (just removing from mapping)
        physical_index_to_remove = int(frame_mapping[frame_id])

        # Remove the mapping entry for this logical position
        # We need to rebuild the mapping without this entry
        remaining_physical_indices = (
            frame_mapping[:frame_id] + frame_mapping[frame_id + 1 :]
        )

        # Clear and rebuild the Redis mapping
        r.delete(indices_key)
        for logical_pos, physical_idx_str in enumerate(remaining_physical_indices):
            r.zadd(indices_key, {physical_idx_str: logical_pos})

        log.info(
            f"Deleted logical frame {frame_id} (physical: {physical_index_to_remove}) from room '{room}'. Physical data preserved."
        )
        return {
            "success": True,
            "deleted_frame": frame_id,
            "physical_preserved": physical_index_to_remove,
        }
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


if __name__ == "__main__":
    main()
