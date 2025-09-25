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
    store_path = get_zarr_store_path(room_id)
    try:
        root = zarr.group(store_path)

        # Check if frame exists in Redis index
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_count = r.zcard(indices_key)

        if frame_count == 0:
            return Response("No frames found in room", status=404)

        if frame_id >= frame_count:
            return Response(f"Frame {frame_id} not found, max frame: {frame_count-1}", status=404)

        # Build response dict with all arrays for this frame
        frame_data = {}
        for key in root.keys():
            dataset = root[key]
            if frame_id < dataset.shape[0]:
                frame_data[key] = dataset[frame_id]

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

@app.route("/rooms/<string:room_id>/frames", methods=["POST"])
def append_frame(room_id):
    """Appends a new frame. Authorized via a short-lived Bearer token."""
    auth_header = request.headers.get("Authorization")
    if not auth_header or not auth_header.startswith("Bearer "):
        return {"error": "Authorization token is missing or invalid"}, 401

    token = auth_header.split(" ")[1]
    token_key = f"room:{room_id}:upload_token:{token}"

    sid_from_token = r.get(token_key)
    if not sid_from_token:
        return {"error": "Token is invalid or has expired"}, 403

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
        last_item = r.zrevrange(indices_key, 0, 0)
        new_index = int(last_item[0]) + 1 if last_item else 0

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
                initial_shape = (1,) + array.shape
                chunks = (1,) + array.shape
                dataset = root.create_array(
                    name=key,
                    shape=initial_shape,
                    chunks=chunks,
                    dtype=array.dtype
                )
                dataset[0] = array
                log.info(f"Created new array '{key}' with shape {initial_shape}")
            else:
                dataset = root[key]
                # Resize array to accommodate new frame
                current_shape = dataset.shape
                new_shape = (new_index + 1,) + current_shape[1:]
                dataset.resize(new_shape)
                dataset[new_index] = array

        r.zadd(indices_key, {str(new_index): new_index})

        socketio.emit("trajectory_updated", {"action": "append", "new_index": new_index}, to=room_id, skip_sid=sid_from_token)

        log.info(f"Appended frame {new_index} to room '{room_id}' with keys: {list(serialized_data.keys())}")
        return {"success": True, "new_index": new_index}
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

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    lock_key = get_lock_key(room, "trajectory:meta")
    if r.get(lock_key) != sid:
        return {"success": False, "error": "Client does not hold the trajectory lock."}
        
    token = str(uuid.uuid4())
    token_key = f"room:{room}:upload_token:{token}"
    r.set(token_key, sid, ex=60)

    log.info(f"Issued upload token for room '{room}' to {sid}")
    return {"success": True, "token": token}

@socketio.on("frames:count")
def handle_len_frames(data):
    sid = request.sid
    room = get_project_room_from_session(sid)

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    try:
        indices_key = f"room:{room}:trajectory:indices"
        frame_count = r.zcard(indices_key)
        return {"success": True, "count": frame_count}
    except Exception as e:
        log.error(f"Failed to get frame count: {e}")
        return {"success": False, "error": "Failed to get frame count"}

if __name__ == '__main__':
    try:
        log.info("Starting ZnDraw Server")
        socketio.run(app, debug=True)
    finally:
        r.flushall()
        shutil.rmtree("data", ignore_errors=True)
