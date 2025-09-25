from flask import Flask, Response, request
from flask_socketio import SocketIO, emit, join_room, rooms, leave_room
import redis
import znh5md
import numpy as np
import zarr
import os
import typing as t
import uuid
import logging
from collections import defaultdict

log = logging.getLogger(__name__)
# attch handler to print the time
log.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)

class LockData(t.TypedDict):
    target: str

app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'
socketio = SocketIO(app)
r = redis.Redis()
locks = {}

upload_key: tuple[str, str]| None = None
room_data = defaultdict(list)

def get_project_room_from_session(sid):
    """
    Finds the project room a client has joined.
    Filters out the client's personal session ID room.
    """
    for room in rooms(sid=sid):
        if room != sid:
            return room
    return None


@app.route("/frame/<string:room_id>/<int:frame_id>")
def get_frame(room_id, frame_id):
    if room_id not in room_data:
        return Response("Room not found", status=404)
    frames = room_data[room_id]
    if frame_id < 0 or frame_id >= len(frames):
        return Response("Frame not found", status=404)
    frame = frames[frame_id]
    return Response(frame.tobytes(), content_type='application/octet-stream')


# @socketio.on('disconnect')
# def handle_disconnect():
#     sid = request.sid
#     print(f"Client disconnected: {sid}")
#     lock_keys = r.scan_iter("room:*:lock:*")
#     for key in lock_keys:
#         if r.get(key) == sid:
#             print(f"Cleaning up orphaned lock '{key}' held by {sid}")
#             r.delete(key)

@socketio.on("join_room")
def handle_join(data):
    room = data["room"]
    if previous_room := get_project_room_from_session(request.sid) is not None:
        leave_room(previous_room)
        print(f"Client left room: {previous_room}")
        
    join_room(room)
    print(f"Client joined room: {room}")


@socketio.on("lock:acquire")
def acquire_lock(data: LockData):
    sid = request.sid
    room = get_project_room_from_session(sid)

    if locks.get(room) is None:
        locks[room] = sid
        log.info(f"Lock acquired for room {room} by {sid}")
        return True
    log.info(f"Lock for room {room} already held by {locks.get(room)}, denied for {sid}")
    return False

@socketio.on("lock:release")
def release_lock(data: LockData):
    sid = request.sid
    room = get_project_room_from_session(sid)
    if locks.get(room) == sid:
        locks[room] = None
        log.info(f"Lock released for room {room} by {sid}")
        return True
    else:
        log.info(f"Lock for room {room} held by {locks.get(room)}, not by {sid}")
    return False

@socketio.on("upload:prepare")
def handle_upload_prepare(data):
    """
    Generates a single-use token for a client to authorize an HTTP upload.
    """
    sid = request.sid
    room = get_project_room_from_session(sid)

    if not room:
        return {"success": False, "error": "Client has not joined a room."}

    if locks.get(room) != sid:
        return {"success": False, "error": "Client does not hold the trajectory lock."}
        
    token = str(uuid.uuid4())
    global upload_key
    upload_key = (token, room)

    # Store the token with a short TTL, associated with the client's SID
    # r.set(token_key, sid, ex=60) # Token is valid for 60 seconds

    log.info(f"Issued upload token for room '{room}' to {sid}")
    return {"success": True, "token": token}


@app.route("/upload", methods=["POST"])
def append_frame():
    auth_header = request.headers.get("Authorization")
    if not auth_header or not auth_header.startswith("Bearer "):
        return {"error": "Authorization token is missing or invalid"}, 401
    
    token = auth_header.split(" ")[1]
    global upload_key
    if upload_key is None:
        return {"error": "No upload token available"}, 403
    _token, _room = upload_key
    if token != _token:
        return {"error": "Invalid or expired token"}, 403
    upload_key = None  # Invalidate the token after use

    new_frame_data = np.frombuffer(request.data, dtype=np.float32)
    room_data[_room].append(new_frame_data)
    log.info(f"Received frame data of size {new_frame_data.nbytes} bytes")
    # np.save("received_frame.npy", new_frame_data)
    return {"success": True}


# @socketio.on("set_frame")
# def handle_set_frame(data):
#     frame_id = data["frame"]

#     r.set("current_frame", frame_id)
#     r.publish("frame_channel", frame_id)
#     emit("set_frame", {"frame": frame_id}, broadcast=True)

# # optional for multiple server instances
# def redis_listener():
#     pubsub = r.pubsub()
#     pubsub.subscribe("frame_channel")
#     for msg in pubsub.listen():
#         if msg["type"] == "message":
#             frame_id = int(msg["data"])
#             # Rebroadcast via socket.io
#             emit("set_frame", {"frame": frame_id}, broadcast=True)

# socketio.start_background_task(redis_listener)

if __name__ == '__main__':
    # create_zarr_store_from_h5(H5_SOURCE_PATH, ZARR_STORE_PATH)
    socketio.run(app)
