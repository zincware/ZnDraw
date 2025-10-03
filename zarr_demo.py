import os

import numpy as np
import redis
import zarr
import znh5md
from flask import Flask, Response
from flask_socketio import SocketIO, emit

app = Flask(__name__)
app.config["SECRET_KEY"] = "secret!"
socketio = SocketIO(app)
r = redis.Redis()

ZARR_STORE_PATH = "structures.zarr"
# The original HDF5 file, only used for the one-time conversion.
H5_SOURCE_PATH = "/Users/fzills/tools/zndraw-communication-testing/structures.h5"


def create_zarr_store_from_h5(h5_path, zarr_path):
    """
    Reads data from an HDF5 file and writes it to a Zarr store.
    This is the "preprocessing" or "ingestion" step.
    """
    if os.path.exists(zarr_path):
        print(f"Zarr store '{zarr_path}' already exists. Skipping creation.")
        return

    print(f"Creating Zarr store from '{h5_path}'...")

    io = znh5md.IO(h5_path)
    total_frames = len(io)

    if total_frames == 0:
        print("Source file has no frames. Exiting.")
        return

    # Get the shape from the first frame
    positions = np.array([x.get_positions() for x in io[:]])

    z = zarr.create_array(
        store=ZARR_STORE_PATH,
        shape=positions.shape,
        chunks="auto",  # Chunk along the first axis (frames)
        dtype="f4",
    )

    print(f"Successfully created Zarr store with {total_frames} frames.")


@app.route("/frame/<int:frame_id>")
def get_frame(frame_id):
    try:
        # 1. Access the 'positions' array dataset inside the Zarr group.
        z = zarr.open(ZARR_STORE_PATH, mode="r")
        frame_data = z[frame_id]

        # 4. Return the binary data. The manual cache is no longer needed.
        return Response(
            np.array(frame_data).tobytes(), mimetype="application/octet-stream"
        )

    except KeyError:
        return Response(
            "Error: 'positions' dataset not found in Zarr store.", status=500
        )
    except Exception as e:
        return Response(f"An unexpected error occurred: {e}", status=500)


@socketio.on("set_frame")
def handle_set_frame(data):
    frame_id = data["frame"]

    r.set("current_frame", frame_id)
    r.publish("frame_channel", frame_id)
    emit("set_frame", {"frame": frame_id}, broadcast=True)


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

if __name__ == "__main__":
    create_zarr_store_from_h5(H5_SOURCE_PATH, ZARR_STORE_PATH)
    socketio.run(app)
