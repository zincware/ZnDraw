import time
import uuid


from dask.distributed import Client, Variable
from flask import Flask, redirect, render_template, request, session, url_for
from flask_socketio import SocketIO, emit

from zndraw.data import DataHandler

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

socketio = SocketIO(app)


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    return render_template("index.html")


@app.route("/exit")
def exit():
    """Exit the session."""
    socketio.stop()
    return "Server shutting down..."


@socketio.on("configuration:id")
def handle_get_id_on_configurations(json):
    start_time = time.time()

    _id = int(json["id"])

    print("received json: " + str(json))
    atoms_dict = DataHandler(Client(app.config["dask-scheduler"])).get_atoms_json(
        slice(_id - 50, _id + 50)
    )
    print("time to get atoms: " + str(time.time() - start_time))

    # emit("configuration:id", atoms_dict)
    return atoms_dict


@socketio.on("config")
def config():
    data_handler = DataHandler(Client(app.config["dask-scheduler"]))

    emit("config", {"n_frames": len(data_handler), "max_fps": 5})

@socketio.on("selection")
def selection(data):
    print(data)