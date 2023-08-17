import uuid

import tqdm
import znh5md
from flask import Flask, redirect, render_template, session, url_for, request
from flask_socketio import SocketIO, emit
from dask.distributed import Client, Variable
from zndraw.data import DataHandler
import numpy as np
import time

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

socketio = SocketIO(app)


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    return render_template("index.html")


@app.route("/read/<filename>")
def read(filename):
    """Read a file."""
    filename = filename.replace(".slash.", "/").replace("\\", ".slash.")
    session["filename"] = filename
    # TODO convert file to H5MD, if not already

    # dataset = znh5md.ASEH5MD(filename)
    # for idx in tqdm.tqdm(range(1000)):
    #     _ = dataset.get_atoms_list(idx)

    return redirect(url_for("index"))


@app.route("/session")
def session_view():
    """View the session."""
    return dict(session)

@app.route("/exit")
def exit():
    """Exit the session."""
    socketio.stop()
    return "Server shutting down..."


@socketio.on("configuration:id")
def handle_get_id_on_configurations(json):
    start_time = time.time()

    _id = int(json["id"])

    print("received json: " + str(json) )
    atoms_dict = DataHandler(Client(app.config["dask-scheduler"])).get_atoms_json(slice(_id-50, _id+50))
    print("time to get atoms: " + str(time.time() - start_time))

    # emit("configuration:id", atoms_dict)
    return atoms_dict

@socketio.on("config")
def config():
    data_handler = DataHandler(Client(app.config["dask-scheduler"]))

    emit("config", {"n_frames": len(data_handler), "max_fps": 5})
