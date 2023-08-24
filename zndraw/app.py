import time
import uuid

from dask.distributed import Client
from flask import Flask, redirect, render_template, request, session, url_for
from flask_socketio import SocketIO, emit

from zndraw.data import DataHandler

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

socketio = SocketIO(app)

######################### Temporary #########################
import enum
from typing import Union

from ase.data import chemical_symbols
from pydantic import BaseModel, Field
from typing_extensions import Annotated
import ase
import numpy as np

Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})

class Duplicate(BaseModel):
    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)
    symbol: Symbols

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom
        return [atoms]

class Methods(BaseModel):
    method: Annotated[Union[Duplicate, None], Field(alias="Method")] = None

######################### Temporary #########################


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
        slice(max(_id - 50, 0), _id + 50)
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
    session["selection"] = data


@socketio.on("interaction:schema")
def get_interaction_scheme():
    schema = Methods.model_json_schema()
    # print(schema)

    emit("interaction:schema", schema)

@socketio.on("interaction:submit")
def run_interaction_scheme(data):
    print(data)
    data_handler = DataHandler(Client(app.config["dask-scheduler"]))

    atoms = Methods(**data).method.run(atom_ids=session["selection"], atoms=data_handler[session["step"]])
    data_handler[session["step"] + 1] = atoms[0]
    emit("cache:reset")
    emit("display", {"index": session["step"] + 1})

@app.route("/display/<int:index>")
def display(index):
    socketio.emit("display", {"index": index})
    return "OK"

@app.route("/cache/<int:index>")
def cache(index):
    socketio.emit("cache:load", {"index": index})
    return "OK"

@app.route("/cache/reset/")
def cache_reset():
    socketio.emit("cache:reset")
    return "OK"

@socketio.on("step")
def step(data):
    session["step"] = int(data["step"])