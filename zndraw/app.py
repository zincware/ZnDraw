# import eventlet

# eventlet.monkey_patch()

import importlib
import uuid

import numpy as np
import tqdm
from flask import Flask, render_template
from flask_socketio import SocketIO, emit

from zndraw.data import atoms_to_json, get_atomsdict_list
from zndraw.settings import GlobalConfig

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
    print("Server shutting down...")
    socketio.stop()
    return "Server shutting down..."


@socketio.on("atoms:request")
def atoms_request(data):
    """Return the atoms."""
    print(f"atoms:request {data = }")
    for atoms_dict in get_atomsdict_list(app.config["filename"]):
        emit("atoms:upload", atoms_dict)


@socketio.on("modifier:schema")
def modifier_schema():
    config = GlobalConfig.load()

    for modifier in config.modify_functions:
        module_name, function_name = modifier.rsplit(".", 1)
        module = importlib.import_module(module_name)
        modifier_cls = getattr(module, function_name)
        schema = modifier_cls.model_json_schema()

        if modifier in config.function_schema:
            kwargs = config.function_schema[modifier]
            for key, value in kwargs.items():
                schema["properties"][key]["default"] = value

        socketio.emit(
            "modifier:schema",
            {"name": modifier, "schema": schema},
        )


@socketio.on("modifier:run")
def modifier_run(data):
    import ase

    points = np.array([[val["x"], val["y"], val["z"]] for val in data["points"]])
    segments = np.array(data["segments"])


    atoms = ase.Atoms(
        numbers=data["atoms"]["numbers"],
        cell=data["atoms"]["cell"],
        pbc=True,
        positions=data["atoms"]["positions"],
    )

    module_name, function_name = data["name"].rsplit(".", 1)
    module = importlib.import_module(module_name)
    modifier_cls = getattr(module, function_name)
    modifier = modifier_cls(**data["params"])
    # available_methods = {x.__name__: x for x in [Explode, Duplicate]}

    # modifier = available_methods[data["name"]](**data["params"])
    print(f"modifier:run {modifier = } from {data['params'] = }")
    atoms_list = modifier.run(atom_ids=data["selection"], atoms=atoms, points=points, segments=segments)
    socketio.emit("atoms:clear", int(data["step"]) + 1)
    for idx, atoms in tqdm.tqdm(enumerate(atoms_list)):
        atoms_dict = atoms_to_json(atoms)
        socketio.emit("atoms:upload", {idx + 1 + int(data["step"]): atoms_dict})


@socketio.on("config")
def config(data):
    pass


# @app.route("/download/<int:idx>")
# def download(idx):
#     socketio.emit("atoms:download", [idx])
#     return "OK"

# @socketio.on("atoms:download")
# def atoms_download(data):
#     """Return the atoms."""
#     print(f"atoms:download {data = }")
