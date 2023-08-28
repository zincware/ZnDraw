import importlib
import uuid
from io import StringIO

import numpy as np
import tqdm
from flask import Flask, render_template
from flask_socketio import SocketIO, emit

from zndraw.data import atoms_from_json, atoms_to_json, get_atomsdict_list
from zndraw.settings import GlobalConfig

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

io = SocketIO(app, max_http_buffer_size=int(1e10))  # , async_mode="threading")
# 10 GB Upload limit


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    return render_template("index.html")


@app.route("/exit")
def exit_route():
    """Exit the session."""
    print("Server shutting down...")
    io.stop()
    return "Server shutting down..."


def _read_file(filename):
    for idx, atoms_dict in enumerate(get_atomsdict_list(filename)):
        io.emit("atoms:upload", atoms_dict)


@io.on("atoms:request")
def atoms_request(data):
    """Return the atoms."""

    if "filename" in app.config:
        io.start_background_task(target=_read_file, filename=app.config["filename"])
    else:
        emit("atoms:upload", {})


@io.on("modifier:schema")
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

        io.emit(
            "modifier:schema",
            {"name": modifier, "schema": schema},
        )


@io.on("modifier:run")
def modifier_run(data):
    import ase

    points = np.array([[val["x"], val["y"], val["z"]] for val in data["points"]])
    segments = np.array(data["segments"])

    if "atoms" in data:
        atoms = atoms_from_json(data["atoms"])
    else:
        atoms = ase.Atoms()

    module_name, function_name = data["name"].rsplit(".", 1)
    module = importlib.import_module(module_name)
    modifier_cls = getattr(module, function_name)
    modifier = modifier_cls(**data["params"])
    # available_methods = {x.__name__: x for x in [Explode, Duplicate]}

    # modifier = available_methods[data["name"]](**data["params"])
    print(f"modifier:run {modifier = }")
    atoms_list = modifier.run(
        atom_ids=data["selection"], atoms=atoms, points=points, segments=segments
    )
    io.emit("atoms:clear", int(data["step"]) + 1)
    for idx, atoms in tqdm.tqdm(enumerate(atoms_list)):
        atoms_dict = atoms_to_json(atoms)
        io.emit("atoms:upload", {idx + 1 + int(data["step"]): atoms_dict})

    io.emit("view:set", int(data["step"]) + 1)
    io.emit("view:play")


@io.on("analysis:schema")
def analysis_schema(data):
    config = GlobalConfig.load()
    if "atoms" not in data:
        return
    atoms = atoms_from_json(data["atoms"])

    for modifier in config.analysis_functions:
        module_name, function_name = modifier.rsplit(".", 1)
        module = importlib.import_module(module_name)
        modifier_cls = getattr(module, function_name)
        schema = modifier_cls.schema_from_atoms(atoms)

        io.emit(
            "analysis:schema",
            {"name": modifier, "schema": schema},
        )


@io.on("selection:schema")
def selection_schema():
    config = GlobalConfig.load()

    for selection in config.selection_functions:
        module_name, function_name = selection.rsplit(".", 1)
        module = importlib.import_module(module_name)
        cls = getattr(module, function_name)

        data = {"name": selection, "schema": cls.model_json_schema()}

        io.emit(
            "selection:schema",
            data,
        )


@io.on("selection:run")
def selection_run(data):
    import ase

    if "atoms" in data:
        atoms = atoms_from_json(data["atoms"])
    else:
        atoms = ase.Atoms()

    try:
        module_name, function_name = data["name"].rsplit(".", 1)
        module = importlib.import_module(module_name)
        selection_cls = getattr(module, function_name)
        selection = selection_cls(**data["params"])

        selected_ids = selection.get_ids(atoms, data["selection"])
        io.emit("selection:run", selected_ids)
    except ValueError as err:
        print(err)


@io.on("analysis:run")
def analysis_run(data):
    atoms_list = [atoms_from_json(x) for x in data["atoms_list"].values()]

    print(f"Analysing {len(atoms_list)} frames")

    module_name, function_name = data["name"].rsplit(".", 1)
    module = importlib.import_module(module_name)
    cls = getattr(module, function_name)
    instance = cls(**data["params"])

    fig = instance.run(atoms_list, data["selection"])
    return fig.to_json()


@io.on("config")
def config(data):
    pass


@io.on("download")
def download(data):
    atoms = [atoms_from_json(x) for x in data["atoms_list"].values()]
    if "selection" in data:
        atoms = [atoms[data["selection"]] for atoms in atoms]
    import ase.io

    file = StringIO()
    ase.io.write(file, atoms, format="extxyz")
    file.seek(0)
    return file.read()


@io.on("atoms:download")
def atoms_download(data):
    """Return the atoms."""
    emit("atoms:download", data, broadcast=True, include_self=False)


@io.on("atoms:upload")
def atoms_upload(data):
    """Return the atoms."""
    emit("atoms:upload", data, broadcast=True, include_self=False)


@io.on("view:set")
def display(data):
    """Display the atoms at the given index"""
    emit("view:set", data, broadcast=True, include_self=False)


@io.on("atoms:size")
def atoms_size(data):
    """Return the atoms."""
    emit("atoms:size", data, broadcast=True, include_self=False)


@io.on("upload")
def upload(data):
    from io import StringIO

    import ase.io
    import tqdm

    # tested with small files only

    format = data["filename"].split(".")[-1]
    if format == "h5":
        print("H5MD format not supported for uploading yet")
        # import znh5md
        # stream = BytesIO(data["content"].encode("utf-8"))
        # atoms = znh5md.ASEH5MD(stream).get_atoms_list()
        # for idx, atoms in tqdm.tqdm(enumerate(atoms)):
        #     atoms_dict = atoms_to_json(atoms)
        #     io.emit("atoms:upload", {idx: atoms_dict})
    else:
        stream = StringIO(data["content"])
        io.emit("atoms:clear", 0)
        for idx, atoms in tqdm.tqdm(enumerate(ase.io.iread(stream, format=format))):
            atoms_dict = atoms_to_json(atoms)
            io.emit("atoms:upload", {idx: atoms_dict})
    emit("view:set", 0)


@io.on("scene:schema")
def scene_schema():
    import enum

    from pydantic import BaseModel, Field

    class Material(str, enum.Enum):
        MeshBasicMaterial = "MeshBasicMaterial"
        MeshLambertMaterial = "MeshLambertMaterial"
        MeshMatcapMaterial = "MeshMatcapMaterial"
        MeshPhongMaterial = "MeshPhongMaterial"
        MeshPhysicalMaterial = "MeshPhysicalMaterial"
        MeshStandardMaterial = "MeshStandardMaterial"
        MeshToonMaterial = "MeshToonMaterial"

    # create a class for the material, resolution, etc.
    class Scene(BaseModel):
        material: Material = Field(Material.MeshPhongMaterial, description="Material")
        resolution: int = Field(10, ge=1, le=50, description="Resolution")
        wireframe: bool = Field(False, description="Wireframe")
        loop: bool = Field(
            False,
            alias="Animation Loop",
            description="Automatically restart animation when finished.",
        )
        simulation_box: bool = Field(
            False,
            description="Show the simulation box.",
        )
        bonds: bool = Field(
            True,
            description="Show bonds.",
        )

    schema = Scene.model_json_schema()

    schema["properties"]["wireframe"]["format"] = "checkbox"
    schema["properties"]["Animation Loop"]["format"] = "checkbox"
    schema["properties"]["simulation_box"]["format"] = "checkbox"
    schema["properties"]["resolution"]["format"] = "range"
    schema["properties"]["bonds"]["format"] = "checkbox"

    # import json

    # print(json.dumps(schema, indent=2))

    return schema
