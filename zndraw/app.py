import multiprocessing as mp
import uuid
from io import StringIO

import ase
import numpy as np
import tqdm
from flask import Flask, render_template
from flask_socketio import SocketIO, emit

from zndraw.analyse import get_analysis_class
from zndraw.data import atoms_from_json, atoms_to_json
from zndraw.draw import Geometry
from zndraw.modify import get_modify_class
from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig
from zndraw.zndraw import ZnDraw

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

io = SocketIO(
    app, max_http_buffer_size=int(1e10), cors_allowed_origins="*"
)  # , async_mode="threading")
# 10 GB Upload limit


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    if "upgrade_insecure_requests" in app.config:
        return render_template(
            "index.html",
            upgrade_insecure_requests=app.config["upgrade_insecure_requests"],
        )
    return render_template("index.html")


@app.route("/exit")
def exit_route():
    """Exit the session."""
    print("Server shutting down...")
    io.stop()
    return "Server shutting down..."


def _read_file(filename, start, stop, step, compute_bonds, url=None):
    if url is None:
        if compute_bonds:
            instance = ZnDraw(socket=io, display_new=False)
        else:
            instance = ZnDraw(socket=io, display_new=False, bonds_calculator=None)
    else:
        if compute_bonds:
            instance = ZnDraw(url=url, display_new=False)
        else:
            instance = ZnDraw(url=url, display_new=False, bonds_calculator=None)

    instance.read(filename, start, stop, step)


@io.on("atoms:request")
def atoms_request(url):
    """Return the atoms."""

    if "filename" in app.config:
        if app.config["multiprocessing"]:
            proc = mp.Process(
                target=_read_file,
                args=(
                    app.config["filename"],
                    app.config["start"],
                    app.config["stop"],
                    app.config["step"],
                    app.config["compute_bonds"],
                    url,
                ),
            )
            proc.start()
        else:
            io.start_background_task(
                target=_read_file,
                filename=app.config["filename"],
                start=app.config["start"],
                stop=app.config["stop"],
                step=app.config["step"],
                compute_bonds=app.config["compute_bonds"],
            )
    else:
        emit("atoms:upload", {0: atoms_to_json(ase.Atoms())})


@io.on("modifier:schema")
def modifier_schema():
    config = GlobalConfig.load()

    cls = get_modify_class(config.get_modify_methods())
    io.emit("modifier:schema", cls.model_json_schema())


@io.on("modifier:run")
def modifier_run(data):
    config = GlobalConfig.load()
    cls = get_modify_class(config.get_modify_methods())

    points = np.array([[val["x"], val["y"], val["z"]] for val in data["points"]])
    segments = np.array(data["segments"])

    if "atoms" in data:
        atoms = atoms_from_json(data["atoms"])
    else:
        atoms = ase.Atoms()
    instance = cls(**data["params"])

    atoms_list = instance.run(
        atom_ids=data["selection"],
        atoms=atoms,
        points=points,
        segments=segments,
        json_data=data["atoms"] if "atoms" in data else None,
        url=data["url"],
    )
    io.emit("atoms:clear", int(data["step"]) + 1)
    for idx, atoms in tqdm.tqdm(enumerate(atoms_list)):
        atoms_dict = atoms_to_json(atoms)
        io.emit("atoms:upload", {idx + 1 + int(data["step"]): atoms_dict})

    io.emit("view:set", int(data["step"]) + 1)
    io.emit("view:play")


@io.on("selection:schema")
def selection_schema():
    config = GlobalConfig.load()
    cls = get_selection_class(config.get_selection_methods())

    io.emit("selection:schema", cls.model_json_schema())


@io.on("selection:run")
def selection_run(data):
    config = GlobalConfig.load()
    cls = get_selection_class(config.get_selection_methods())

    if "atoms" in data:
        atoms = atoms_from_json(data["atoms"])
    else:
        atoms = ase.Atoms()

    try:
        selection = cls(**data["params"])
        selected_ids = selection.get_ids(atoms, data["selection"])
        io.emit("selection:run", selected_ids)
    except ValueError as err:
        print(err)


@io.on("analysis:schema")
def analysis_schema(data):
    config = GlobalConfig.load()

    cls = get_analysis_class(config.get_analysis_methods())

    if "atoms" not in data:
        return
    atoms = atoms_from_json(data["atoms"])

    io.emit("analysis:schema", cls.model_json_schema_from_atoms(atoms))


@io.on("analysis:run")
def analysis_run(data):
    config = GlobalConfig.load()
    cls = get_analysis_class(config.get_analysis_methods())

    atoms_list = [atoms_from_json(x) for x in data["atoms_list"].values()]

    print(f"Analysing {len(atoms_list)} frames")

    try:
        instance = cls(**data["params"])
        fig = instance.run(atoms_list, data["selection"])
        return fig.to_json()
    except ValueError as err:
        print(err)


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
        fps: int = Field(60, ge=1, le=120, description="Frames per second")
        material: Material = Field(Material.MeshPhongMaterial, description="Material")
        resolution: int = Field(10, ge=1, le=50, description="Resolution")
        particle_size: float = Field(1.0, ge=0.1, le=5, description="Particle Size")
        bonds_size: float = Field(1.0, ge=0.1, le=5, description="Bonds Size")
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
        label_offset: int = Field(
            0,
            ge=-7,
            le=7,
            description="Move the label to the left or right (keypress i).",
        )

    schema = Scene.model_json_schema()

    schema["properties"]["wireframe"]["format"] = "checkbox"
    schema["properties"]["Animation Loop"]["format"] = "checkbox"
    schema["properties"]["simulation_box"]["format"] = "checkbox"
    schema["properties"]["resolution"]["format"] = "range"
    schema["properties"]["label_offset"]["format"] = "range"
    schema["properties"]["particle_size"]["format"] = "range"
    schema["properties"]["fps"]["format"] = "range"
    schema["properties"]["particle_size"]["step"] = 0.1
    schema["properties"]["bonds_size"]["format"] = "range"
    schema["properties"]["bonds_size"]["step"] = 0.1
    schema["properties"]["bonds"]["format"] = "checkbox"

    # import json

    # print(json.dumps(schema, indent=2))

    return schema


@io.on("draw:schema")
def draw_schema():
    io.emit("draw:schema", Geometry.updated_schema())


@io.on("atoms:delete")
def delete_atoms(data):
    emit("atoms:delete", data, broadcast=True, include_self=False)


@io.on("atoms:insert")
def insert_atoms(data):
    emit("atoms:insert", data, broadcast=True, include_self=False)


@io.on("message:log")
def log(data):
    emit("message:log", data, broadcast=True, include_self=False)
