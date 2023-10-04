import importlib
import multiprocessing as mp
import uuid

from flask import Flask, render_template, request, session
from flask_socketio import SocketIO, emit, join_room

from zndraw.data import atoms_from_json
from zndraw.draw import Geometry
from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig
from zndraw.zndraw import FileIO, ZnDraw

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

io = SocketIO(
    app, max_http_buffer_size=int(1e10), cors_allowed_origins="*"
)  # , async_mode="threading")
# 10 GB Upload limit


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    if "uuid" not in session:
        session["uuid"] = str(uuid.uuid4())

        kwargs = {
            "url": request.url_root,
            "token": session["uuid"],
            "file": FileIO(
                name=app.config.get("filename"),
                start=app.config.get("start"),
                stop=app.config.get("stop"),
                step=app.config.get("step"),
            ),
            "wait": True,
        }

        proc = mp.Process(
            target=ZnDraw,
            kwargs=kwargs,
        )
        proc.start()

    return render_template(
        "index.html",
        upgrade_insecure_requests=app.config["upgrade_insecure_requests"],
        uuid=session["uuid"],
    )


@io.on("join")
def on_join(data):
    """Join a room."""
    try:
        join_room(session["uuid"])
    except KeyError:
        session["uuid"] = data["uuid"]
        join_room(session["uuid"])


@app.route("/exit")
def exit_route():
    """Exit the session."""
    print("Server shutting down...")
    io.stop()
    return "Server shutting down..."


@io.on("modifier:schema")
def modifier_schema():
    config = GlobalConfig.load()

    for modifier in config.modify_functions:
        module_name, function_name = modifier.rsplit(".", 1)
        try:
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
        except ImportError:
            print(f"Could not import {modifier}")


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
    io.emit("selection:schema", get_selection_class().model_json_schema())


@io.on("draw:schema")
def draw_schema():
    io.emit("draw:schema", Geometry.updated_schema())


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


@io.on("atoms:request")
def atoms_request(url):
    """Return the atoms."""
    emit("atoms:request", url, include_self=False, to=session["uuid"])


@io.on("modifier:run")
def modifier_run(data):
    emit("modifier:run", data, include_self=False, to=session["uuid"])


@io.on("selection:run")
def selection_run(data):
    emit("selection:run", data, include_self=False, to=session["uuid"])


@io.on("analysis:run")
def analysis_run(data):
    emit("analysis:run", data, include_self=False, to=session["uuid"])


@io.on("download:request")
def download_request(data):
    emit("download:request", data, include_self=False, to=session["uuid"])


@io.on("download:response")
def download_response(data):
    emit("download:response", data, include_self=False, to=session["uuid"])


@io.on("atoms:download")
def atoms_download(data):
    """Return the atoms."""
    emit("atoms:download", data, include_self=False, to=session["uuid"])


@io.on("atoms:upload")
def atoms_upload(data):
    """Return the atoms."""
    emit("atoms:upload", data, include_self=False, to=session["uuid"])


@io.on("view:set")
def display(data):
    """Display the atoms at the given index"""
    emit("view:set", data, include_self=False, to=session["uuid"])


@io.on("atoms:size")
def atoms_size(data):
    """Return the atoms."""
    emit("atoms:size", data, include_self=False, to=session["uuid"])


@io.on("upload")
def upload(data):
    emit("upload", data, include_self=False, to=session["uuid"])


@io.on("atoms:delete")
def delete_atoms(data):
    emit("atoms:delete", data, include_self=False, to=session["uuid"])


@io.on("atoms:insert")
def insert_atoms(data):
    emit("atoms:insert", data, include_self=False, to=session["uuid"])


@io.on("message:log")
def log(data):
    emit("message:log", data, include_self=False, to=session["uuid"])


@io.on("selection:set")
def selection_get(data):
    emit("selection:set", data, include_self=False, to=session["uuid"])


@io.on("draw:get_line")
def draw_points(data):
    emit("draw:get_line", data, include_self=False, to=session["uuid"])


@io.on("analysis:figure")
def analysis_figure(data):
    emit("analysis:figure", data, include_self=False, to=session["uuid"])
