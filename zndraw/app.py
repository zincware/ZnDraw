import contextlib
import uuid

from flask import Flask, redirect, render_template, request, session
from flask_socketio import SocketIO, call, emit, join_room

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())
app.config["ROOM_HOSTS"] = {}
app.config["DEFAULT_PYCLIENT"] = None

io = SocketIO(
    app, max_http_buffer_size=int(1e10), cors_allowed_origins="*"
)  # , async_mode="threading")
# 10 GB Upload limit


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    try:
        token = session["token"]
    except KeyError:
        token = uuid.uuid4().hex
        session["token"] = token

    return render_template(
        "index.html",
        upgrade_insecure_requests=app.config["upgrade_insecure_requests"],
        token=session["token"],
    )


@io.on("connect")
def connect():
    with contextlib.suppress(KeyError):
        # If you connect through Python, you don't have a token.

        token = session["token"]
        join_room(token)
        join_room("webclients")
        # who ever connected latest is the HOST of the room
        try:
            app.config["ROOM_HOSTS"][token].append(request.sid)
        except KeyError:
            app.config["ROOM_HOSTS"][token] = [request.sid]

        emit("webclient:available", request.sid, to="default")

        print(
            f"connected {request.sid} and updated HOSTS to {app.config['ROOM_HOSTS']}"
        )


@io.on("disconnect")
def disconnect():
    with contextlib.suppress(KeyError):
        token = session["token"]
        # leave_room(token) # I guess if disconnect, it will automatically leave the room?
        # leave_room("webclients") # wrap in try..except?
        app.config["ROOM_HOSTS"][token].remove(request.sid)
        if not app.config["ROOM_HOSTS"][token]:
            del app.config["ROOM_HOSTS"][token]
    print(f'disconnect {request.sid} and updated HOSTS to {app.config["ROOM_HOSTS"]}')


@app.route("/token/<token>")
def token(token):
    session["token"] = token
    return redirect("/")


@io.on("join")
def join(token):
    # only used by pyclients that only connect via socket (no HTML)
    session["token"] = token
    print(f"pyclient {request.sid} joined room {token}")
    join_room(token)
    join_room("pyclients")
    if token == "default":
        app.config["DEFAULT_PYCLIENT"] = request.sid


@app.route("/exit")
def exit_route():
    """Exit the session."""
    print("Server shutting down...")
    io.stop()
    return "Server shutting down..."


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


@io.on("modifier:run")
def modifier_run(data):
    if "sid" in data:
        sid = data.pop("sid")
        data["sid"] = request.sid
        return call("modifier:run", data, to=sid)
    else:
        data["sid"] = request.sid
        return call("modifier:run", data, to=app.config["DEFAULT_PYCLIENT"])


@io.on("analysis:run")
def analysis_run(data):
    if "sid" in data:
        sid = data.pop("sid")
        data["sid"] = request.sid
        emit("analysis:run", data, to=sid)
    else:
        data["sid"] = request.sid
        emit("analysis:run", data, to=app.config["DEFAULT_PYCLIENT"])


@io.on("analysis:figure")
def analysis_figure(data):
    if "sid" in data:
        sid = data.pop("sid")
        emit("analysis:figure", data["figure"], to=sid)
    else:
        emit("analysis:figure", data["figure"], to=session["token"])


@io.on("scene:set")
def scene_set(data):
    if "sid" in data:
        emit("scene:set", data["index"], include_self=False, to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            emit("scene:set", data["index"], include_self=False, to=session["token"])
        except KeyError:
            return "No host found."


@io.on("scene:step")
def scene_step(data):
    if "sid" in data:
        return call("scene:step", to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            return call("scene:step", to=app.config["ROOM_HOSTS"][session["token"]][0])
        except KeyError:
            return "No host found."


@io.on("atoms:download")
def atoms_download(data):
    if "sid" in data:
        sid = data.pop("sid")
        return call("atoms:download", data["indices"], to=sid)
    else:
        try:
            return call(
                "atoms:download",
                data["indices"],
                to=app.config["ROOM_HOSTS"][session["token"]][0],
            )
        except KeyError:
            return "No host found."


@io.on("atoms:upload")
def atoms_upload(data: dict):
    if "sid" in data:
        # if the data is sent from the default pyclient, it will have a sid
        sid = data.pop("sid")
        emit("atoms:upload", data, include_self=False, to=sid)
    else:
        emit("atoms:upload", data, include_self=False, to=session["token"])


@io.on("atoms:delete")
def atoms_delete(data: dict):
    if "sid" in data:
        # if the data is sent from the default pyclient, it will have a sid
        sid = data.pop("sid")
        emit("atoms:delete", data["index"], include_self=False, to=sid)
    else:
        emit("atoms:delete", data["index"], include_self=False, to=session["token"])


@io.on("atoms:length")
def atoms_length(data: dict):
    print(f"atoms:length for {data}")
    if "sid" in data:
        return call("atoms:length", to=data["sid"])
    else:
        try:
            return call(
                "atoms:length", to=app.config["ROOM_HOSTS"][session["token"]][0]
            )
        except KeyError:
            return "No host found."


@io.on("analysis:schema")
def analysis_schema(data: dict):
    if "sid" in data:
        emit("analysis:schema", data["schema"], include_self=False, to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            emit(
                "analysis:schema",
                data["schema"],
                include_self=False,
                to=session["token"],
            )
        except KeyError:
            return "No host found."


@io.on("modifier:schema")
def modifier_schema(data: dict):
    if "sid" in data:
        emit("modifier:schema", data["schema"], include_self=False, to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            emit(
                "modifier:schema",
                data["schema"],
                include_self=False,
                to=session["token"],
            )
        except KeyError:
            return "No host found."


@io.on("selection:schema")
def selection_schema(data: dict):
    if "sid" in data:
        emit("selection:schema", data["schema"], include_self=False, to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            emit(
                "selection:schema",
                data["schema"],
                include_self=False,
                to=session["token"],
            )
        except KeyError:
            return "No host found."


@io.on("draw:schema")
def draw_schema(data: dict):
    if "sid" in data:
        emit("draw:schema", data["schema"], include_self=False, to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            emit("draw:schema", data["schema"], include_self=False, to=session["token"])
        except KeyError:
            return "No host found."


@io.on("scene:points")
def scene_points(data: dict):
    if "sid" in data:
        return call("scene:points", to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            return call(
                "scene:points", to=app.config["ROOM_HOSTS"][session["token"]][0]
            )
        except KeyError:
            return "No host found."


@io.on("scene:segments")
def scene_segments(data: dict):
    if "sid" in data:
        return call("scene:segments", to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            return call(
                "scene:segments", to=app.config["ROOM_HOSTS"][session["token"]][0]
            )
        except KeyError:
            return "No host found."


@io.on("selection:get")
def selection_get(data: dict):
    if "sid" in data:
        return call("selection:get", to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            return call(
                "selection:get", to=app.config["ROOM_HOSTS"][session["token"]][0]
            )
        except KeyError:
            return "No host found."


@io.on("selection:set")
def selection_set(data: dict):
    if "sid" in data:
        emit("selection:set", data["selection"], include_self=False, to=data["sid"])
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            emit(
                "selection:set",
                data["selection"],
                include_self=False,
                to=session["token"],
            )
        except KeyError:
            return "No host found."


@io.on("selection:run")
def selection_run(data: dict):
    sid = data.pop("sid", None)
    data["sid"] = request.sid
    if sid is not None:
        emit("selection:run", data, to=sid)
    else:
        try:
            # emit to all webclients in the group, if no sid is provided
            emit("selection:run", data, to=app.config["DEFAULT_PYCLIENT"])
        except KeyError:
            return "No host found."


@io.on("upload")
def upload(data):
    if "sid" in data:  # currently not expected to happen
        raise ValueError("Uploading to specific pyclient currently not supported.")
    else:
        emit(
            "upload",
            {"data": data, "sid": session["token"]},
            to=app.config["DEFAULT_PYCLIENT"],
        )


@io.on("atoms:insert")
def insert_atoms(data):
    raise NotImplementedError("This feature is not implemented yet.")


@io.on("message:log")
def log(data):
    sid = data.pop("sid", session["token"])
    emit("message:log", data["message"], to=sid)


@io.on("download:request")
def download_request(data):
    emit(
        "download:request",
        {"data": data, "sid": session["token"]},
        to=app.config["DEFAULT_PYCLIENT"],
    )


@io.on("download:response")
def download_response(data):
    emit("download:response", data["data"], to=data["sid"])
