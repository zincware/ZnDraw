import contextlib
import datetime
import logging
import traceback
from threading import Lock
from uuid import uuid4

from flask import current_app as app
from flask import request, session
from flask_socketio import call, emit, join_room

from ..app import socketio as io

log = logging.getLogger(__name__)

modifier_lock = Lock()


def _webclients_room(data: dict) -> str:
    """Return the room name for the webclients."""
    if "sid" in data:
        return data["sid"]
    return f"webclients_{data['token']}"


def _webclients_default(data: dict) -> str:
    """Return the SID of the default webclient."""
    if "sid" in data:
        return data["sid"]
    # TODO: if there is a keyerror, it will not be properly handled and the
    #  python interface is doomed to wait for TimeoutError.
    try:
        return f"webclients_{data['token']}"
    except KeyError:
        log.critical("No webclient connected.")


def _pyclients_room(data: dict) -> str:
    """All pyclients run via get, so this is not used."""
    return f"pyclients_{data['token']}"


def _pyclients_default(data: dict) -> str:
    """Return the SID of the default pyclient."""
    if "sid" in data:
        return data["sid"]
    return app.config["DEFAULT_PYCLIENT"]


def _get_uuid_for_sid(sid) -> str:
    """Given a sid, return the UUID that is associated with it.
    The SID is given by flask, the UUID is defined by zndraw
    and can be used to reconnect.
    """
    inv_clients = {v: k for k, v in app.config["pyclients"].items()}
    return inv_clients[sid]


def _get_queue_position(job_id) -> int:
    """Return the position of the job_id in the queue."""
    try:
        return app.config["MODIFIER"]["queue"].index(job_id)
    except ValueError:
        return -1


@io.on("connect")
def connect():
    if app.config["DEFAULT_PYCLIENT"] is None and "token" in session:
        # refuse connection if there is no default pyclient
        return False
    with contextlib.suppress(KeyError):
        # If you connect through Python, you don't have a token.
        token = session["token"]
        join_room(f"webclients_{token}")
        # who ever connected latest is the HOST of the room
        try:
            app.config["ROOM_HOSTS"][token].append(request.sid)
        except KeyError:
            app.config["ROOM_HOSTS"][token] = [request.sid]

        emit("webclient:available", request.sid, to=app.config["DEFAULT_PYCLIENT"])

        data = {"modifiers": []}  # {schema: ..., name: ...}
        for name, schema in app.config["MODIFIER"]["default_schema"].items():
            data["modifiers"].append({"schema": schema, "name": name})
        data["token"] = token

        emit("modifier:register", data, to=app.config["DEFAULT_PYCLIENT"])

        # TODO emit("modifier:register", _all modifiers_, to=app.config["DEFAULT_PYCLIENT"]')

        log.debug(
            f"connected {request.sid} and updated HOSTS to {app.config['ROOM_HOSTS']}"
        )
        emit("message:log", "Connection established", to=request.sid)
        if token not in app.config["PER-TOKEN-DATA"]:
            app.config["PER-TOKEN-DATA"][token] = {}

        # append to zndraw.log a line isoformat() + " " + token
        if "token" not in app.config:
            with open("zndraw.log", "a") as f:
                f.write(
                    datetime.datetime.now().isoformat() + " " + token + " connected \n"
                )


@io.on("disconnect")
def disconnect():
    with contextlib.suppress(KeyError):
        token = session["token"]
        try:
            del app.config["pyclients"][_get_uuid_for_sid(request.sid)]
        except KeyError:
            pass
        if "token" not in app.config:
            with open("zndraw.log", "a") as f:
                f.write(
                    datetime.datetime.now().isoformat()
                    + " "
                    + token
                    + " disconnected \n"
                )
        try:
            app.config["ROOM_HOSTS"][token].remove(request.sid)
        except ValueError:
            pass  # SID not in the list
        if not app.config["ROOM_HOSTS"][token]:
            del app.config["ROOM_HOSTS"][token]
        # remove the pyclient from the dict
    log.debug(
        f'disconnect {request.sid} and updated HOSTS to {app.config["ROOM_HOSTS"]}'
    )


@io.on("join")
def join(data: dict):
    """
    Arguments:
        data: {"token": str, "uuid": str}
    """
    # only used by pyclients that only connect via socket (no HTML)
    token = data["token"]
    uuid = data["uuid"]
    auth_token = data["auth_token"]
    session["authenticated"] = auth_token == app.config["AUTH_TOKEN"]
    log.debug(f"Client {request.sid} is {session['authenticated'] = }")
    if uuid in app.config["pyclients"]:
        log.critical(f"UUID {uuid} is already registered in {app.config['pyclients']}.")
        emit("message:log", f"UUID {uuid} is already registered", to=request.sid)
    app.config["pyclients"][uuid] = request.sid
    session["token"] = token
    join_room(f"pyclients_{token}")
    if token not in app.config["PER-TOKEN-DATA"]:
        app.config["PER-TOKEN-DATA"][token] = {}
    if token == "default":
        # this would be very easy to exploit
        app.config["DEFAULT_PYCLIENT"] = request.sid


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
        fps: int = Field(60, ge=1, le=120, description="Maxiumm frames per second")
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
        line_label: bool = Field(
            True,
            description="Show the length of the line.",
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
    schema["properties"]["line_label"]["format"] = "checkbox"

    return schema


@io.on("modifier:run")
def modifier_run(data):
    # emit entered the queue
    JOB_ID = uuid4()
    TIMEOUT = 60
    app.config["MODIFIER"]["queue"].append(JOB_ID)

    while True:
        if app.config["MODIFIER"]["queue"][0] == JOB_ID:
            acquired = modifier_lock.acquire(blocking=False)
        else:
            acquired = False
        if acquired:
            print("modifier_lock acquired")
            break
        else:
            emit(
                "modifier:run:enqueue",
                _get_queue_position(JOB_ID),
                to=_webclients_room({"token": session["token"]}),
            )
            print("waiting for modifier_lock")
            io.sleep(1)

    name = data["params"]["method"]["discriminator"]
    # handle custom modifiers (not default)
    if name in app.config["PER-TOKEN-DATA"][session["token"]].get("modifier", {}):
        token = app.config["PER-TOKEN-DATA"][session["token"]]["modifier"][name]
        data["sid"] = app.config["pyclients"][token]
    # handle custom modifiers (default)
    elif name in app.config["MODIFIER"]:
        data["sid"] = app.config["pyclients"][app.config["MODIFIER"][name]]

    # need to set the target of the modifier to the webclients room
    data["target"] = session["token"]
    # This should not go to request.sid but all webclients in the room
    emit("modifier:run:submitted", {}, to=_webclients_room({"token": session["token"]}))
    emit("modifier:run", data, to=_pyclients_default(data))

    io.sleep(TIMEOUT)
    if JOB_ID in app.config["MODIFIER"]["queue"]:
        # modifier failed
        app.config["MODIFIER"]["queue"].pop(0)
        app.config["MODIFIER"]["active"] = None
        modifier_lock.release()
        log.critical("Modifier failed - releasing lock.")
        # TODO: emit a error message that the modifier failed to the webclients
        emit(
            "modifier:run:criticalfail",
            {},
            to=_webclients_room({"token": session["token"]}),
        )


@io.on("analysis:run")
def analysis_run(data):
    data["target"] = session["token"]
    emit("analysis:run", data, include_self=False, to=_pyclients_default(data))


@io.on("analysis:figure")
def analysis_figure(data):
    emit(
        "analysis:figure", data["figure"], include_self=False, to=_webclients_room(data)
    )


@io.on("scene:set")
def scene_set(data):
    emit("scene:set", data["index"], include_self=False, to=_webclients_room(data))


@io.on("scene:step")
def scene_step(data):
    return call("scene:step", to=_webclients_room(data))


@io.on("atoms:download")
def atoms_download(data):
    return call("atoms:download", data["indices"], to=_webclients_default(data))


@io.on("atoms:upload")
def atoms_upload(data: dict):
    to = _webclients_default(data)
    # remove token and sid from the data, because JavaScript does not expect it
    data.pop("token", None)
    data.pop("sid", None)
    emit("atoms:upload", data, include_self=False, to=to)


@io.on("atoms:delete")
def atoms_delete(data: dict):
    emit("atoms:delete", data["index"], include_self=False, to=_webclients_room(data))


@io.on("atoms:length")
def atoms_length(data: dict):
    return call("atoms:length", to=_webclients_default(data))


@io.on("analysis:schema")
def analysis_schema(data: dict):
    if "sid" in data:
        emit("analysis:schema", data["schema"], include_self=False, to=data["sid"])
    else:
        raise ValueError


@io.on("modifier:schema")
def modifier_schema(data: dict):
    emit(
        "modifier:schema", data["schema"], include_self=False, to=_webclients_room(data)
    )


@io.on("selection:schema")
def selection_schema(data: dict):
    if "sid" in data:
        emit("selection:schema", data["schema"], include_self=False, to=data["sid"])
    else:
        raise ValueError


@io.on("draw:schema")
def draw_schema(data: dict):
    if "sid" in data:
        emit("draw:schema", data["schema"], include_self=False, to=data["sid"])
    else:
        raise ValueError


@io.on("points:get")
def scene_points(data: dict):
    return call("points:get", to=_webclients_default(data))


@io.on("scene:segments")
def scene_segments(data: dict):
    return call("scene:segments", to=_webclients_default(data))


@io.on("selection:get")
def selection_get(data: dict):
    return call("selection:get", to=_webclients_default(data))


@io.on("selection:set")
def selection_set(data: dict):
    emit(
        "selection:set",
        data["selection"],
        include_self=False,
        to=_webclients_room(data),
    )


@io.on("selection:run")
def selection_run(data: dict):
    data["target"] = session["token"]
    emit("selection:run", data, include_self=False, to=_pyclients_default(data))


@io.on("upload")
def upload(data):
    data = {"data": data, "target": session["token"]}
    emit("upload", data, include_self=False, to=_pyclients_default(data))


@io.on("atoms:insert")
def insert_atoms(data):
    raise NotImplementedError("This feature is not implemented yet.")


@io.on("message:log")
def message_log(data):
    emit("message:log", data["message"], to=_webclients_room(data))


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


@io.on("scene:play")
def scene_play(data):
    log.debug(f"scene:play {data}")
    emit("scene:play", to=_webclients_room(data))


@io.on("scene:pause")
def scene_pause(data):
    log.debug(f"scene:pause {data}")
    emit("scene:pause", to=_webclients_room(data))


@io.on("scene:trash")
def scene_trash(data):
    data["target"] = session["token"]
    emit("scene:trash", data, to=_pyclients_default(data))


@io.on("modifier:register")
def modifier_register(data):
    data["token"] = session["token"]
    if "modifier" not in app.config["PER-TOKEN-DATA"][session["token"]]:
        # this is a dict for the pyclient, but it has to be for every webclient ...
        app.config["PER-TOKEN-DATA"][session["token"]]["modifier"] = {}

    try:
        # we can only register one modifier at a time
        name = data["modifiers"][0]["name"]
        if name in app.config["MODIFIER"]["default_schema"]:
            msg = f"'{name}' is already registered as a default modifier and therefore reserved. Choose another name for your modifier!"
            log.critical(msg)
            emit("message:log", msg, to=request.sid)

        if name in app.config["PER-TOKEN-DATA"][session["token"]]["modifier"]:
            msg = f"Modifier {name} is already registered."
            log.critical(msg)
            emit("message:log", msg, to=request.sid)
        # get the key from the value request.sid by inverting the dict
        app.config["PER-TOKEN-DATA"][session["token"]]["modifier"][
            name
        ] = _get_uuid_for_sid(request.sid)
        log.critical(f'{app.config["PER-TOKEN-DATA"] = }')
        log.critical(f'{app.config["pyclients"] = }')

        if data["modifiers"][0]["default"]:
            if not session["authenticated"]:
                msg = "Unauthenticated users cannot register default modifiers."
                log.critical(msg)
                emit("message:log", msg, to=request.sid)
            app.config["MODIFIER"]["default_schema"][name] = data["modifiers"][0][
                "schema"
            ]
            app.config["MODIFIER"][name] = _get_uuid_for_sid(request.sid)
    except KeyError:
        print("Could not identify the modifier name.")
        traceback.print_exc()

    emit("modifier:register", data, to=app.config["DEFAULT_PYCLIENT"])


@io.on("bookmarks:get")
def bookmarks_get(data: dict):
    return call("bookmarks:get", to=_webclients_default(data))


@io.on("bookmarks:set")
def bookmarks_set(data: dict):
    emit(
        "bookmarks:set",
        data["bookmarks"],
        include_self=False,
        to=_webclients_room(data),
    )


@io.on("points:set")
def points_set(data: dict):
    emit("points:set", data["value"], include_self=False, to=_webclients_room(data))


@io.on("debug")
def debug(data: dict):
    emit("debug", data, include_self=False, to=_webclients_room(data))


@io.on("modifier:run:running")
def modifier_run_running(data: dict):
    app.config["MODIFIER"]["active"] = data.get("name", "unknown")
    emit("modifier:run:running", data, include_self=False, to=_webclients_room(data))


@io.on("modifier:run:finished")
def modifier_run_finished(data: dict):
    # remove 0th element from queue
    app.config["MODIFIER"]["queue"].pop(0)
    app.config["MODIFIER"]["active"] = None
    modifier_lock.release()
    print("modifier_lock released")
    emit("modifier:run:finished", data, include_self=False, to=_webclients_room(data))


@io.on("modifier:run:failed")
def modifier_run_failed():
    """Take care if the modifier does not respond."""
    # remove 0th element from queue
    app.config["MODIFIER"]["queue"].pop(0)

    app.config["MODIFIER"]["active"] = None
    modifier_lock.release()
    log.critical("Modifier failed - releasing lock.")
