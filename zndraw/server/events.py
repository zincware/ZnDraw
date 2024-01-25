import contextlib
import dataclasses
import datetime
import logging
import traceback
from threading import Lock
from uuid import uuid4

from flask import current_app as app
from flask import request, session
from flask_socketio import call, emit, join_room

from zndraw.utils import typecast

from ..app import cache
from ..app import socketio as io
from .data import (
    AnalysisFigureData,
    AnalysisRunData,
    AnalysisSchemaData,
    AtomsDownloadData,
    AtomsLengthData,
    BookmarksSetData,
    DeleteAtomsData,
    JoinData,
    MessageData,
    ModifierRegisterData,
    ModifierRunData,
    ModifierRunRunningData,
    ModifierSchemaData,
    PlayData,
    PointsSetData,
    SceneSetData,
    SceneStepData,
    SceneUpdateData,
    SelectionRunData,
    SelectionSetData,
    SubscribedUserData,
)

log = logging.getLogger(__name__)

modifier_lock = Lock()


def _webclients_room(data: dict) -> str:
    """Return the room name for the webclients."""
    if isinstance(data, dict):
        if "sid" in data:
            return data["sid"]
        return f"webclients_{data['token']}"
    elif hasattr(data, "sid"):
        if data.sid is not None:
            return data.sid
    return f"webclients_{data.token}"


def _webclients_default(data: dict) -> str:
    """Return the SID of the default webclient."""
    if isinstance(data, dict):
        if "sid" in data:
            return data["sid"]
        # TODO: if there is a keyerror, it will not be properly handled and the
        #  python interface is doomed to wait for TimeoutError.
        try:
            return f"webclients_{data['token']}"
        except KeyError:
            log.critical("No webclient connected.")
    else:
        if hasattr(data, "sid"):
            if data.sid is not None:
                return data.sid
        try:
            return f"webclients_{data.token}"
        except KeyError:
            log.critical("No webclient connected.")


def _pyclients_room(data: dict) -> str:
    """All pyclients run via get, so this is not used."""
    return f"pyclients_{data['token']}"


def _pyclients_default(data: dict) -> str:
    """Return the SID of the default pyclient."""
    if isinstance(data, dict):
        if "sid" in data:
            return data["sid"]
    elif hasattr(data, "sid"):
        if data.sid is not None:
            return data.sid
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


def _subscribe_user(data: SubscribedUserData, subscription_type: str):
    """
    Subscribe to user updates for a given subscription type.

    data: {user: str}
    subscription_type: str (e.g., "STEP" or "CAMERA")
    """
    token = session.get("token")
    if token is None:
        return

    cache_key = f"PER-TOKEN-{subscription_type}-SUBSCRIPTIONS:{token}"
    per_token_subscriptions = cache.get(cache_key) or {}

    names = cache.get(f"PER-TOKEN-NAME:{token}") or {}

    # Get the SID from data["user"] and add it to the list of subscribers
    for sid, name in names.items():
        if name == data.user:
            if (
                subscription_type == "CAMERA"
                and per_token_subscriptions.get(sid) == request.sid
            ):
                print("Cannot subscribe to a user that is subscribed to you")
                return

            per_token_subscriptions[request.sid] = sid
            break

    cache.set(cache_key, per_token_subscriptions)


def _get_subscribers(token: str, subscription_type: str):
    """
    Get subscribers for a given subscription type.

    token: str
    subscription_type: str (e.g., "STEP" or "CAMERA")
    """
    cache_key = f"PER-TOKEN-{subscription_type}-SUBSCRIPTIONS:{token}"
    per_token_subscriptions = cache.get(cache_key) or {}

    subscribers = [
        sid for sid, this in per_token_subscriptions.items() if this == request.sid
    ]

    return subscribers


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

        data = {"sid": request.sid, "token": token}
        data["host"] = app.config["ROOM_HOSTS"][token][0] == request.sid
        names = cache.get(f"PER-TOKEN-NAME:{session['token']}") or {}
        names[request.sid] = uuid4().hex[:8].upper()
        cache.set(f"PER-TOKEN-NAME:{session['token']}", names)

        emit("webclient:available", data, to=app.config["DEFAULT_PYCLIENT"])

        connected_users = [
            {"name": names[sid]} for sid in app.config["ROOM_HOSTS"][token]
        ]

        emit(
            "connectedUsers",
            list(reversed(connected_users)),
            to=_webclients_room({"token": token}),
        )

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
        names = cache.get(f"PER-TOKEN-NAME:{session['token']}") or {}
        connected_users = [
            {"name": names[sid]} for sid in app.config["ROOM_HOSTS"][token]
        ]
        emit(
            "connectedUsers",
            list(reversed(connected_users)),
            to=_webclients_room({"token": token}),
        )

    log.debug(
        f'disconnect {request.sid} and updated HOSTS to {app.config["ROOM_HOSTS"]}'
    )


@io.on("join")
@typecast
def join(data: JoinData):
    """
    Arguments:
        data: {"token": str, "uuid": str}
    """
    # only used by pyclients that only connect via socket (no HTML)
    session["authenticated"] = data.auth_token == app.config["AUTH_TOKEN"]
    log.debug(f"Client {request.sid} is {session['authenticated'] = }")
    if data.uuid in app.config["pyclients"]:
        log.critical(
            f"UUID {data.uuid} is already registered in {app.config['pyclients']}."
        )
        emit("message:log", f"UUID {data.uuid} is already registered", to=request.sid)
    app.config["pyclients"][data.uuid] = request.sid
    session["token"] = data.token
    join_room(f"pyclients_{data.token}")
    if data.token not in app.config["PER-TOKEN-DATA"]:
        app.config["PER-TOKEN-DATA"][data.token] = {}
    if data.token == "default":
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
@typecast
def modifier_run(data: ModifierRunData):
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

    name = data.name
    # handle custom modifiers (not default)
    if name in app.config["PER-TOKEN-DATA"][session["token"]].get("modifier", {}):
        token = app.config["PER-TOKEN-DATA"][session["token"]]["modifier"][name]
        data.sid = app.config["pyclients"][token]
    # handle custom modifiers (default)
    elif name in app.config["MODIFIER"]:
        data.sid = app.config["pyclients"][app.config["MODIFIER"][name]]

    # need to set the target of the modifier to the webclients room
    data.target = session["token"]
    # This should not go to request.sid but all webclients in the room
    emit("modifier:run:submitted", {}, to=_webclients_room({"token": session["token"]}))
    emit("modifier:run", dataclasses.asdict(data), to=_pyclients_default(data))

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
@typecast
def analysis_run(data: AnalysisRunData):
    data.target = session["token"]
    emit(
        "analysis:run",
        dataclasses.asdict(data),
        include_self=False,
        to=_pyclients_default(data),
    )


@io.on("analysis:figure")
@typecast
def analysis_figure(data: AnalysisFigureData):
    emit("analysis:figure", data.figure, include_self=False, to=_webclients_room(data))


@io.on("scene:set")
@typecast
def scene_set(data: SceneSetData):
    emit("scene:set", data.index, include_self=False, to=_webclients_room(data))


@io.on("scene:step")
@typecast
def scene_step(data: SceneStepData):
    return call("scene:step", to=_webclients_room(data))


@io.on("atoms:download")
@typecast
def atoms_download(data: AtomsDownloadData):
    return call("atoms:download", data.indices, to=_webclients_default(data))


@io.on("atoms:upload")
def atoms_upload(data: dict):
    # TODO: this has a bad structure and should be updates to fixed keys!
    to = _webclients_default(data)
    # remove token and sid from the data, because JavaScript does not expect it
    data.pop("token", None)
    data.pop("sid", None)
    emit("atoms:upload", data, include_self=False, to=to)


@io.on("atoms:delete")
@typecast
def atoms_delete(data: DeleteAtomsData):
    emit("atoms:delete", data.index, include_self=False, to=_webclients_room(data))


@io.on("atoms:length")
@typecast
def atoms_length(data: AtomsLengthData):
    return call("atoms:length", to=_webclients_default(data))


@io.on("analysis:schema")
@typecast
def analysis_schema(data: AnalysisSchemaData):
    emit("analysis:schema", data.schema, include_self=False, to=data.sid)


@io.on("modifier:schema")
@typecast
def modifier_schema(data: ModifierSchemaData):
    emit("modifier:schema", data.schema, include_self=False, to=_webclients_room(data))


@io.on("selection:schema")
@typecast
def selection_schema(data: AnalysisSchemaData):
    emit("selection:schema", data.schema, include_self=False, to=data.sid)


@io.on("draw:schema")
@typecast
def draw_schema(data: AnalysisSchemaData):
    emit("draw:schema", data.schema, include_self=False, to=data.sid)


@io.on("points:get")
@typecast
def scene_points(data: AtomsLengthData):
    return call("points:get", to=_webclients_default(data))


@io.on("scene:segments")
@typecast
def scene_segments(data: AtomsLengthData):
    return call("scene:segments", to=_webclients_default(data))


@io.on("selection:get")
@typecast
def selection_get(data: AtomsLengthData):
    return call("selection:get", to=_webclients_default(data))


@io.on("selection:set")
@typecast
def selection_set(data: SelectionSetData):
    if data.token is None:
        data.token = session["token"]
    emit(
        "selection:set",
        data.selection,
        include_self=False,
        to=_webclients_room(data),
    )


@io.on("selection:run")
@typecast
def selection_run(data: SelectionRunData):
    data.target = session["token"]
    emit(
        "selection:run",
        dataclasses.asdict(data),
        include_self=False,
        to=_pyclients_default(data),
    )


@io.on("upload")
def upload(data: dict):
    # TODO: this has a bad structure and should be updates to fixed keys!
    data = {"data": data, "target": session["token"]}
    emit("upload", data, include_self=False, to=_pyclients_default(data))


@io.on("atoms:insert")
def insert_atoms(data):
    raise NotImplementedError("This feature is not implemented yet.")


@io.on("message:log")
@typecast
def message_log(data: MessageData):
    emit("message:log", data.message, to=_webclients_room(data))


@io.on("download:request")
@typecast
def download_request(data: dict):
    # TODO: this has an optional key "selection"
    emit(
        "download:request",
        {"data": data, "sid": session["token"]},
        to=app.config["DEFAULT_PYCLIENT"],
    )


@io.on("download:response")
def download_response(data):
    emit("download:response", data["data"], to=data["sid"])


@io.on("scene:play")
@typecast
def scene_play(data: PlayData):
    log.debug(f"scene:play {data}")
    emit("scene:play", to=_webclients_room(data))


@io.on("scene:pause")
@typecast
def scene_pause(data: PlayData):
    log.debug(f"scene:pause {data}")
    emit("scene:pause", to=_webclients_room(data))


@io.on("scene:trash")
def scene_trash():
    data = {"target": session["token"]}
    emit("scene:trash", data, to=_pyclients_default(data))


@io.on("modifier:register")
@typecast
def modifier_register(data: ModifierRegisterData):
    data.token = session["token"]
    if "modifier" not in app.config["PER-TOKEN-DATA"][session["token"]]:
        # this is a dict for the pyclient, but it has to be for every webclient ...
        app.config["PER-TOKEN-DATA"][session["token"]]["modifier"] = {}

    try:
        # we can only register one modifier at a time
        name = data.name
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

        if data.is_default:
            if not session["authenticated"]:
                msg = "Unauthenticated users cannot register default modifiers."
                log.critical(msg)
                emit("message:log", msg, to=request.sid)
            app.config["MODIFIER"]["default_schema"][name] = data.schema
            app.config["MODIFIER"][name] = _get_uuid_for_sid(request.sid)
    except KeyError:
        print("Could not identify the modifier name.")
        traceback.print_exc()

    emit(
        "modifier:register", dataclasses.asdict(data), to=app.config["DEFAULT_PYCLIENT"]
    )


@io.on("bookmarks:get")
@typecast
def bookmarks_get(data: AtomsLengthData):
    return call("bookmarks:get", to=_webclients_default(data))


@io.on("bookmarks:set")
@typecast
def bookmarks_set(data: BookmarksSetData):
    emit(
        "bookmarks:set",
        data.bookmarks,
        include_self=False,
        to=_webclients_room(data),
    )


@io.on("points:set")
@typecast
def points_set(data: PointsSetData):
    if data.token is None:
        data.token = session["token"]
    emit("points:set", data.value, include_self=False, to=_webclients_room(data))


@io.on("debug")
def debug(data: dict):
    emit("debug", data, include_self=False, to=_webclients_room(data))


@io.on("modifier:run:running")
@typecast
def modifier_run_running(data: ModifierRunRunningData):
    app.config["MODIFIER"]["active"] = data.name
    emit(
        "modifier:run:running",
        dataclasses.asdict(data),
        include_self=False,
        to=_webclients_room(data),
    )


@io.on("modifier:run:finished")
@typecast
def modifier_run_finished(data: AtomsLengthData):
    # remove 0th element from queue
    app.config["MODIFIER"]["queue"].pop(0)
    app.config["MODIFIER"]["active"] = None
    modifier_lock.release()
    print("modifier_lock released")
    emit(
        "modifier:run:finished",
        dataclasses.asdict(data),
        include_self=False,
        to=_webclients_room(data),
    )


@io.on("modifier:run:failed")
def modifier_run_failed():
    """Take care if the modifier does not respond."""
    # remove 0th element from queue
    app.config["MODIFIER"]["queue"].pop(0)

    app.config["MODIFIER"]["active"] = None
    modifier_lock.release()
    log.critical("Modifier failed - releasing lock.")


@io.on("connectedUsers:subscribe:step")
@typecast
def connected_users_subscribe_step(data: SubscribedUserData):
    """
    Subscribe to step updates for connected users.

    data: {user: str}
    """
    _subscribe_user(data, "STEP")


@io.on("connectedUsers:subscribe:camera")
@typecast
def connected_users_subscribe_camera(data: SubscribedUserData):
    """
    Subscribe to camera updates for connected users.

    data: {user: str}
    """
    _subscribe_user(data, "CAMERA")


@io.on("scene:update")
@typecast
def scene_update(data: SceneUpdateData):
    """Update the scene.

    data: {step: int, camera: {position: [float, float, float], rotation: [float, float, float]}}
    """
    token = session.get("token")
    if token is None:
        return

    if data.step is not None:
        step_subscribers = _get_subscribers(token, "STEP")
        emit(
            "scene:update",
            dataclasses.asdict(data),
            include_self=False,
            to=step_subscribers,
        )

    if data.camera is not None:
        camera_subscribers = _get_subscribers(token, "CAMERA")
        emit(
            "scene:update",
            dataclasses.asdict(data),
            include_self=False,
            to=camera_subscribers,
        )
