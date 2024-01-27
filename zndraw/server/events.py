import dataclasses
import datetime
import logging
from threading import Lock
from uuid import uuid4

from celery import chain
from flask import request, session
from flask_socketio import call, emit, join_room

from zndraw.server import tasks
from zndraw.utils import typecast

from ..app import cache
from ..app import socketio as io
from ..data import (
    AnalysisFigureData,
    AtomsDownloadData,
    CeleryTaskData,
    DeleteAtomsData,
    JoinData,
    MessageData,
    SceneSetData,
    SceneUpdateData,
    SchemaData,
    SubscribedUserData,
)

log = logging.getLogger(__name__)

modifier_lock = Lock()


def get_main_room_host(token: str) -> str:
    ROOM_HOSTS = cache.get("ROOM_HOSTS")
    return ROOM_HOSTS[token][0]


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


def _get_uuid_for_sid(sid) -> str:
    """Given a sid, return the UUID that is associated with it.
    The SID is given by flask, the UUID is defined by zndraw
    and can be used to reconnect.
    """
    inv_clients = {v: k for k, v in cache.get("pyclients").items()}
    return inv_clients[sid]


def _get_queue_position(job_id) -> int:
    """Return the position of the job_id in the queue."""
    try:
        return cache.get("MODIFIER")["queue"].index(job_id)
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
    try:
        token = session["token"]
        # if you connect through Python, you don't have a token
        read_file_chain = chain(
            tasks.read_file.s(request.url_root, request.sid),
            tasks.analysis_schema.si(request.url_root, token),
        )
        read_file_chain.delay()
        tasks.get_selection_schema.delay(request.url_root, request.sid)
        tasks.scene_schema.delay(request.url_root, request.sid)
        tasks.geometries_schema.delay(request.url_root, request.sid)
        tasks.modifier_schema.delay(request.url_root, request.sid)

        join_room(f"webclients_{token}")
        # who ever connected latest is the HOST of the room
        ROOM_HOSTS = cache.get("ROOM_HOSTS")

        if token not in ROOM_HOSTS:
            ROOM_HOSTS[token] = [request.sid]
        else:
            ROOM_HOSTS[token].append(request.sid)

        cache.set("ROOM_HOSTS", ROOM_HOSTS)

        data = {"sid": request.sid, "token": token}
        data["host"] = ROOM_HOSTS[token][0] == request.sid
        names = cache.get(f"PER-TOKEN-NAME:{session['token']}") or {}
        names[request.sid] = uuid4().hex[:8].upper()
        cache.set(f"PER-TOKEN-NAME:{session['token']}", names)

        connected_users = [{"name": names[sid]} for sid in ROOM_HOSTS[token]]

        emit(
            "connectedUsers",
            list(reversed(connected_users)),
            to=_webclients_room({"token": token}),
        )

        # TODO: modifier registry
        # data = {"modifiers": []}  # {schema: ..., name: ...}
        # MODIFIER = cache.get("MODIFIER")
        # for name, schema in MODIFIER["default_schema"].items():
        #     data["modifiers"].append({"schema": schema, "name": name})
        # data["token"] = token

        # emit("modifier:register", data, to=DEFAULT_PYCLIENT)

        # TODO emit("modifier:register", _all modifiers_, to=app.config["DEFAULT_PYCLIENT"]')

        log.debug(f"connected {request.sid} and updated HOSTS to {ROOM_HOSTS}")
        emit("message:log", "Connection established", to=request.sid)
        PER_TOKEN_DATA = cache.get("PER-TOKEN-DATA")
        if token not in PER_TOKEN_DATA:
            PER_TOKEN_DATA[token] = {}
        cache.set("PER-TOKEN-DATA", PER_TOKEN_DATA)

        # append to zndraw.log a line isoformat() + " " + token
        log.info(
            datetime.datetime.now().isoformat() + " " + token
            if token
            else "client" + " connected"
        )

    except KeyError:
        # clients that connect directly via socketio do not call "join" to
        # register their token
        session["token"] = None


@io.on("celery:task:results")
@typecast
def celery_task_results(msg: CeleryTaskData):
    emit(msg.event, msg.data, to=msg.target)


@io.on("disconnect")
def disconnect():
    token = session["token"]
    ROOM_HOSTS = cache.get("ROOM_HOSTS")
    if token not in ROOM_HOSTS:
        return
    if request.sid in ROOM_HOSTS[token]:
        ROOM_HOSTS[token].remove(request.sid)
    if len(ROOM_HOSTS[token]) == 0:
        del ROOM_HOSTS[token]
    cache.set("ROOM_HOSTS", ROOM_HOSTS)
    # with contextlib.suppress(KeyError):

    #     PYCLIENTS = cache.get("pyclients")
    #     PYCLIENTS.pop(_get_uuid_for_sid(request.sid), None)

    #     cache.set("pyclients", PYCLIENTS)
    #     log.info(f"disconnect {request.sid} and updated PYCLIENTS to {PYCLIENTS}")

    #     # remove the pyclient from the dict
    #     names = cache.get(f"PER-TOKEN-NAME:{session['token']}") or {}
    #     connected_users = [{"name": names[sid]} for sid in ROOM_HOSTS[token]]
    #     emit(
    #         "connectedUsers",
    #         list(reversed(connected_users)),
    #         to=_webclients_room({"token": token}),
    #     )
    # log.debug(f"disconnect {request.sid} and updated HOSTS to {ROOM_HOSTS}")


@io.on("join")
@typecast
def join(data: JoinData):
    """
    Arguments:
        data: {"token": str, "uuid": str}
    """
    # only used by pyclients that only connect via socket (no HTML)
    # session["authenticated"] = data.auth_token == app.config["AUTH_TOKEN"]
    # log.debug(f"Client {request.sid} is {session['authenticated'] = }")
    # PYCLIENTS = cache.get("pyclients")
    # if data.uuid in PYCLIENTS:
    #     log.critical(
    #         f"UUID {data.uuid} is already registered in {app.config['pyclients']}."
    #     )
    #     emit("message:log", f"UUID {data.uuid} is already registered", to=request.sid)
    # PYCLIENTS[data.uuid] = request.sid
    # cache.set("pyclients", PYCLIENTS)
    session["token"] = data.token
    # join_room(f"pyclients_{data.token}")
    # PER_TOKEN_DATA = cache.get("PER-TOKEN-DATA")
    # if data.token not in PER_TOKEN_DATA:
    #     PER_TOKEN_DATA[data.token] = {}
    # cache.set("PER-TOKEN-DATA", PER_TOKEN_DATA)


@io.on("analysis:figure")
@typecast
def analysis_figure(data: AnalysisFigureData):
    emit(
        "analysis:figure",
        data.figure,
        include_self=False,
        to=f"webclients_{session['token']}",
    )


@io.on("scene:set")
@typecast
def scene_set(data: SceneSetData):
    emit(
        "scene:set", data.index, include_self=False, to=f"webclients_{session['token']}"
    )


@io.on("scene:step")
@typecast
def scene_step():
    return call("scene:step", to=f"webclients_{session['token']}")


@io.on("atoms:download")
@typecast
def atoms_download(data: AtomsDownloadData):
    return call("atoms:download", data.indices, to=get_main_room_host(session["token"]))


@io.on("atoms:upload")
def atoms_upload(data: dict):
    emit("atoms:upload", data, include_self=False, to=f"webclients_{session['token']}")


@io.on("atoms:delete")
@typecast
def atoms_delete(data: DeleteAtomsData):
    emit(
        "atoms:delete",
        data.index,
        include_self=False,
        to=f"webclients_{session['token']}",
    )


@io.on("atoms:length")
def atoms_length():
    return call("atoms:length", to=get_main_room_host(session["token"]))


@io.on("modifier:schema")
@typecast
def modifier_schema(data: SchemaData):
    emit("modifier:schema", data.schema, include_self=False, to=_webclients_room(data))


@io.on("selection:schema")
@typecast
def selection_schema(data: SchemaData):
    emit("selection:schema", data.schema, include_self=False, to=data.sid)


@io.on("scene:schema")
@typecast
def scene_schema(data: SchemaData):
    emit("scene:schema", data.schema, include_self=False, to=data.sid)


@io.on("draw:schema")
@typecast
def draw_schema(data: SchemaData):
    emit("draw:schema", data.schema, include_self=False, to=data.sid)


@io.on("points:get")
def scene_points():
    return call("points:get", to=get_main_room_host(session["token"]))


@io.on("scene:segments")
def scene_segments():
    return call("scene:segments", to=get_main_room_host(session["token"]))


@io.on("selection:get")
def selection_get():
    return call("selection:get", to=get_main_room_host(session["token"]))


@io.on("selection:set")
def selection_set(data: list[int]):
    emit(
        "selection:set",
        data,
        include_self=False,
        to=f"webclients_{session['token']}",
    )


@io.on("message:log")
@typecast
def message_log(data: MessageData):
    emit("message:log", data.message, to=f"webclients_{session['token']}")


@io.on("download:response")
def download_response(data):
    emit("download:response", data["data"], to=data["sid"])


@io.on("scene:play")
@typecast
def scene_play():
    emit("scene:play", to=f"webclients_{session['token']}")


@io.on("scene:pause")
@typecast
def scene_pause():
    emit("scene:pause", to=f"webclients_{session['token']}")


@io.on("bookmarks:get")
@typecast
def bookmarks_get():
    return call("bookmarks:get", to=get_main_room_host(session["token"]))


@io.on("bookmarks:set")
@typecast
def bookmarks_set(data: dict):
    emit(
        "bookmarks:set",
        data,
        include_self=False,
        to=f"webclients_{session['token']}",
    )


@io.on("points:set")
def points_set(data: list[list[float]]):
    emit("points:set", data, include_self=False, to=f"webclients_{session['token']}")


@io.on("debug")
def debug(data: dict):
    emit("debug", data, include_self=False, to=_webclients_room(data))


# @io.on("modifier:run:running")
# @typecast
# def modifier_run_running(data: ModifierRunRunningData):
#     MODIFIER = cache.get("MODIFIER")
#     MODIFIER["active"] = data.name
#     emit(
#         "modifier:run:running",
#         dataclasses.asdict(data),
#         include_self=False,
#         to=_webclients_room(data),
#     )


# @io.on("modifier:run:finished")
# @typecast
# def modifier_run_finished(data: AtomsLengthData):
#     # remove 0th element from queue
#     MODIFIER = cache.get("MODIFIER")
#     MODIFIER["queue"].pop(0)
#     MODIFIER["active"] = None
#     cache.set("MODIFIER", MODIFIER)
#     modifier_lock.release()
#     print("modifier_lock released")
#     emit(
#         "modifier:run:finished",
#         dataclasses.asdict(data),
#         include_self=False,
#         to=_webclients_room(data),
#     )


# @io.on("modifier:run:failed")
# def modifier_run_failed():
#     """Take care if the modifier does not respond."""
#     # remove 0th element from queue
#     MODIFIER = cache.get("MODIFIER")
#     MODIFIER["queue"].pop(0)
#     MODIFIER["active"] = None
#     cache.set("MODIFIER", MODIFIER)
#     modifier_lock.release()
#     log.critical("Modifier failed - releasing lock.")


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


@io.on("selection:run")
def selection_run(data: dict):
    """Run the selection."""
    tasks.run_selection.delay(request.url_root, session["token"], data)


@io.on("scene:trash")
def scene_trash():
    tasks.scene_trash.delay(request.url_root, session["token"])


@io.on("analysis:run")
def analysis_run(data: dict):
    """Run the analysis."""
    io.emit("analysis:run:enqueue", to=f"webclients_{session['token']}", include_self=False)
    tasks.run_analysis.delay(request.url_root, session["token"], data)


@io.on("modifier:run")
def modifier_run(data: dict):
    """Run the modifier."""
    io.emit("modifier:run:enqueue", to=f"webclients_{session['token']}", include_self=False)
    tasks.run_modifier.delay(request.url_root, session["token"], data)
