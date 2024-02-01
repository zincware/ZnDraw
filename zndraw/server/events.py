import dataclasses
import datetime
import logging
from threading import Lock
from uuid import uuid4

from celery import chain
from flask import current_app, request, session
from flask_socketio import call, emit, join_room
from socketio.exceptions import TimeoutError
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from zndraw.db import schema as db_schema
from zndraw.server import tasks
from zndraw.settings import GlobalConfig
from zndraw.utils import typecast

from ..app import cache
from ..app import socketio as io
from ..data import (
    AnalysisFigureData,
    CeleryTaskData,
    DeleteAtomsData,
    FrameData,
    JoinData,
    MessageData,
    ModifierRegisterData,
    SceneSetData,
    SceneUpdateData,
    SchemaData,
    SubscribedUserData,
)

log = logging.getLogger(__name__)

modifier_lock = Lock()

DB_PATH = GlobalConfig.load().database.get_path()
engine = create_engine(f"sqlite:///{DB_PATH}")


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
        token = str(session["token"])
        URL = f"http://127.0.0.1:{current_app.config['PORT']}"
        # if you connect through Python, you don't have a token
        read_file_chain = chain(
            tasks.read_file.s(URL, request.sid, token),
            tasks.analysis_schema.si(URL, token),
        )
        read_file_chain.delay()
        tasks.get_selection_schema.delay(URL, request.sid)
        tasks.scene_schema.delay(URL, request.sid)
        tasks.geometries_schema.delay(URL, request.sid)
        tasks.modifier_schema.delay(URL, token)

        join_room(f"webclients_{token}")
        join_room(f"{token}")
        # who ever connected latest is the HOST of the room

        # with Session(engine) as session:
        #     session.add(some_object)
        #     session.add(some_other_object)
        #     session.commit()

        # check if there is a db_schema.Room with the given toke, if not create one
        with Session(engine) as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            if room is None:
                room = db_schema.Room(
                    token=token, currentStep=0, points=[], selection=[]
                )
                ses.add(room)
                ses.commit()

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


@io.on("celery:task:emit")
@typecast
def celery_task_results(msg: CeleryTaskData):
    token = str(session["token"])
    if msg.event == "atoms:upload":
        with Session(engine) as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            # check if the index is already in the db
            frame = (
                ses.query(db_schema.Frame)
                .filter_by(index=msg.data["index"], room=room)
                .first()
            )
            if frame is not None:
                # if so, update the data
                frame.data = msg.data["data"]
            else:
                # create a db_schema.Frame
                frame = db_schema.Frame(
                    index=msg.data["index"], data=msg.data["data"], room=room
                )
                ses.add(frame)
            ses.commit()
    emit(msg.event, msg.data, to=msg.target)
    if msg.disconnect:
        io.server.disconnect(request.sid)


@io.on("celery:task:call")
@typecast
def celery_task_call(msg: CeleryTaskData):
    try:
        return call(msg.event, msg.data, to=msg.target, timeout=msg.timeout)
    except TimeoutError:
        log.critical(f"TimeoutError for {msg.event} with {msg.data}")
        return None


@io.on("disconnect")
def disconnect():
    log.critical(f"-------------->>> disconnecting {request.sid}")
    token = session["token"]
    ROOM_HOSTS = cache.get("ROOM_HOSTS")
    if token in ROOM_HOSTS:
        if request.sid in ROOM_HOSTS[token]:
            ROOM_HOSTS[token].remove(request.sid)
        if len(ROOM_HOSTS[token]) == 0:
            del ROOM_HOSTS[token]
        cache.set("ROOM_HOSTS", ROOM_HOSTS)

    # pyclients only
    # check if the SID is in MODIFIER_HOSTS or ROOM_MODIFIER_HOSTS
    MODIFIER_HOSTS = cache.get("MODIFIER_HOSTS")
    for name, sids in MODIFIER_HOSTS.items():
        while request.sid in sids:
            MODIFIER_HOSTS[name].remove(request.sid)
    cache.set("MODIFIER_HOSTS", MODIFIER_HOSTS)

    ROOM_MODIFIER_HOSTS = cache.get("ROOM_MODIFIER_HOSTS")
    if token in ROOM_MODIFIER_HOSTS:
        for name, sids in ROOM_MODIFIER_HOSTS[token].items():
            while request.sid in sids:
                ROOM_MODIFIER_HOSTS[token][name].remove(request.sid)

    # !!!!!!!!!!!!!!
    # TODO: remove schema if no more hosts are connected and send to webclients

    cache.set("ROOM_MODIFIER_HOSTS", ROOM_MODIFIER_HOSTS)


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
    join_room(f"{data.token}")
    join_room(f"pyclients_{data.token}")
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
def scene_step():
    token = str(session["token"])
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return 0
        return room.currentStep


@io.on("atoms:download")
def atoms_download(indices: list[int]):
    token = str(session["token"])

    data = {}
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return
        # get all db_schema.Frame with the given indices
        frames = ses.query(db_schema.Frame).filter(
            db_schema.Frame.index.in_(indices), db_schema.Frame.room == room
        )
        for frame in frames:
            data[frame.index] = frame.data

    return data


@io.on("atoms:upload")
@typecast
def atoms_upload(data: FrameData):
    token = str(session["token"])
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        frame = (
            ses.query(db_schema.Frame).filter_by(index=data.index, room=room).first()
        )
        if frame is not None:
            # if so, update the data
            frame.data = data.data
        else:
            # create a db_schema.Frame
            frame = db_schema.Frame(index=data.index, data=data.data, room=room)
            ses.add(frame)
        ses.commit()
    emit(
        "atoms:upload",
        dataclasses.asdict(data),
        include_self=False,
        to=f"webclients_{session['token']}",
    )


@io.on("atoms:delete")
@typecast
def atoms_delete(data: DeleteAtomsData):
    token = str(session["token"])
    emit(
        "atoms:delete",
        data.index,
        include_self=False,
        to=f"webclients_{session['token']}",
    )
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        ses.query(db_schema.Frame).filter(
            db_schema.Frame.index.in_(data.index), db_schema.Frame.room == room
        ).delete(synchronize_session=False)
        ses.commit()



@io.on("atoms:length")
def atoms_length():
    token = str(session["token"])
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        return ses.query(db_schema.Frame).filter_by(room=room).count()


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
    token = str(session["token"])
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return []
        return room.points


@io.on("scene:segments")
def scene_segments():
    return call("scene:segments", to=get_main_room_host(session["token"]))


@io.on("selection:get")
def selection_get():
    token = str(session["token"])
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return []
        return room.selection


@io.on("selection:set")
def selection_set(data: list[int]):
    token = str(session["token"])
    emit(
        "selection:set",
        data,
        include_self=False,
        to=f"webclients_{token}",
    )
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return
        room.selection = data
        ses.commit()


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
def bookmarks_get():
    token = str(session["token"])
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return []
        return {bm.step: bm.text for bm in room.bookmarks}


@io.on("bookmarks:set")
@typecast
def bookmarks_set(data: dict):
    token = str(session["token"])
    emit(
        "bookmarks:set",
        data,
        include_self=False,
        to=f"webclients_{token}",
    )
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return
        # remove all bookmarks for the given token
        ses.query(db_schema.Bookmark).filter_by(room=room).delete()
        # add all bookmarks from the data
        for k, v in data.items():
            bookmark = db_schema.Bookmark(step=k, text=v, room=room)
            ses.add(bookmark)

        ses.commit()


@io.on("points:set")
def points_set(data: list[list[float]]):
    token = str(session["token"])
    emit("points:set", data, include_self=False, to=f"webclients_{token}")
    with Session(engine) as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            return
        room.points = data
        ses.commit()


@io.on("debug")
def debug(data: dict):
    emit("debug", data, include_self=False, to=_webclients_room(data))


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
    token = str(session["token"])

    if data.step is not None:
        step_subscribers = _get_subscribers(token, "STEP")
        emit(
            "scene:update",
            dataclasses.asdict(data),
            include_self=False,
            to=step_subscribers,
        )
        with Session(engine) as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            if room is None:
                return
            room.currentStep = data.step
            ses.commit()

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
    tasks.run_selection.delay(
        f"http://127.0.0.1:{current_app.config['PORT']}", session["token"], data
    )


@io.on("scene:trash")
def scene_trash():
    tasks.scene_trash.delay(
        f"http://127.0.0.1:{current_app.config['PORT']}", session["token"]
    )


@io.on("analysis:run")
def analysis_run(data: dict):
    """Run the analysis."""
    io.emit(
        "analysis:run:enqueue", to=f"webclients_{session['token']}", include_self=False
    )
    tasks.run_analysis.delay(
        f"http://127.0.0.1:{current_app.config['PORT']}", session["token"], data
    )


@io.on("modifier:run")
def modifier_run(data: dict):
    """Run the modifier."""
    io.emit(
        "modifier:run:enqueue", to=f"webclients_{session['token']}", include_self=False
    )
    tasks.run_modifier.apply_async(
        (f"http://127.0.0.1:{current_app.config['PORT']}", session["token"], data)
    )


@io.on("modifier:register")
@typecast
def modifier_register(data: ModifierRegisterData):
    """Register the modifier."""

    if data.default:
        MODIFIER_HOSTS = cache.get("MODIFIER_HOSTS")
        MODIFIER_SCHEMA = cache.get("MODIFIER_SCHEMA")
    else:
        MODIFIER_HOSTS = cache.get("ROOM_MODIFIER_HOSTS").get(session["token"], {})
        MODIFIER_SCHEMA = cache.get("ROOM_MODIFIER_SCHEMA").get(session["token"], {})

    # TODO: do not allow modifiers that are already in defaults, handle duplicates better!
    if data.name in MODIFIER_HOSTS:
        if MODIFIER_SCHEMA[data.name] != data.schema:
            log.critical(
                f"Modifier {data.name} is already registered with a different schema."
            )
            return
        log.info(f"Register additional handlfer for already registered {data.name}.")
        MODIFIER_HOSTS[data.name].append(request.sid)
    else:
        log.info(f"Register new handler for {data.name}.")
        MODIFIER_HOSTS[data.name] = [request.sid]
        MODIFIER_SCHEMA[data.name] = data.schema

    if data.default:
        cache.set("MODIFIER_HOSTS", MODIFIER_HOSTS)
        cache.set("MODIFIER_SCHEMA", MODIFIER_SCHEMA)
    else:
        ROOM_MODIFIER_HOSTS = cache.get("ROOM_MODIFIER_HOSTS")
        ROOM_MODIFIER_SCHEMA = cache.get("ROOM_MODIFIER_SCHEMA")

        ROOM_MODIFIER_HOSTS[session["token"]] = MODIFIER_HOSTS
        ROOM_MODIFIER_SCHEMA[session["token"]] = MODIFIER_SCHEMA

        cache.set("ROOM_MODIFIER_HOSTS", ROOM_MODIFIER_HOSTS)
        cache.set("ROOM_MODIFIER_SCHEMA", ROOM_MODIFIER_SCHEMA)

        # emit new modifier schema to all webclients
        tasks.modifier_schema.delay(
            f"http://127.0.0.1:{current_app.config['PORT']}", session["token"]
        )


@io.on("modifier:available")
def modifier_available(available: bool):
    """Update the modifier availability."""
    MODIFIER_AVAILABLE = cache.get("MODIFIER_AVAILABLE")
    MODIFIER_AVAILABLE[request.sid] = available
    cache.set("MODIFIER_AVAILABLE", MODIFIER_AVAILABLE)


@io.on("modifier:timeout")
def modifier_timeout(timeout: int):
    """Update the modifier timeout."""
    MODIFIER_TIMEOUT = cache.get("MODIFIER_TIMEOUT")
    MODIFIER_TIMEOUT[request.sid] = timeout
    cache.set("MODIFIER_TIMEOUT", MODIFIER_TIMEOUT)
