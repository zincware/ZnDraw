import dataclasses
import datetime
import logging
from threading import Lock
from typing import List
from uuid import uuid4

from celery import chain
from flask import current_app, request, session
from flask_socketio import call, emit, join_room
from socketio.exceptions import TimeoutError

from zndraw.db import Session
from zndraw.db import schema as db_schema
from zndraw.server import tasks
from zndraw.utils import typecast

from ..app import socketio as io
from ..data import (
    AnalysisFigureData,
    CeleryTaskData,
    DeleteAtomsData,
    FrameData,
    JoinData,
    MessageData,
    ModifierRegisterData,
    RoomGetData,
    SceneSetData,
    SceneUpdateData,
    SchemaData,
    SubscribedUserData,
)
from .utils import get_queue_position

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

        # check if there is a db_schema.Room with the given token, if not create one
        with Session() as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            if room is None:
                room = db_schema.Room(
                    token=token, currentStep=0, points=[], selection=[]
                )
                ses.add(room)

            client = db_schema.Client(
                sid=request.sid, room=room, name=uuid4().hex[:8].upper()
            )
            # check if any client in the room is host
            if (
                ses.query(db_schema.Client).filter_by(room=room, host=True).first()
                is None
            ):
                client.host = True
            ses.add(client)
            ses.commit()

        # ROOM_HOSTS = cache.get("ROOM_HOSTS")

        # if token not in ROOM_HOSTS:
        #     ROOM_HOSTS[token] = [request.sid]
        # else:
        #     ROOM_HOSTS[token].append(request.sid)

        # cache.set("ROOM_HOSTS", ROOM_HOSTS)

        # # data = {"sid": request.sid, "token": token}
        # # data["host"] = ROOM_HOSTS[token][0] == request.sid
        # names = cache.get(f"PER-TOKEN-NAME:{session['token']}") or {}
        # names[request.sid] = uuid4().hex[:8].upper()
        # cache.set(f"PER-TOKEN-NAME:{session['token']}", names)

        # connected_users = [{"name": names[sid]} for sid in ROOM_HOSTS[token]]

        # get all clients in the room
        with Session() as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            clients = ses.query(db_schema.Client).filter_by(room=room).all()
            connected_users = [{"name": client.name} for client in clients]

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

        # log.debug(f"connected {request.sid} and updated HOSTS to {ROOM_HOSTS}")
        # emit("message:log", "Connection established", to=request.sid)
        # PER_TOKEN_DATA = cache.get("PER-TOKEN-DATA")
        # if token not in PER_TOKEN_DATA:
        #     PER_TOKEN_DATA[token] = {}
        # cache.set("PER-TOKEN-DATA", PER_TOKEN_DATA)

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
        if msg.data[0]["update_database"]:
            tasks.update_atoms(token, msg.data)
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
    token = str(session["token"])

    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        client = ses.query(db_schema.Client).filter_by(sid=request.sid).first()
        if client is not None:
            ses.delete(client)
            if client.host:
                new_host = ses.query(db_schema.Client).filter_by(room=room).first()
                if new_host is not None:
                    new_host.host = True
            ses.commit()

    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        clients = ses.query(db_schema.Client).filter_by(room=room).all()
        connected_users = [{"name": client.name} for client in clients]

    emit(
        "connectedUsers",
        list(reversed(connected_users)),
        to=_webclients_room({"token": token}),
    )

    # ROOM_HOSTS = cache.get("ROOM_HOSTS")
    # if token in ROOM_HOSTS:
    #     if request.sid in ROOM_HOSTS[token]:
    #         ROOM_HOSTS[token].remove(request.sid)
    #     if len(ROOM_HOSTS[token]) == 0:
    #         del ROOM_HOSTS[token]
    #     cache.set("ROOM_HOSTS", ROOM_HOSTS)

    # pyclients only
    # check if the SID is in MODIFIER_HOSTS or ROOM_MODIFIER_HOSTS
    with Session() as ses:
        global_modifier_client = (
            ses.query(db_schema.GlobalModifierClient).filter_by(sid=request.sid).all()
        )
        for gmc in global_modifier_client:
            ses.delete(gmc)

        room_modifier_client = (
            ses.query(db_schema.RoomModifierClient).filter_by(sid=request.sid).all()
        )
        for rmc in room_modifier_client:
            ses.delete(rmc)
        ses.commit()

    # MODIFIER_HOSTS = cache.get("MODIFIER_HOSTS")
    # for name, sids in MODIFIER_HOSTS.items():
    #     while request.sid in sids:
    #         MODIFIER_HOSTS[name].remove(request.sid)
    # cache.set("MODIFIER_HOSTS", MODIFIER_HOSTS)

    # ROOM_MODIFIER_HOSTS = cache.get("ROOM_MODIFIER_HOSTS")
    # if token in ROOM_MODIFIER_HOSTS:
    #     for name, sids in ROOM_MODIFIER_HOSTS[token].items():
    #         while request.sid in sids:
    #             ROOM_MODIFIER_HOSTS[token][name].remove(request.sid)

    # # !!!!!!!!!!!!!!
    # # TODO: remove schema if no more hosts are connected and send to webclients

    # cache.set("ROOM_MODIFIER_HOSTS", ROOM_MODIFIER_HOSTS)


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
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        return room.currentStep


@io.on("atoms:download")
def atoms_download(indices: list[int]):
    token = str(session["token"])

    data = {}
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        # get all db_schema.Frame with the given indices
        frames = ses.query(db_schema.Frame).filter(
            db_schema.Frame.index.in_(indices), db_schema.Frame.room == room
        )
        for frame in frames:
            data[frame.index] = frame.data

    return data


@io.on("atoms:upload")
def atoms_upload(data: List[FrameData]):
    token = str(session["token"])
    if data[0]["update_database"]:
        tasks.update_atoms(token, data)
    emit(
        "atoms:upload",
        data,
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
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        ses.query(db_schema.Frame).filter(
            db_schema.Frame.index.in_(data.index), db_schema.Frame.room == room
        ).delete(synchronize_session=False)
        ses.commit()


@io.on("atoms:length")
def atoms_length():
    token = str(session["token"])
    with Session() as ses:
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
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        return room.points


@io.on("scene:segments")
def scene_segments():
    token = str(session["token"])
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        host = ses.query(db_schema.Client).filter_by(room=room, host=True).first()
    return call("scene:segments", to=host.sid)


@io.on("selection:get")
def selection_get():
    token = str(session["token"])
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
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
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
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
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
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
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
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
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
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
    token = str(session["token"])
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        current_client = ses.query(db_schema.Client).filter_by(sid=request.sid).first()
        controller_client = (
            ses.query(db_schema.Client).filter_by(name=data.user).first()
        )
        if (
            controller_client.sid == request.sid
            or controller_client.step_controller == current_client
        ):
            current_client.step_controller = None
        else:
            current_client.step_controller = controller_client
        ses.commit()


@io.on("connectedUsers:subscribe:camera")
@typecast
def connected_users_subscribe_camera(data: SubscribedUserData):
    """
    Subscribe to camera updates for connected users.

    data: {user: str}
    """
    token = str(session["token"])
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        current_client = ses.query(db_schema.Client).filter_by(sid=request.sid).first()
        controller_client = (
            ses.query(db_schema.Client).filter_by(name=data.user).first()
        )
        if (
            controller_client.sid == request.sid
            or controller_client.camera_controller == current_client
        ):
            current_client.camera_controller = None
        else:
            current_client.camera_controller = controller_client
        ses.commit()


@io.on("scene:update")
@typecast
def scene_update(data: SceneUpdateData):
    """Update the scene.

    data: {step: int, camera: {position: [float, float, float], rotation: [float, float, float]}}
    """
    token = str(session["token"])

    if data.step is not None:
        with Session() as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            current_client = (
                ses.query(db_schema.Client).filter_by(sid=request.sid).first()
            )
            if room is None:
                raise ValueError("No room found for token.")
            if current_client.host:
                room.currentStep = data.step
                ses.commit()

            step_subscribers = (
                ses.query(db_schema.Client)
                .filter_by(step_controller=current_client)
                .all()
            )
            step_subscribers = [client.sid for client in step_subscribers]

        emit(
            "scene:update",
            dataclasses.asdict(data),
            include_self=False,
            to=step_subscribers,
        )

    if data.camera is not None:
        with Session() as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            if room is None:
                raise ValueError("No room found for token.")
            current_client = (
                ses.query(db_schema.Client).filter_by(sid=request.sid).first()
            )
            camera_subscribers = (
                ses.query(db_schema.Client)
                .filter_by(camera_controller=current_client)
                .all()
            )
            camera_subscribers = [client.sid for client in camera_subscribers]

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
    # split into separate streams based on the modifier name
    tasks.run_modifier(
        f"http://127.0.0.1:{current_app.config['PORT']}", session["token"], data
    )


def _register_global_modifier(data):
    with Session() as ses:
        global_modifier = (
            ses.query(db_schema.GlobalModifier).filter_by(name=data.name).first()
        )
        if global_modifier is None:
            global_modifier = db_schema.GlobalModifier(
                name=data.name, schema=data.schema
            )
            ses.add(global_modifier)
        elif global_modifier.schema != data.schema:
            log.critical(
                f"Modifier {data.name} is already registered with a different schema."
            )
            return
        # attach GlobalModifierClient
        global_modifier_client = db_schema.GlobalModifierClient(
            sid=request.sid,
            timeout=data.timeout,
            available=False,
            global_modifier=global_modifier,
        )
        ses.add(global_modifier_client)
        ses.commit()


def _register_room_modifier(data):
    # TODO: do we want one table for global and room modifiers?
    token = str(session["token"])
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        room_modifier = (
            ses.query(db_schema.RoomModifier).filter_by(name=data.name).first()
        )
        if room_modifier is None:
            room_modifier = db_schema.RoomModifier(
                name=data.name, schema=data.schema, room=room
            )
            ses.add(room_modifier)
        elif room_modifier.schema != data.schema:
            log.critical(
                f"Modifier {data.name} is already registered with a different schema."
            )
            return
        # attach RoomModifierClient
        room_modifier_client = db_schema.RoomModifierClient(
            sid=request.sid,
            timeout=data.timeout,
            available=False,
            room_modifier=room_modifier,
        )
        ses.add(room_modifier_client)
        ses.commit()


@io.on("modifier:register")
@typecast
def modifier_register(data: ModifierRegisterData):
    """Register the modifier."""
    if data.default:
        # TODO: authenticattion
        _register_global_modifier(data)
    else:
        _register_room_modifier(data)

    # emit new modifier schema to all webclients
    tasks.modifier_schema.delay(
        f"http://127.0.0.1:{current_app.config['PORT']}", session["token"]
    )


@io.on("modifier:available")
def modifier_available(available: bool):
    """Update the modifier availability."""
    with Session() as ses:
        global_modifier_client = (
            ses.query(db_schema.GlobalModifierClient).filter_by(sid=request.sid).first()
        )
        room_modifier_client = (
            ses.query(db_schema.RoomModifierClient).filter_by(sid=request.sid).first()
        )
        if global_modifier_client is not None:
            global_modifier_client.available = available
        elif room_modifier_client is not None:
            room_modifier_client.available = available
        ses.commit()


@io.on("modifier:queue:update")
def modifier_queue_update(data: dict):
    """Update the modifier queue."""
    queue_position = get_queue_position(data["queue_name"], data["job_id"])
    emit("modifier:queue:update", queue_position, to=f"webclients_{session['token']}")


@io.on("ping")
def ping() -> str:
    return "pong"


@io.on("room:get")
def room_get(data: RoomGetData):
    tasks.handle_room_get.delay(data, session["token"], request.url_root, request.sid)
