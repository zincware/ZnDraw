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
    MessageData,
    ModifierRegisterData,
    RoomGetData,
    RoomSetData,
    SceneSetData,
    SceneUpdateData,
    SchemaData,
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


@io.on("connect")
def connect():
    try:
        token = str(session["token"])
        log.debug(f"connecting (webclient) {request.sid}")
        URL = f"http://127.0.0.1:{current_app.config['PORT']}"
        # if you connect through Python, you don't have a token
        read_file_chain = chain(
            tasks.read_file.s(
                url=URL,
                target=request.sid,
                token=token,
                fileio=current_app.config["FileIO"],
            ),
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

            client = db_schema.WebClient(
                sid=request.sid, room=room, name=uuid4().hex[:8].upper()
            )
            # check if any client in the room is host
            if (
                ses.query(db_schema.WebClient).filter_by(room=room, host=True).first()
                is None
            ):
                client.host = True
            ses.add(client)
            ses.commit()

        # get all clients in the room
        with Session() as ses:
            room = ses.query(db_schema.Room).filter_by(token=token).first()
            clients = ses.query(db_schema.WebClient).filter_by(room=room).all()
            connected_users = [{"name": client.name} for client in clients]

        emit(
            "connectedUsers",
            list(reversed(connected_users)),
            to=_webclients_room({"token": token}),
        )

        # append to zndraw.log a line isoformat() + " " + token
        log.info(
            datetime.datetime.now().isoformat() + " " + token
            if token
            else "client" + " connected"
        )

    except KeyError:
        # clients that connect directly via socketio do not call "join" to
        # register their token
        log.debug(f"connecting (pyclient) {request.sid}")
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
    token = str(session.get("token"))

    url = request.url_root
    if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
        url = url.replace("http://", "https://")
    url = url.replace("http", "ws")

    tasks.on_disconnect.delay(token=token, sid=request.sid, url=url)


@io.on("join")
def join(data: dict):
    """
    Arguments:
        data: {"token": str, "auth_token": str}
    """
    if current_app.config["AUTH_TOKEN"] is None:
        session["authenticated"] = True
    else:
        session["authenticated"] = (
            data["auth_token"] == current_app.config["AUTH_TOKEN"]
        )
    token = data["token"]
    session["token"] = token
    join_room(f"{token}")
    join_room(f"pyclients_{token}")

    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            room = db_schema.Room(token=token, currentStep=0, points=[], selection=[])
            ses.add(room)
        ses.commit()


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
        host = ses.query(db_schema.WebClient).filter_by(room=room, host=True).first()
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
        current_client = (
            ses.query(db_schema.WebClient).filter_by(sid=request.sid).first()
        )
        controller_client = (
            ses.query(db_schema.WebClient).filter_by(name=data.user).first()
        )
        if (
            controller_client.sid == request.sid
            or controller_client.camera_controller == current_client
        ):
            current_client.camera_controller = None
        else:
            current_client.camera_controller = controller_client
        ses.commit()


@io.on("camera:update")
def camera_update(data: dict):
    token = str(session["token"])
    timestamp = datetime.datetime.utcnow().isoformat()
    session["camera-update"] = timestamp
    # TODO: store this in the session
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        if room is None:
            raise ValueError("No room found for token.")
        current_client = (
            ses.query(db_schema.WebClient).filter_by(sid=request.sid).first()
        )
        camera_subscribers = (
            ses.query(db_schema.WebClient)
            .filter_by(camera_controller=current_client)
            .all()
        )
        camera_subscribers = [client.sid for client in camera_subscribers]

    emit(
        "camera:update",
        data,
        include_self=False,
        to=camera_subscribers,
    )
    io.sleep(1)
    if session["camera-update"] == timestamp:
        data = RoomSetData(camera=data)
        url = request.url_root
        if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
            url = url.replace("http://", "https://")
        tasks.handle_room_set.delay(data.to_dict(), session["token"], url, request.sid)


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
                ses.query(db_schema.WebClient).filter_by(sid=request.sid).first()
            )
            if room is None:
                raise ValueError("No room found for token.")
            if current_client.host:
                room.currentStep = data.step
                ses.commit()

            step_subscribers = (
                ses.query(db_schema.WebClient)
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
                ses.query(db_schema.WebClient).filter_by(sid=request.sid).first()
            )
            camera_subscribers = (
                ses.query(db_schema.WebClient)
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
    url = f"http://127.0.0.1:{current_app.config['PORT']}"
    tasks.run_modifier(url, session["token"], data)


def _register_global_modifier(data):
    with Session() as ses:
        global_modifier = (
            ses.query(db_schema.Modifier)
            .filter_by(name=data.name, room_token=None)
            .first()
        )
        if global_modifier is None:
            global_modifier = db_schema.Modifier(
                name=data.name, schema=data.schema, room_token=None
            )
            ses.add(global_modifier)
        elif global_modifier.schema != data.schema:
            log.critical(
                f"Modifier {data.name} is already registered with a different schema."
            )
            return
        # attach GlobalModifierClient
        global_modifier_client = db_schema.ModifierClient(
            sid=request.sid,
            timeout=data.timeout,
            available=False,
            modifier=global_modifier,
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
            ses.query(db_schema.Modifier).filter_by(name=data.name, room=room).first()
        )
        if room_modifier is None:
            room_modifier = db_schema.Modifier(
                name=data.name, schema=data.schema, room=room
            )
            ses.add(room_modifier)
        elif room_modifier.schema != data.schema:
            log.critical(
                f"Modifier {data.name} is already registered with a different schema."
            )
            return
        # attach RoomModifierClient
        room_modifier_client = db_schema.ModifierClient(
            sid=request.sid,
            timeout=data.timeout,
            available=False,
            modifier=room_modifier,
        )
        ses.add(room_modifier_client)
        ses.commit()


@io.on("modifier:register")
@typecast
def modifier_register(data: ModifierRegisterData):
    """Register the modifier."""
    if data.default and not session["authenticated"]:
        log.critical("Unauthenticated user tried to register a default modifier.")
        return
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
    tasks.activate_modifier.delay(request.sid, available)


@io.on("ping")
def ping() -> str:
    return "pong"


@io.on("room:get")
def room_get(data: RoomGetData):
    url = request.url_root
    if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
        url = url.replace("http://", "https://")
    tasks.handle_room_get.delay(data, session["token"], url, request.sid)


@io.on("room:set")
@typecast
def room_set(data: RoomSetData):
    emit(
        "room:set",
        data.to_dict(),
        include_self=False,
        to=f"webclients_{session['token']}",
    )
    url = request.url_root
    if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
        url = url.replace("http://", "https://")

    if data.update_database:
        # TODO: we need to differentiate, if the data comes from a pyclient or a webclient
        # TODO: for fast updates, e.g. points, step during play this is not fast enough
        tasks.handle_room_set.delay(data.to_dict(), session["token"], url, request.sid)


@io.on("step:update")
def step_update(step: int):
    timestamp = datetime.datetime.utcnow().isoformat()
    session["step-update"] = timestamp

    data = RoomSetData(step=step)

    emit(
        "room:set",
        data.to_dict(),
        include_self=False,
        to=f"webclients_{session['token']}",
    )
    io.sleep(1)
    if session["step-update"] == timestamp:
        url = request.url_root
        if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
            url = url.replace("http://", "https://")
        tasks.handle_room_set.delay(data.to_dict(), session["token"], url, request.sid)


@io.on("points:update")
def points_update(points: list[list[float]]):
    timestamp = datetime.datetime.utcnow().isoformat()
    session["points-update"] = timestamp

    data = RoomSetData(points=points)

    emit(
        "room:set",
        data.to_dict(),
        include_self=False,
        to=f"webclients_{session['token']}",
    )
    io.sleep(1)
    if session["points-update"] == timestamp:
        url = request.url_root
        if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
            url = url.replace("http://", "https://")
        tasks.handle_room_set.delay(data.to_dict(), session["token"], url, request.sid)


@io.on("file:upload")
def file_upload(data: dict):
    url = request.url_root
    if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
        url = url.replace("http://", "https://")
    tasks.upload_file.delay(
        url=url,
        token=str(session["token"]),
        filename=data["filename"],
        content=data["content"],
    )


@io.on("file:download")
def file_download():
    url = request.url_root
    if current_app.config["upgrade_insecure_requests"] and not "127.0.0.1" in url:
        url = url.replace("http://", "https://")
    tasks.download_file.delay(url=url, token=str(session["token"]), sid=request.sid)
