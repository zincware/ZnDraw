import dataclasses
import importlib
import importlib.metadata
import importlib.util
import logging
import uuid

import znsocket
from flask import current_app, request, session
from flask_socketio import SocketIO, emit, join_room
from redis import Redis

from zndraw.tasks import (
    inspect_zntrack_node,
    load_zntrack_figures,
    load_zntrack_frames,
    run_room_copy,
    run_room_worker,
    run_scene_dependent_schema,
    run_schema,
    run_upload_file,
)

log = logging.getLogger(__name__)
__version__ = importlib.metadata.version("zndraw")


@dataclasses.dataclass
class DummyClient:
    """Dummy replacement for znsocket.Client."""

    sio: SocketIO
    refresh_callbacks: list = dataclasses.field(default_factory=list)


def init_socketio_events(io: SocketIO):
    @io.on("connect")
    def connect():
        emit("version", __version__)

    @io.on("shutdown")
    def shutdown():
        if "AUTH_TOKEN" not in current_app.config or session["authenticated"]:
            log.critical("Shutting down server")
            current_app.extensions["redis"].flushall()

            current_app.extensions["celery"].control.purge()
            current_app.extensions["celery"].control.broadcast("shutdown")

            io.stop()
        else:
            log.critical("Unauthenticated user tried to shut down the server.")

    @io.on("disconnect")
    def disconnect():
        try:
            room = str(session["token"])
        except KeyError:
            log.critical(f"disconnecting {request.sid}")
            return
        r = current_app.extensions["redis"]

        if "name" in session:
            log.critical(f"disconnecting (webclient) {request.sid} from room {room}")
            r.srem(f"room:{room}:webclients", session["name"])
            emit(
                "room:users:refresh", list(r.smembers(f"room:{room}:webclients")), to=room
            )
        else:
            log.critical(f"disconnecting (pyclient) {request.sid}")

        for ref_token in [room, "default"]:
            modifier_registry = znsocket.Dict(
                r=r, key=f"registry:{ref_token}:modifier", repr_type="full"
            )

            modifier_schema = znsocket.Dict(
                r=r,
                key=f"schema:{ref_token}:modifier",
                repr_type="full",
                socket=DummyClient(sio=io),
            )

            # TODO: default room modifier
            for modifier in modifier_registry.pop(request.sid, []):
                for other in modifier_registry.values():
                    if modifier in other:
                        break
                else:
                    log.debug(f"Remove {modifier} from room {room}")
                    modifier_schema.pop(
                        modifier
                    )  # TODO this does not work well with public yet.

    @io.on("webclient:connect")
    def webclient_connect():
        try:
            token = session["token"]
        except KeyError:
            token = uuid.uuid4().hex[:8]
            session["token"] = token

        if "AUTH_TOKEN" not in current_app.config:
            session["authenticated"] = True
        else:
            if "authenticated" not in session:
                session["authenticated"] = False

        room = str(session["token"])
        join_room(room)  # rename token to room or room_token

        run_schema.delay(room)
        run_scene_dependent_schema.delay(room)

        session["name"] = uuid.uuid4().hex[:8]

        r = current_app.extensions["redis"]
        r.sadd(f"room:{room}:webclients", session["name"])

        # TODO: this is currently not used afaik
        emit("room:users:refresh", list(r.smembers(f"room:{room}:webclients")), to=room)

        log.critical(f"connecting (webclient) {request.sid} to {room}")

        if "TUTORIAL" in current_app.config:
            emit("tutorial:url", current_app.config["TUTORIAL"])
        if "SIMGEN" in current_app.config:
            emit("showSiMGen", True)

        return {
            "name": session["name"],
            "room": room,
            "authenticated": session["authenticated"],
        }

    @io.on("join")  # rename pyclient:connect
    def join(data: dict):
        """
        Arguments:
            data: {"token": str, "auth_token": str}
        """
        # TODO: prohibit "token" to be "default"

        if "AUTH_TOKEN" not in current_app.config:
            session["authenticated"] = True
        else:
            session["authenticated"] = (
                data["auth_token"] == current_app.config["AUTH_TOKEN"]
            )
        token = data["token"]
        session["token"] = token
        room = str(session["token"])

        join_room(room)
        log.critical(f"connecting (pyclient) {request.sid} to {room}")
        # join_room(f"pyclients_{token}")

    @io.on("room:worker:run")
    def room_worker_run():
        """Start a worker to process all (available) queued tasks."""
        room = session.get("token")
        run_room_worker.delay(room)

    @io.on("schema:refresh")
    def schema_refresh():
        room = session.get("token")
        run_scene_dependent_schema.delay(room)

    @io.on("room:alert")
    def room_alert(msg: str):
        """Forward the alert message to every client in the room"""
        # TODO: identify the source client.
        room = session.get("token")
        emit("room:alert", msg, to=room)

    @io.on("room:upload:file")
    def room_upload_file(data: dict):
        room = session.get("token")
        run_upload_file.delay(room, data)

    @io.on("room:lock:set")
    def room_lock_set(locked: bool):
        room = session.get("token")
        r: Redis = current_app.extensions["redis"]
        r.set(f"room:{room}:locked", str(locked))
        emit("room:lock:set", locked, to=room)

    @io.on("room:lock:get")
    def room_lock_get() -> bool:
        room = session.get("token")
        r: Redis = current_app.extensions["redis"]
        locked = r.get(f"room:{room}:locked")
        return locked == "True"

    @io.on("room:token:get")
    def room_token_get() -> str:
        return session.get("token")

    @io.on("zntrack:available")
    def check_zntrack_available() -> bool:
        return importlib.util.find_spec("zntrack") is not None

    @io.on("zntrack:list-stages")
    def zntrack_list_stages(data: dict):
        try:
            import dvc.api

            fs = dvc.api.DVCFileSystem(url=data.get("remote"), rev=data.get("rev"))
            return [x.name for x in fs.repo.stage.collect() if hasattr(x, "name")]
        except Exception:
            return []

    @io.on("zntrack:inspect-stage")
    def zntrack_inspect_stage(data: dict):
        return inspect_zntrack_node(**data)

    @io.on("zntrack:load-frames")
    def zntrack_load_frames(data: dict):
        load_zntrack_frames.delay(room=session.get("token"), **data)

    @io.on("zntrack:load-figure")
    def zntrack_load_figure(data: dict):
        load_zntrack_figures.delay(room=session.get("token"), **data)

    @io.on("room:copy")
    def room_copy():
        room = session.get("token")
        run_room_copy.delay(room)
