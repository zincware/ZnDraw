import importlib.metadata
import json
import logging
import uuid

import znsocket
from flask import current_app, request, session
from flask_socketio import SocketIO, emit, join_room
from redis import Redis

from zndraw.analyse import Analysis
from zndraw.draw import Geometry
from zndraw.modify import Modifier
from zndraw.scene import Scene
from zndraw.selection import Selection
from zndraw.tasks import (
    run_analysis,
    run_geometry,
    run_modifier,
    run_selection,
    run_upload_file,
)
from zndraw.utils import get_cls_from_json_schema, get_schema_with_instance_defaults

log = logging.getLogger(__name__)
__version__ = importlib.metadata.version("zndraw")


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

            # Get all modifier names associated with the client from Redis
            room_modifier_names = r.hgetall(f"room:{room}:modifiers")
            default_modifier_names = r.hgetall("room:default:modifiers")

            for modifier in room_modifier_names:
                r.srem(f"room:{room}:modifiers:{modifier}", request.sid)
                if r.scard(f"room:{room}:modifiers:{modifier}") == 0:
                    r.hdel(f"room:{room}:modifiers", modifier)
            for modifier in default_modifier_names:
                r.srem(f"room:default:modifiers:{modifier}", request.sid)
                if r.scard(f"room:default:modifiers:{modifier}") == 0:
                    r.hdel("room:default:modifiers", modifier)

            # Emit an event to refresh the modifier schema
            emit("modifier:schema:refresh", broadcast=True)

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

    @io.on("room:all:frames:refresh")
    def room_all_frames_refresh(indices: list[int]):
        emit("room:frames:refresh", [int(x) for x in indices], broadcast=True)

    @io.on("room:frames:insert")
    def room_frames_insert(data: dict):
        index = data.pop("index")
        value = data.pop("value")
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")

        lst = znsocket.List(r, f"room:{room}:frames")
        if not r.exists(f"room:{room}:frames"):
            default_lst = znsocket.List(r, "room:default:frames")
            # TODO: using a redis copy action would be faster
            lst.extend(default_lst)
        lst.insert(index, value)

        # not sure how to update, insert requires everything to be updated after the insertion
        # can be done custom on the client side to avoid resending everything
        emit("room:frames:refresh", [int(x) for x in data], to=room)
        # This method should be called, because it can move frames from the default
        # room to the current room. Doing so in the background
        # can cause issues with further operations on the frames. (see room frames set)
        return "OK"

    @io.on("room:length:get")
    def room_frames_length_get() -> int:
        room = session.get("token")
        r: Redis = current_app.extensions["redis"]
        room_key = (
            f"room:{room}:frames"
            if r.exists(f"room:{room}:frames")
            else "room:default:frames"
        )
        return len(znsocket.List(r, room_key))

    @io.on("modifier:register")
    def modifier_register(data: dict):
        """Register the modifier."""
        if data["public"] and not session["authenticated"]:
            log.critical("Unauthenticated user tried to register a default modifier.")
            return

        r: Redis = current_app.extensions["redis"]
        room = session.get("token")

        data["ZNDRAW_CLIENT_SID"] = request.sid

        if data.pop("public"):
            r.hset("room:default:modifiers", data["name"], json.dumps(data))
            r.sadd(f"room:default:modifiers:{data['name']}", request.sid)
        else:
            r.hset(f"room:{room}:modifiers", data["name"], json.dumps(data))
            r.sadd(f"room:{room}:modifiers:{data['name']}", request.sid)

        # only if public this is required
        emit("modifier:schema:refresh", broadcast=True)

    @io.on("modifier:schema")
    def modifier_schema():
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")

        modifiers: dict = r.hgetall(f"room:{room}:modifiers")
        modifiers |= r.hgetall("room:default:modifiers")
        # reconstruct them to get the schema
        classes = []
        for modifier in modifiers.values():
            modifier = json.loads(modifier)
            cls = get_cls_from_json_schema(modifier["schema"], modifier["name"])
            classes.append(cls)

        emit(
            "modifier:schema",
            Modifier.get_updated_schema(extensions=classes),
            to=request.sid,
        )

    @io.on("draw:schema")
    def draw_schema():
        return Geometry.get_updated_schema()

    @io.on("scene:schema")
    def scene_schema():
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")
        config = znsocket.Dict(r, f"room:{room}:config")
        try:
            scene = Scene(**config["scene"])
            emit("scene:schema", get_schema_with_instance_defaults(scene), to=room)
        except KeyError:
            emit("scene:schema", Scene.get_updated_schema(), to=request.sid)

    @io.on("selection:schema")
    def selection_schema():
        emit("selection:schema", Selection.get_updated_schema(), to=request.sid)

    @io.on("analysis:schema")
    def analysis_schema():
        emit("analysis:schema", Analysis.get_updated_schema(), to=request.sid)

    @io.on("geometry:schema")
    def geometry_schema():
        emit("geometry:schema", Geometry.get_updated_schema(), to=request.sid)

    @io.on("modifier:run")
    def modifier_run(data: dict):
        room = session.get("token")
        r: Redis = current_app.extensions["redis"]

        name = data["method"]["discriminator"]

        public = r.smembers(f"room:default:modifiers:{name}")
        privat = r.smembers(f"room:{room}:modifiers:{name}")

        data["ZNDRAW_CLIENT_ROOM"] = room

        queue_position = 1

        if len(public):
            # The modifier was registered with public=True
            queue_position = r.rpush(f"modifier:queue:{name}", json.dumps(data))
        elif len(privat):
            # The modifier was registered with public=False
            queue_position = r.rpush(f"modifier:queue:{room}:{name}", json.dumps(data))
        else:
            # This would be the queue for default modifiers.
            # but they are queued using celery directly.
            # so no need for redis queue.
            pass

        emit("room:modifier:queue", queue_position, to=room)

        clients: set[str] = public | privat
        if len(clients):
            for sid in clients:
                # kindly ask every client if they are available
                emit("modifier:wakeup", to=sid)
        else:
            run_modifier.delay(room, data)

    @io.on("room:modifier:queue")
    def room_modifier_run(data: int):
        """Forwarding finished message."""
        room = session.get("token")
        emit("room:modifier:queue", data, to=room)

    @io.on("modifier:available")
    def modifier_available(modifier_names: list[str]) -> None:
        """Update state of registered modifier classes."""
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")  # TODO: Why use get, token should always be set

        for name in modifier_names:
            public = r.smembers(f"room:default:modifiers:{name}")
            privat = r.smembers(f"room:{room}:modifiers:{name}")

            if request.sid in public:
                # TODO: !!! there will be an issue if you register a privat modifier "a"
                # and then someone else registers a public modifier "a"
                public_task = r.lpop(f"modifier:queue:{name}")
                if public_task:
                    log.debug(f"running public task {public_task} on {request.sid}")
                    emit("modifier:run", json.loads(public_task), to=request.sid)
                    return

            if request.sid in privat:
                privat_task = r.lpop(f"modifier:queue:{room}:{name}")
                if privat_task:
                    log.debug(f"running private task {privat_task} on {request.sid}")
                    emit("modifier:run", json.loads(privat_task), to=request.sid)
                    return

        log.debug(f"No task available for {modifier_names}")

    @io.on("analysis:run")
    def analysis_run(data: dict):
        room = session.get("token")
        emit("room:analysis:queue", 1, to=room)  # TODO: find the correct queue position
        run_analysis.delay(room, data)

    @io.on("room:analysis:queue")
    def room_analysis_run(data: int):
        """Forwarding finished message."""
        room = session.get("token")
        emit("room:analysis:queue", data, to=room)

    @io.on("room:geometry:queue")
    def room_geometry_run(data: int):
        """Forwarding finished message."""
        room = session.get("token")
        emit("room:geometry:queue", data, to=room)

    @io.on("geometry:run")
    def geometry_run(data: dict):
        room = session.get("token")
        emit("geometry:run:enqueue", to=room)
        run_geometry.delay(room, data)

    # @io.on("room:geometry:set")
    # def room_geometry_set(data: list | None = None):
    #     r: Redis = current_app.extensions["redis"]
    #     room = session.get("token")

    #     lst = znsocket.List(r, f"room:{room}:geometries")

    #     if data is not None:
    #         del lst[:]
    #         lst.extend(data)
    #     print(f"room:geometry:set: {lst}")

    #     emit(
    #         "room:geometry:set",
    #         [x.model_dump() for x in lst],
    #         to=room,
    #         include_self=False,
    #     )

    # @io.on("room:geometry:get")
    # def room_geometry_get():
    #     r: Redis = current_app.extensions["redis"]
    #     room = session.get("token")
    #     # TODO: do not send the data but request it on the JS side
    #     return [x.model_dump() for x in znsocket.List(r, f"room:{room}:geometries")]

    @io.on("room:config:get")
    def room_config_get():
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")
        data = znsocket.Dict(r, f"room:{room}:config")
        if set(data.keys()) != {"arrows", "scene"}:
            from zndraw.config import ZnDrawConfig

            data.update(ZnDrawConfig(vis=None).to_dict())
        return dict(data)

    @io.on("room:config:set")
    def room_config_set(data: dict):
        print(f"room:config:set: {data}")
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")
        config = znsocket.Dict(r, f"room:{room}:config")
        config.update(data)
        emit("room:config:set", data, to=room, include_self=False)
        scene = Scene(**config["scene"])
        emit("scene:schema", get_schema_with_instance_defaults(scene), to=room)

    @io.on("analysis:figure:set")
    def analysis_figure_set(data: dict):
        # This is currently using push and the figure is not stored
        room = session.get("token")
        r: Redis = current_app.extensions["redis"]
        if not r.exists(f"room:{room}:analysis:figures"):
            default_figures = znsocket.Dict(r, "room:default:analysis:figures")
            figures = znsocket.Dict(r, f"room:{room}:analysis:figures")
            figures.update(default_figures)
        else:
            figures = znsocket.Dict(r, f"room:{room}:analysis:figures")
        # remove all keys that are not in data
        for key in list(figures.keys()):
            if key not in data:
                del figures[key]
        # add all keys that are in data
        for key, value in data.items():
            figures[key] = value
        emit(
            "analysis:figure:refresh",
            list(data.keys()),
            to=room,
        )

    @io.on("analysis:figure:get")
    def analysis_figure_get(key: str | None = None) -> dict | None:
        room = session.get("token")
        r: Redis = current_app.extensions["redis"]
        # check if room exists, else copy from default
        if not r.exists(f"room:{room}:analysis:figures"):
            default_figures = znsocket.Dict(r, "room:default:analysis:figures")
            figures = znsocket.Dict(r, f"room:{room}:analysis:figures")
            figures.update(default_figures)
        else:
            figures = znsocket.Dict(r, f"room:{room}:analysis:figures")
        if key is None:
            return dict(figures)
        return figures.get(key)

    @io.on("analysis:figure:keys")
    def analysis_figure_keys() -> list[str]:
        room = session.get("token")
        r: Redis = current_app.extensions["redis"]
        if not r.exists(f"room:{room}:analysis:figures"):
            return znsocket.Dict(r, "room:default:analysis:figures").keys()
        else:
            return znsocket.Dict(r, f"room:{room}:analysis:figures").keys()

    @io.on("selection:run")
    def selection_run(data: dict):
        room = session.get("token")
        emit("room:selection:queue", 1, to=room)  # TODO: find the correct queue position
        run_selection.delay(room, data)

    # @io.on("selection:run:finished")
    # def selection_run_finished():
    #     """Forwarding finished message."""
    #     room = session.get("token")
    #     emit("selection:run:finished", to=room)

    # @io.on("selection:run:running")
    # def selection_run_running():
    #     """Forwarding running method."""
    #     room = session.get("token")
    #     emit("selection:run:running", to=room)

    @io.on("room:selection:queue")
    def room_selection_run(data: int):
        """Forwarding finished message."""
        room = session.get("token")
        emit("room:selection:queue", data, to=room)

    @io.on("room:alert")
    def room_alert(msg: str):
        """Forward the alert message to every client in the room"""
        # TODO: identify the source client.
        room = session.get("token")
        emit("room:alert", msg, to=room)

    @io.on("room:selection:set")
    def room_selection_set(data: dict[str, list[int]]):
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")
        r.hset(
            f"room:{room}:selection", mapping={k: json.dumps(v) for k, v in data.items()}
        )
        emit("room:selection:set", data, to=room, include_self=False)

    @io.on("room:selection:get")
    def room_selection_get() -> dict[str, list[int]]:
        r: Redis = current_app.extensions["redis"]
        room = session.get("token")
        result = r.hgetall(f"room:{room}:selection")
        if "0" in result:
            return {k: json.loads(v) for k, v in result.items()}
        return {"0": []}

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

    @io.on("room:frames:refresh")
    def room_frames_refresh(frames: list[int]):
        room = session.get("token")
        emit("room:frames:refresh", frames, to=room)

    @io.on("room:token:get")
    def room_token_get() -> str:
        return session.get("token")
