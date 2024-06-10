import json
import logging
import uuid

import znsocket
from flask import current_app, request, session
from flask_socketio import emit, join_room
from redis import Redis

from zndraw.analyse import Analysis
from zndraw.draw import Geometry
from zndraw.modify import Modifier
from zndraw.scene import Scene
from zndraw.select import Selection
from zndraw.utils import get_cls_from_json_schema

from ..app import socketio as io
from ..tasks import run_analysis, run_geometry, run_modifier, run_selection

log = logging.getLogger(__name__)


@io.on("connect")
def connect():
    pass


@io.on("shutdown")
def shutdown():
    if current_app.config["AUTH_TOKEN"] is None or session["authenticated"]:
        log.critical("Shutting down server")
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
    r = current_app.config["redis"]

    if "name" in session:
        log.critical(f"disconnecting (webclient) {request.sid} from room {room}")
        r.srem(f"room:{room}:webclients", session["name"])
        emit("room:users:refresh", list(r.smembers(f"room:{room}:webclients")), to=room)
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

    room = str(session["token"])
    join_room(room)  # rename token to room or room_token

    session["name"] = uuid.uuid4().hex[:8]

    r = current_app.config["redis"]
    r.sadd(f"room:{room}:webclients", session["name"])

    emit("room:users:refresh", list(r.smembers(f"room:{room}:webclients")), to=room)
    # set step, camera, bookmarks, points

    log.critical(f"connecting (webclient) {request.sid} to {room}")

    emit("selection:schema", Selection.updated_schema())
    emit("modifier:schema:refresh", to=request.sid)
    emit("scene:schema", Scene.updated_schema())
    emit("geometry:schema", Geometry.updated_schema())
    emit("analysis:schema:refresh", to=request.sid)

    if current_app.config["TUTORIAL"] is not None:
        emit("tutorial:url", current_app.config["TUTORIAL"], to=request.sid)
    if current_app.config["SIMGEN"]:
        emit("showSiMGen", True, to=request.sid)

    return {
        "name": session["name"],
        "room": room,
        "needs_auth": current_app.config["AUTH_TOKEN"] is not None,
    }


@io.on("join")  # rename pyclient:connect
def join(data: dict):
    """
    Arguments:
        data: {"token": str, "auth_token": str}
    """
    # TODO: prohibit "token" to be "default"

    if current_app.config["AUTH_TOKEN"] is None:
        session["authenticated"] = True
    else:
        session["authenticated"] = data["auth_token"] == current_app.config["AUTH_TOKEN"]
    token = data["token"]
    session["token"] = token
    room = str(session["token"])

    join_room(room)
    log.critical(f"connecting (pyclient) {request.sid} to {room}")
    # join_room(f"pyclients_{token}")


@io.on("room:frames:get")
def room_frames_get(frames: list[int]) -> dict[int, dict]:
    # print(f"requesting frames: {frames}")
    if len(frames) == 0:
        return {}
    r: Redis = current_app.config["redis"]
    room = session.get("token")

    if r.exists(f"room:{room}:frames"):
        data = znsocket.List(r, f"room:{room}:frames")[frames]

    else:
        try:
            data = znsocket.List(r, "room:default:frames")[frames]
        except IndexError:
            data = []

    return {idx: json.loads(d) for idx, d in zip(frames, data) if d is not None}


@io.on("room:frames:set")
def room_frames_set(data: dict[int, str]):
    r: Redis = current_app.config["redis"]
    room = session.get("token")

    # add = {}
    # remove = []
    lst = znsocket.List(r, f"room:{room}:frames")

    if not r.exists(f"room:{room}:frames"):
        default_lst = znsocket.List(r, "room:default:frames")
        # TODO: using a redis copy action would be faster
        lst.extend(default_lst)

    data = {int(k): v for k, v in sorted(data.items())}
    for key, value in sorted(data.items()):
        if int(key) == len(lst):
            lst.append(value)
        elif int(key) > len(lst):
            log.critical(f"frame {key} is out of bounds and is being ignored.")
            pass
        else:
            lst[int(key)] = value

    emit("room:frames:refresh", [int(x) for x in data], to=room)


@io.on("room:all:frames:refresh")
def room_all_frames_refresh(indices: list[int]):
    emit("room:frames:refresh", [int(x) for x in indices], broadcast=True)


@io.on("room:frames:delete")
def room_frames_delete(frames: list[int]):
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    lst = znsocket.List(r, f"room:{room}:frames")
    if not r.exists(f"room:{room}:frames"):
        default_lst = znsocket.List(r, "room:default:frames")
        # TODO: using a redis copy action would be faster
        lst.extend(default_lst)
    del lst[frames]
    # TODO how to update here?
    emit("room:frames:refresh", [int(x) for x in frames], to=room)


@io.on("room:frames:insert")
def room_frames_insert(data: dict):
    index = data.pop("index")
    value = data.pop("value")
    r: Redis = current_app.config["redis"]
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


@io.on("room:length:get")
def room_frames_length_get() -> int:
    room = session.get("token")
    r: Redis = current_app.config["redis"]
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

    r: Redis = current_app.config["redis"]
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
    r: Redis = current_app.config["redis"]
    room = session.get("token")

    modifiers: dict = r.hgetall(f"room:{room}:modifiers")
    modifiers |= r.hgetall("room:default:modifiers")
    # reconstruct them to get the schema
    classes = []
    for modifier in modifiers.values():
        modifier = json.loads(modifier)
        cls = get_cls_from_json_schema(modifier["schema"], modifier["name"])
        classes.append(cls)

    emit("modifier:schema", Modifier.updated_schema(extensions=classes), to=request.sid)


@io.on("draw:schema")
def draw_schema():
    return Geometry.updated_schema()


@io.on("scene:schema")
def scene_schema():
    emit("scene:schema", Scene.updated_schema(), to=request.sid)


@io.on("selection:schema")
def selection_schema():
    emit("selection:schema", Selection.updated_schema(), to=request.sid)


@io.on("analysis:schema")
def analysis_schema():
    emit("analysis:schema", Analysis.updated_schema(), to=request.sid)


@io.on("modifier:run")
def modifier_run(data: dict):
    room = session.get("token")
    emit("room:modifier:queue", 1, to=room)  # TODO: find the correct queue position
    url = f"http://127.0.0.1:{current_app.config['PORT']}"

    r: Redis = current_app.config["redis"]

    name = data["method"]["discriminator"]

    public = r.smembers(f"room:default:modifiers:{name}")
    privat = r.smembers(f"room:{room}:modifiers:{name}")

    data["ZNDRAW_CLIENT_ROOM"] = room

    if len(public):
        # The modifier was registered with public=True
        r.rpush(f"modifier:queue:{name}", json.dumps(data))
    elif len(privat):
        # The modifier was registered with public=False
        r.rpush(f"modifier:queue:{room}:{name}", json.dumps(data))
    else:
        # This would be the queue for default modifiers.
        # but they are queued using celery directly.
        # so no need for redis queue.
        pass

    clients: set[str] = public | privat
    if len(clients):
        for sid in clients:
            # kindly ask every client if they are available
            emit("modifier:wakeup", to=sid)
    else:
        run_modifier.delay(url, room, data)


@io.on("room:modifier:queue")
def room_modifier_run(data: int):
    """Forwarding finished message."""
    room = session.get("token")
    emit("room:modifier:queue", data, to=room)


@io.on("modifier:available")
def modifier_available(modifier_names: list[str]) -> None:
    """Update state of registered modifier classes."""
    r: Redis = current_app.config["redis"]
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
    url = f"http://127.0.0.1:{current_app.config['PORT']}"
    run_analysis.delay(url, room, data)


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
    url = f"http://127.0.0.1:{current_app.config['PORT']}"
    run_geometry.delay(url, room, data)


@io.on("room:geometry:set")
def room_geometry_set(data: list):
    r: Redis = current_app.config["redis"]
    room = session.get("token")

    # # add = {}
    # # remove = []
    lst = znsocket.List(r, f"room:{room}:geometries")
    del lst[:]
    lst.extend(data)
    emit("room:geometry:set", data, to=room, include_self=False)


@io.on("room:geometry:get")
def room_geometry_get():
    r: Redis = current_app.config["redis"]
    room = session.get("token")

    return list(znsocket.List(r, f"room:{room}:geometries"))


@io.on("analysis:figure:set")
def analysis_figure_set(data: dict):
    # This is currently using push and the figure is not stored
    room = session.get("token")
    r: Redis = current_app.config["redis"]
    r.set(f"room:{room}:analysis:figure", json.dumps(data))
    emit(
        "analysis:figure:set",
        data,
        to=room,
    )


@io.on("analysis:figure:get")
def analysis_figure_get() -> dict:
    room = session.get("token")
    r: Redis = current_app.config["redis"]
    return json.loads(r.get(f"room:{room}:analysis:figure"))


@io.on("selection:run")
def selection_run(data: dict):
    room = session.get("token")
    emit("room:selection:queue", 1, to=room)  # TODO: find the correct queue position
    url = f"http://127.0.0.1:{current_app.config['PORT']}"
    run_selection.delay(url, room, data)


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


@io.on("room:log")
def room_log(msg: str):
    """Forward the alert message to every client in the room"""
    room = session.get("token")
    emit("room:log", msg, to=room)


@io.on("room:selection:set")
def room_selection_set(data: dict[str, list[int]]):
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    r.hmset(f"room:{room}:selection", {k: json.dumps(v) for k, v in data.items()})
    emit("room:selection:set", data, to=room, include_self=False)


@io.on("room:selection:get")
def room_selection_get() -> dict[str, list[int]]:
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    result = r.hgetall(f"room:{room}:selection")
    if "0" in result:
        return {k: json.loads(v) for k, v in result.items()}
    return {"0": []}


@io.on("room:step:set")
def room_step_set(step: int):
    print(f"setting step to {step}")
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    r.set(f"room:{room}:step", step)

    emit("room:step:set", step, to=room, include_self=False)


@io.on("room:step:get")
def room_step_get() -> int:
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    return r.get(f"room:{room}:step") or 0


@io.on("room:points:set")
def room_points_set(data: dict):
    print(f"setting points to {data}")
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    r.hmset(f"room:{room}:points", {k: json.dumps(v) for k, v in data.items()})

    emit("room:points:set", data, to=room, include_self=False)
    # TODO: add rotation! save position and rotation and scale?


@io.on("room:points:get")
def room_points_get() -> dict[str, list[list]]:
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    result: dict[int, list[list]] = r.hgetall(f"room:{room}:points")
    # TODO: type consistency!!
    result = {k: json.loads(v) for k, v in result.items()}
    if "0" not in result:
        return {"0": []}
    return result


@io.on("room:bookmarks:set")
def room_bookmarks_set(data: dict):
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    r.hmset(f"room:{room}:bookmarks", data)
    emit("room:bookmarks:set", data, to=room, include_self=False)


@io.on("room:bookmarks:get")
def room_bookmarks_get() -> dict:
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    return r.hgetall(f"room:{room}:bookmarks")


@io.on("room:camera:set")
def room_camera_set(data: dict):
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    r.set(f"room:{room}:camera", json.dumps(data))
    emit("room:camera:set", data, to=room, include_self=False)


@io.on("room:camera:get")
def room_camera_get() -> dict:
    r: Redis = current_app.config["redis"]
    room = session.get("token")
    return json.loads(r.get(f"room:{room}:camera"))
