import logging
from threading import Lock

from flask import request, session
from flask_socketio import emit

from zndraw.server import tasks
from zndraw.utils import typecast

from ..app import cache
from ..app import socketio as io
from .data import (
    CeleryTaskData,
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
    return cache.get("DEFAULT_PYCLIENT")


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
        print("Submitting jobs .....................")
        tasks.get_selection_schema.delay(request.url_root, request.sid)
        tasks.read_file.delay(request.url_root, request.sid)
    except KeyError:
        pass


@io.on("celery:task:results")
@typecast
def celery_task_results(msg: CeleryTaskData):
    emit(msg.event, msg.data, to=msg.target)

    # DEFAULT_PYCLIENT = cache.get("DEFAULT_PYCLIENT")
    # if DEFAULT_PYCLIENT is None and "token" in session:
    #     # refuse connection if there is no default pyclient
    #     return False
    # with contextlib.suppress(KeyError):
    #     # If you connect through Python, you don't have a token.
    #     token = session["token"]
    #     join_room(f"webclients_{token}")
    #     # who ever connected latest is the HOST of the room
    #     ROOM_HOSTS = cache.get("ROOM_HOSTS")
    #     try:
    #         ROOM_HOSTS[token].append(request.sid)
    #     except KeyError:
    #         ROOM_HOSTS[token] = [request.sid]
    #     cache.set("ROOM_HOSTS", ROOM_HOSTS)

    #     data = {"sid": request.sid, "token": token}
    #     data["host"] = ROOM_HOSTS[token][0] == request.sid
    #     names = cache.get(f"PER-TOKEN-NAME:{session['token']}") or {}
    #     names[request.sid] = uuid4().hex[:8].upper()
    #     cache.set(f"PER-TOKEN-NAME:{session['token']}", names)

    #     emit("webclient:available", data, to=DEFAULT_PYCLIENT)

    #     connected_users = [{"name": names[sid]} for sid in ROOM_HOSTS[token]]

    #     emit(
    #         "connectedUsers",
    #         list(reversed(connected_users)),
    #         to=_webclients_room({"token": token}),
    #     )

    #     data = {"modifiers": []}  # {schema: ..., name: ...}
    #     MODIFIER = cache.get("MODIFIER")
    #     for name, schema in MODIFIER["default_schema"].items():
    #         data["modifiers"].append({"schema": schema, "name": name})
    #     data["token"] = token

    #     emit("modifier:register", data, to=DEFAULT_PYCLIENT)

    #     # TODO emit("modifier:register", _all modifiers_, to=app.config["DEFAULT_PYCLIENT"]')

    #     log.debug(f"connected {request.sid} and updated HOSTS to {ROOM_HOSTS}")
    #     emit("message:log", "Connection established", to=request.sid)
    #     PER_TOKEN_DATA = cache.get("PER-TOKEN-DATA")
    #     if token not in PER_TOKEN_DATA:
    #         PER_TOKEN_DATA[token] = {}
    #     cache.set("PER-TOKEN-DATA", PER_TOKEN_DATA)

    #     # append to zndraw.log a line isoformat() + " " + token
    #     if "token" not in app.config:
    #         with open("zndraw.log", "a") as f:
    #             f.write(
    #                 datetime.datetime.now().isoformat() + " " + token + " connected \n"
    #             )
