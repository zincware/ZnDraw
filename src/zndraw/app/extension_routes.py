"""Extension and worker management routes.

Handles worker status tracking and extension discovery.
Extension registration is handled via Socket.IO in events.py.
"""

import json
import logging

from flask import Blueprint, Response, current_app, jsonify, request

from zndraw.auth import get_current_user, require_auth
from zndraw.server import socketio

from .constants import SocketEvents
from .redis_keys import ExtensionKeys, JobKeys, RoomKeys
from .worker_stats import WorkerStats

log = logging.getLogger(__name__)

extensions = Blueprint("extensions", __name__)


def _submit_extension_impl(
    category: str, extension: str, public: bool, room_id: str
) -> tuple[Response, int]:
    """Shared implementation for extension submit endpoints.

    Parameters
    ----------
    category : str
        Extension category (e.g., "modifiers", "analysis")
    extension : str
        Extension name
    public : bool
        Whether this is a global (public=True) or room-scoped (public=False) extension
    room_id : str
        Target room ID

    Returns
    -------
    tuple[dict, int]
        JSON response and HTTP status code
    """
    from .job_manager import JobManager
    from .queue_manager import emit_queue_update

    # Get authenticated user
    user_name = get_current_user()

    # Extract request data
    json_data = request.json
    if json_data is None:
        json_data = {}

    data = json_data.pop("data", None)
    if data is None:
        # If no "data" key, treat the entire JSON body as data
        data = json_data

    log.info(
        f"Extension submit: room={room_id}, category={category}, extension={extension}, "
        f"public={public}, user={user_name}, data={json.dumps(data)}"
    )

    redis_client = current_app.extensions["redis"]

    # Check if this is a server-side (Celery) extension
    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "analysis": analysis,
        "settings": settings,
    }

    is_celery_extension = (
        category in category_map and extension in category_map[category]
    )

    # Celery extensions must use public endpoint
    if is_celery_extension:
        if not public:
            return jsonify({
                "error": f"Server-side extension '{extension}' must use public endpoint",
                "code": "WRONG_ENDPOINT",
                "extension": extension,
                "expectedEndpoint": "public",
            }), 400
        # Celery extensions are inherently global
    else:
        # Validate client extension exists in correct namespace
        if public:
            keys = ExtensionKeys.for_global_extension(category, extension)
            schema_key = ExtensionKeys.global_schema_key(category)
            fs_type = "global"
        else:
            keys = ExtensionKeys.for_extension(room_id, category, extension)
            schema_key = ExtensionKeys.schema_key(room_id, category)
            fs_type = "room-scoped"

        extension_schema = redis_client.hget(schema_key, extension)
        if extension_schema is None:
            return jsonify({
                "error": f"{fs_type.capitalize()} extension '{extension}' not found",
                "code": "EXTENSION_NOT_FOUND",
                "extension": extension,
                "namespace": fs_type,
            }), 404

    # Store the extension data
    room_keys = RoomKeys(room_id)
    redis_client.hset(
        room_keys.user_extension_data(user_name, category),
        extension,
        json.dumps(data),
    )

    # Handle settings differently - no job queue, just update and notify
    if category == "settings":
        socketio.emit(
            SocketEvents.INVALIDATE,
            {
                "userName": user_name,
                "category": category,
                "extension": extension,
                "roomId": room_id,
            },
            to=f"room:{room_id}",
        )
        log.info(
            f"Updated settings for room {room_id}, user {user_name}: "
            f"{extension} = {json.dumps(data)}"
        )
        return jsonify({"status": "success", "message": "Settings updated"}), 200

    # Create job
    provider = "celery" if is_celery_extension else "client"
    job_id = JobManager.create_job(
        redis_client, room_id, category, extension, data, user_name, provider, public
    )

    queue_position = 0

    if is_celery_extension:
        # Celery extensions use global keys (shared across all rooms)
        keys = ExtensionKeys.for_global_extension(category, extension)

        # Add to queue
        redis_client.rpush(
            keys.queue,
            json.dumps({
                "user_name": user_name,
                "data": data,
                "room": room_id,
                "jobId": job_id,
                "provider": "celery",
                "public": public,  # Always True for Celery extensions
            }),
        )
        log.info(
            f"Queued Celery task for user {user_name}, category {category}, "
            f"extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1

        # Trigger celery worker
        from zndraw.app.tasks import celery_job_worker

        config = current_app.extensions["config"]
        server_url = config.server_url
        _ = celery_job_worker.delay(room_id, server_url)

        # Notify clients in room
        emit_queue_update(redis_client, room_id, category, extension, socketio)
    else:
        # Client extensions: use public or room-scoped keys
        if public:
            keys = ExtensionKeys.for_global_extension(category, extension)
            log.info(
                f"Queuing job for global extension {extension} in category {category}"
            )
        else:
            keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add to queue
        redis_client.rpush(
            keys.queue,
            json.dumps({
                "user_name": user_name,
                "data": data,
                "room": room_id,
                "jobId": job_id,
                "public": public,
            }),
        )
        log.info(
            f"Queued task for user {user_name}, category {category}, "
            f"extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1

        # Notify clients (globally for public, room-scoped for private)
        emit_queue_update(
            redis_client,
            None if public else room_id,
            category,
            extension,
            socketio,
        )

    # Notify user
    log.info(
        f"Emitting invalidate for user {user_name}, category {category}, "
        f"extension {extension}, room {room_id} to user:{user_name}"
    )
    socketio.emit(
        SocketEvents.INVALIDATE,
        {
            "userName": user_name,
            "category": category,
            "extension": extension,
            "roomId": room_id,
        },
        to=f"user:{user_name}",
    )
    return jsonify({"status": "success", "queuePosition": queue_position, "jobId": job_id}), 200


@extensions.route(
    "/api/rooms/<room_id>/extensions/private/<category>/<extension>/submit",
    methods=["POST"],
)
@require_auth
def submit_private_extension(room_id: str, category: str, extension: str):
    """Submit a room-scoped extension action.

    Only accesses extensions registered with public=False for this specific room.
    """
    return _submit_extension_impl(category, extension, public=False, room_id=room_id)


@extensions.route(
    "/api/rooms/<room_id>/extensions/public/<category>/<extension>/submit",
    methods=["POST"],
)
@require_auth
def submit_public_extension(room_id: str, category: str, extension: str):
    """Submit a global extension action.

    Only accesses extensions registered with public=True.
    Includes server-side Celery extensions.
    """
    return _submit_extension_impl(category, extension, public=True, room_id=room_id)


@extensions.route(
    "/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/data",
    methods=["GET"],
)
def get_extension_data(room_id: str, category: str, extension: str):
    """Gets the cached data for a user's extension."""
    from zndraw.auth import AuthError, get_current_user

    # Authenticate and get user from JWT token
    try:
        user_name = get_current_user()
    except AuthError as e:
        return {"error": e.message}, e.status_code

    log.info(
        f"get_extension_data called with userName={user_name}, category={category}, extension={extension} for room {room_id}"
    )

    redis_client = current_app.extensions["redis"]
    room_keys = RoomKeys(room_id)
    extension_data = redis_client.hget(
        room_keys.user_extension_data(user_name, category), extension
    )
    if extension_data is None:
        return {"data": None}, 200
    extension_data = json.loads(extension_data)
    return {"data": extension_data}, 200


@extensions.route(
    "/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/workers",
    methods=["GET"],
)
def get_extension_workers(room_id: str, category: str, extension: str):
    """Get worker details for a specific extension.

    Returns:
        {
            "idleWorkers": ["worker_id_1", "worker_id_2"],
            "progressingWorkers": ["worker_id_3"],
            "queueLength": 5
        }
    """

    redis_client = current_app.extensions["redis"]
    keys = ExtensionKeys.for_extension(room_id, category, extension)

    idle_workers = list(redis_client.smembers(keys.idle_workers))
    progressing_workers = list(redis_client.smembers(keys.progressing_workers))
    queue_length = redis_client.llen(keys.queue)

    return {
        "idleWorkers": idle_workers,
        "progressingWorkers": progressing_workers,
        "queueLength": queue_length,
        "totalWorkers": len(idle_workers) + len(progressing_workers),
    }, 200


@extensions.route("/api/workers/<string:worker_id>", methods=["GET"])
def get_worker_state(worker_id: str):
    """Get the current state of a worker.

    Returns:
        {
            "idle": true/false,
            "currentJob": job_id or null
        }
    """
    redis_client = current_app.extensions["redis"]

    # Check if there's any job assigned to this worker
    # by looking for jobs where worker_id matches
    current_job_id = None
    is_idle = True

    for key in redis_client.scan_iter(match="job:*"):
        job_data = redis_client.hgetall(key)
        if (
            job_data.get("worker_id") == worker_id
            and job_data.get("status") == "running"
        ):
            current_job_id = job_data.get("id")
            is_idle = False
            break

    return {"idle": is_idle, "currentJob": current_job_id}, 200
