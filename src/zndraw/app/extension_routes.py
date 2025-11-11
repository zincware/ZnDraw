"""Extension and worker management routes.

Handles worker status tracking and extension discovery.
Extension registration is handled via Socket.IO in events.py.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio

from .constants import SocketEvents
from .redis_keys import ExtensionKeys, JobKeys, RoomKeys
from .worker_stats import WorkerStats

log = logging.getLogger(__name__)

extensions = Blueprint("extensions", __name__)


@extensions.route(
    "/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/submit",
    methods=["POST"],
)
def log_room_extension(room_id: str, category: str, extension: str):
    """Submits a user extension action to create a job."""
    from zndraw.auth import AuthError, get_current_user

    from .job_manager import JobManager
    from .queue_manager import emit_queue_update

    # Authenticate and get user from JWT token
    try:
        user_name = get_current_user()
    except AuthError as e:
        return {"error": e.message}, e.status_code

    json_data = request.json
    if json_data is None:
        json_data = {}

    data = json_data.pop("data", None)
    if data is None:
        # If no "data" key, treat the entire JSON body as data
        data = json_data

    log.info(
        f"Logging extension for room {room_id}: category={category}, extension={extension}, data={json.dumps(data)}"
    )

    # store in redis
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

    # Check if extension exists (either as server-side or client-registered)
    if not is_celery_extension:
        # Check if any client has registered this extension (room-scoped or global)

        # First check room-scoped extension
        schema_key = ExtensionKeys.schema_key(room_id, category)
        extension_schema = redis_client.hget(schema_key, extension)

        # If not found in room, check global extensions
        is_global = False
        if extension_schema is None:
            global_schema_key = ExtensionKeys.global_schema_key(category)
            extension_schema = redis_client.hget(global_schema_key, extension)
            is_global = extension_schema is not None

        if extension_schema is None:
            return {"error": f"No workers available for extension {extension}"}, 400

    # Store the entire extension data as a JSON string
    room_keys = RoomKeys(room_id)
    redis_client.hset(
        room_keys.user_extension_data(user_name, category), extension, json.dumps(data)
    )

    # Handle settings differently - no job queue, just update and notify
    if category == "settings":
        # Emit socket event to notify all clients in room
        socketio.emit(
            SocketEvents.INVALIDATE,
            {
                "userName": user_name,
                "category": category,
                "extension": extension,
                "roomId": room_id,
            },
            to=f"room:{room_id}",  # Broadcast to all users in room
        )
        log.info(
            f"Updated settings for room {room_id}, user {user_name}: {extension} = {json.dumps(data)}"
        )
        return {"status": "success", "message": "Settings updated"}, 200

    # Create job

    provider = "celery" if is_celery_extension else "client"
    job_id = JobManager.create_job(
        redis_client, room_id, category, extension, data, user_name, provider
    )

    queue_position = 0

    if is_celery_extension:
        # Queue job for Celery workers to poll via /jobs/next endpoint
        # Celery workers will check provider="celery" and pick up these jobs
        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add to queue with provider info
        redis_client.rpush(
            keys.queue,
            json.dumps(
                {
                    "user_name": user_name,
                    "data": data,
                    "room": room_id,
                    "jobId": job_id,
                    "provider": "celery",
                }
            ),
        )
        log.info(
            f"Queued Celery task for user {user_name}, category {category}, extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1

        # Trigger a celery worker task to pick up the job
        from zndraw.app.tasks import celery_job_worker

        config = current_app.extensions["config"]
        server_url = config.server_url
        _ = celery_job_worker.delay(room_id, server_url)

        # Notify all clients in room about queue update
        emit_queue_update(redis_client, room_id, category, extension, socketio)
    else:
        # Queue job for client workers to poll via /jobs/next endpoint
        # Use global or room-scoped keys based on whether extension is global
        if is_global:
            keys = ExtensionKeys.for_global_extension(category, extension)
            log.info(
                f"Queuing job for global extension {extension} in category {category}"
            )
        else:
            keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add to queue
        redis_client.rpush(
            keys.queue,
            json.dumps(
                {"user_name": user_name, "data": data, "room": room_id, "jobId": job_id}
            ),
        )
        log.info(
            f"Queued task for user {user_name}, category {category}, extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1  # Zero-indexed position

        # Notify all clients about queue update
        # For global extensions, pass None as room_id to broadcast globally
        emit_queue_update(
            redis_client,
            None if is_global else room_id,
            category,
            extension,
            socketio,
        )

    log.info(
        f"Emitting invalidate for user {user_name}, category {category}, extension {extension}, room {room_id} to user:{user_name}"
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
    return {"status": "success", "queuePosition": queue_position, "jobId": job_id}, 200


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
