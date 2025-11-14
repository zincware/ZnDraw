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
    from datetime import datetime

    from .job_dispatcher import assign_pending_jobs_for_extension
    from .job_manager import JobManager, JobStatus

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

    # Store the extension data (per-user)
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

    # Create job with PENDING status
    provider = "celery" if is_celery_extension else "client"
    job_id = JobManager.create_job(
        redis_client,
        room_id,
        category,
        extension,
        data,
        user_name,
        provider,
        public,
        initial_status=JobStatus.PENDING,
        socketio=socketio,  # Pass socketio for automatic event emission
    )

    # Determine keys based on public/private
    if public or is_celery_extension:
        keys = ExtensionKeys.for_global_extension(category, extension)
        ext_room_id = None  # Global extension
    else:
        keys = ExtensionKeys.for_extension(room_id, category, extension)
        ext_room_id = room_id  # Room-scoped extension

    # Add job to pending queue
    timestamp = datetime.utcnow().timestamp()
    redis_client.zadd(keys.pending_jobs, {job_id: timestamp})
    log.info(f"Added job {job_id} to pending queue for {category}/{extension}")

    # Check for available workers and try to assign immediately
    from .job_dispatcher import get_available_workers

    available_workers = get_available_workers(redis_client, keys)

    if available_workers:
        # Available worker exists - trigger assignment
        log.info(
            f"Available workers for {category}/{extension}, "
            f"attempting immediate assignment"
        )

        assigned = assign_pending_jobs_for_extension(
            redis_client,
            socketio,
            ext_room_id,
            category,
            extension,
        )

        if assigned > 0:
            log.info(f"Job {job_id} assigned to worker immediately")
            queue_position = 0  # Assigned, not in queue
        else:
            # Still in pending queue - position is number of jobs ahead of this one
            queue_position = redis_client.zcard(keys.pending_jobs) - 1
    else:
        # No available workers - job stays in pending queue
        log.info(
            f"No available workers for {category}/{extension}, "
            f"job {job_id} remains pending"
        )
        # Position is number of jobs ahead of this one
        queue_position = redis_client.zcard(keys.pending_jobs) - 1

    # For Celery extensions, trigger task
    if is_celery_extension:
        from zndraw.app.tasks import celery_job_worker

        config = current_app.extensions["config"]
        server_url = config.server_url

        # Dispatch Celery task and get task ID
        job_data = JobManager.get_job(redis_client, job_id)
        celery_task = celery_job_worker.delay(job_data, server_url)

        # Assign job to Celery worker (PENDING → ASSIGNED)
        # Worker ID format: "celery:{task_id}"
        worker_id = f"celery:{celery_task.id}"
        success = JobManager.assign_job(redis_client, job_id, worker_id, socketio=socketio)

        if success:
            log.info(
                f"Assigned job {job_id} to Celery worker {worker_id} and triggered task"
            )
            queue_position = 0  # Assigned, not in queue
        else:
            log.error(
                f"Failed to assign job {job_id} to Celery worker {worker_id}"
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
            "queueLength": 5,
            "totalWorkers": 3
        }
    """

    redis_client = current_app.extensions["redis"]
    keys = ExtensionKeys.for_extension(room_id, category, extension)

    # Get all workers registered for this extension
    all_workers = redis_client.hkeys(keys.workers)

    # Separate idle vs busy workers by checking capacity
    idle_workers = []
    busy_workers = []

    if all_workers:
        pipe = redis_client.pipeline()
        for worker_id in all_workers:
            capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
            pipe.get(capacity_key)
        capacities = pipe.execute()

        for worker_id, capacity in zip(all_workers, capacities):
            capacity_val = int(capacity) if capacity else 0
            if capacity_val >= 1:
                idle_workers.append(worker_id)
            else:
                busy_workers.append(worker_id)

    pending_jobs_count = redis_client.zcard(keys.pending_jobs)

    return {
        "idleWorkers": idle_workers,
        "progressingWorkers": busy_workers,
        "queueLength": pending_jobs_count,
        "totalWorkers": len(idle_workers) + len(busy_workers),
    }, 200


@extensions.route("/api/workers/<string:worker_id>", methods=["GET"])
@require_auth
def get_worker_state(worker_id: str):
    """Get the current state of a worker.

    Returns:
        {
            "idle": true/false,
            "currentJob": job_id or null
        }
    """
    redis_client = current_app.extensions["redis"]

    # Check worker's active jobs using the reverse index
    active_jobs_key = ExtensionKeys.worker_active_jobs_key(worker_id)
    active_job_ids = redis_client.smembers(active_jobs_key)

    # Worker is idle if they have no active jobs
    is_idle = len(active_job_ids) == 0
    current_job_id = None

    if active_job_ids:
        # Return the first active job (there should typically be only one)
        current_job_id = list(active_job_ids)[0]

    return {"idle": is_idle, "currentJob": current_job_id}, 200
