"""Extension and worker management routes.

Handles extension registration, worker status tracking, and extension discovery.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio

from .constants import SocketEvents
from .redis_keys import ExtensionKeys
from .worker_stats import WorkerStats

log = logging.getLogger(__name__)

extensions = Blueprint("extensions", __name__)


@extensions.route("/api/rooms/<string:room_id>/extensions/register", methods=["POST"])
def register_extension(room_id: str):
    data = request.get_json()
    try:
        name = data["name"]
        category = data["category"]
        schema = data["schema"]
        worker_id = data["userName"]  # Actually contains worker_id (sid) sent by client
    except KeyError as e:
        return {"error": f"Missing required field: {e}"}, 400

    redis_client = current_app.extensions["redis"]

    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
        "analysis": analysis,
    }

    if category in category_map and name in category_map[category]:
        log.warning(
            f"Blocked attempt to register extension '{name}' in category '{category}' "
            f"- name conflicts with server-side extension (security violation)"
        )
        return {
            "error": f"Cannot register extension '{name}': name is reserved for server-side extensions"
        }

    log.info(
        f"Registering extension for room {room_id}: name={name}, category={category}"
    )

    keys = ExtensionKeys.for_extension(room_id, category, name)
    worker_extensions_key = ExtensionKeys.user_extensions_key(
        room_id, category, worker_id
    )
    existing_schema = redis_client.hget(keys.schema, name)

    if existing_schema is not None:
        existing_schema = json.loads(existing_schema)
        if existing_schema != schema:
            return {
                "error": "Extension with this name already exists with a different schema"
            }, 400
        redis_client.sadd(keys.idle_workers, worker_id)
        redis_client.sadd(
            worker_extensions_key, name
        )  # Keep this for disconnect cleanup

        # Notify clients about worker count change
        log.info(
            f"Worker {worker_id} re-registered for extension '{name}' "
            f"in category '{category}', invalidating schema"
        )
        socketio.emit(
            SocketEvents.INVALIDATE_SCHEMA,
            {"roomId": room_id, "category": category},
            to=f"room:{room_id}",
        )

        # Check if there are queued tasks that this worker can handle
        # if room_id:  # Null check for type safety
        #     dispatch_next_task(redis_client, socketio, worker_id, room_id, category)

        return {
            "status": "success",
            "message": "Extension already registered with same schema. Worker marked as idle.",
        }
    else:
        # This is a brand new extension for the room
        with redis_client.pipeline() as pipe:
            # Set the schema
            pipe.hset(keys.schema, name, json.dumps(schema))
            # Add the current worker to the set of idle workers for this extension
            pipe.sadd(keys.idle_workers, worker_id)
            # Maintain the reverse mapping for easy cleanup on disconnect
            pipe.sadd(worker_extensions_key, name)
            pipe.execute()

        # Invalidate schema on all clients so they can see the new extension
        socketio.emit(
            SocketEvents.INVALIDATE_SCHEMA,
            {"roomId": room_id, "category": category},
            to=f"room:{room_id}",
        )

        # Check if there are queued tasks that this worker can handle
        # if room_id:  # Null check for type safety
        #     dispatch_next_task(redis_client, socketio, worker_id, room_id, category)

    return {"status": "success"}


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
        # Check if any client has registered this extension

        schema_key = ExtensionKeys.schema_key(room_id, category)
        extension_schema = redis_client.hget(schema_key, extension)

        if extension_schema is None:
            return {"error": f"No workers available for extension {extension}"}, 400

    # Store the entire extension data as a JSON string
    redis_client.hset(
        f"room:{room_id}:user:{user_name}:{category}", extension, json.dumps(data)
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

        server_url = current_app.config.get("SERVER_URL", "http://localhost:5000")
        _ = celery_job_worker.delay(room_id, server_url)

        # Notify all clients in room about queue update
        emit_queue_update(redis_client, room_id, category, extension, socketio)
    else:
        # Queue job for client workers to poll via /jobs/next endpoint
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

        # Notify all clients in room about queue update
        emit_queue_update(redis_client, room_id, category, extension, socketio)

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
    extension_data = redis_client.hget(
        f"room:{room_id}:user:{user_name}:{category}", extension
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

    # Check both modifiers and selections categories to find if worker is progressing
    current_job_id = None
    is_idle = True

    for category in ["modifiers", "selections", "analysis"]:
        # Get all extensions this worker might be registered for
        # We need to check all extensions to find if the worker is in any progressing set
        # First, get the worker's registered extensions
        user_extensions_key = ExtensionKeys.user_extensions_key(
            "*", category, worker_id
        )

        # Since we don't know the room, we need to scan for worker state
        # Check if worker has a current job by scanning all rooms
        # For simplicity, we'll check the most common pattern

        # Alternative approach: check if there's any job assigned to this worker
        # by looking for jobs where worker_id matches
        for key in redis_client.scan_iter(match="job:*"):
            job_data = redis_client.hgetall(key)
            if (
                job_data.get("worker_id") == worker_id
                and job_data.get("status") == "running"
            ):
                current_job_id = job_data.get("id")
                is_idle = False
                break

        if current_job_id:
            break

    return {"idle": is_idle, "currentJob": current_job_id}, 200


@extensions.route("/api/rooms/<string:room_id>/extensions/overview", methods=["GET"])
def get_room_extensions_overview(room_id: str):
    """Get comprehensive extension overview for a room.

    Query Parameters:
        category: Optional category filter ("modifiers", "selections", "analysis")
        search: Optional search string

    Returns:
        {
            extensions: [...],
            summary: {total_extensions, active_workers, total_jobs_24h}
        }
    """
    from datetime import datetime

    from dateutil.parser import isoparse

    from .analytics import aggregate_job_stats
    from .extension_utils import get_server_extensions

    redis_client = current_app.extensions["redis"]
    category_filter = request.args.get("category")
    search = request.args.get("search", "").lower()

    extensions = []
    total_workers = 0
    total_jobs_24h = 0

    # Get all extension categories to scan
    categories = (
        [category_filter]
        if category_filter
        else ["modifiers", "selections", "analysis"]
    )

    for category in categories:
        # Get active extensions (with schemas - can accept new jobs)
        schema_key = ExtensionKeys.schema_key(room_id, category)
        client_extensions = redis_client.hkeys(schema_key)
        client_extensions = [
            name.decode() if isinstance(name, bytes) else name
            for name in client_extensions
        ]

        # Get server-side extensions (always available)
        server_extensions = get_server_extensions(category)

        # Get historical extensions (disconnected clients with job history)
        # Scan for job keys: room:{room_id}:extension:{category}:*:jobs
        historical_extensions = set()
        pattern = f"room:{room_id}:extension:{category}:*:jobs"
        for key in redis_client.scan_iter(pattern):
            key_str = key.decode() if isinstance(key, bytes) else key
            # Parse: room:{room_id}:extension:{category}:{extension}:jobs
            parts = key_str.split(":")
            if len(parts) == 6 and parts[5] == "jobs":
                ext_name = parts[4]
                historical_extensions.add(ext_name)

        # Combine all three sets of extensions
        all_extension_names = (
            set(client_extensions) | server_extensions | historical_extensions
        )

        for ext_name in all_extension_names:
            # Search filter
            if search and search not in ext_name.lower():
                continue

            # Get schema (only exists for client-side extensions in Redis)
            schema_json = redis_client.hget(schema_key, ext_name)
            schema = json.loads(schema_json) if schema_json else {}

            # Get worker stats
            keys = ExtensionKeys.for_extension(room_id, category, ext_name)
            worker_stats = WorkerStats.fetch(redis_client, keys)
            total_workers += worker_stats.total_workers

            # Get jobs for this extension
            job_ids = redis_client.smembers(
                f"room:{room_id}:extension:{category}:{ext_name}:jobs"
            )
            jobs = []
            for job_id in job_ids:
                job_id_str = job_id.decode() if isinstance(job_id, bytes) else job_id
                job_data = redis_client.hgetall(f"job:{job_id_str}")
                if job_data:
                    jobs.append(
                        {
                            k.decode() if isinstance(k, bytes) else k: v.decode()
                            if isinstance(v, bytes)
                            else v
                            for k, v in job_data.items()
                        }
                    )

            # Aggregate stats
            stats = aggregate_job_stats(jobs)

            # Count last 24h jobs
            now = datetime.utcnow()
            jobs_24h = 0
            for j in jobs:
                if j.get("created_at"):
                    try:
                        if (now - isoparse(j["created_at"])).total_seconds() < 86400:
                            jobs_24h += 1
                    except Exception:
                        pass
            total_jobs_24h += jobs_24h

            # Get recent jobs (last 5)
            recent_jobs = sorted(
                [j for j in jobs if j.get("created_at")],
                key=lambda x: x.get("created_at", ""),
                reverse=True,
            )[:5]

            # Only include extensions with at least one job
            if stats.total_jobs == 0:
                continue

            # Determine provider (celery if server-side, client otherwise)
            server_extensions = get_server_extensions(category)
            provider = "celery" if ext_name in server_extensions else "client"

            extensions.append(
                {
                    "name": ext_name,
                    "category": category,
                    "provider": provider,
                    "schema": schema,
                    "workers": {
                        "idle_count": worker_stats.idle_count,
                        "progressing_count": worker_stats.progressing_count,
                        "queue_length": worker_stats.queue_length,
                    },
                    "analytics": {
                        "total_jobs": stats.total_jobs,
                        "success_rate": stats.success_rate,
                        "avg_wait_time_ms": stats.avg_wait_time_ms,
                        "avg_execution_time_ms": stats.avg_execution_time_ms,
                        "last_used": stats.last_used,
                    },
                    "recent_jobs": [
                        {
                            "id": j.get("id"),
                            "status": j.get("status"),
                            "created_at": j.get("created_at"),
                            "wait_time_ms": int(j["wait_time_ms"])
                            if j.get("wait_time_ms")
                            and j["wait_time_ms"] not in ("", "0")
                            else None,
                            "execution_time_ms": int(j["execution_time_ms"])
                            if j.get("execution_time_ms")
                            and j["execution_time_ms"] not in ("", "0")
                            else None,
                        }
                        for j in recent_jobs
                    ],
                }
            )

    return {
        "extensions": extensions,
        "summary": {
            "total_extensions": len(extensions),
            "active_workers": total_workers,
            "total_jobs_24h": total_jobs_24h,
        },
    }


@extensions.route("/api/extensions", methods=["GET"])
def get_global_extensions_overview():
    """Get global extensions overview across all rooms.

    Query Parameters:
        category: Optional category filter
        search: Optional search string

    Returns:
        {extensions: [...]}
    """
    from .analytics import aggregate_job_stats
    from .extension_utils import get_server_extensions

    redis_client = current_app.extensions["redis"]
    category_filter = request.args.get("category")
    search = request.args.get("search", "").lower()

    # Scan all extension job keys to find extensions with jobs
    # Pattern: room:{room_id}:extension:{category}:{extension}:jobs
    extensions_map = {}  # {category:extension: {rooms: set, jobs: []}}
    valid_categories = ["modifiers", "selections", "analysis"]

    for key in redis_client.scan_iter("room:*:extension:*:jobs"):
        key_str = key.decode() if isinstance(key, bytes) else key
        # Parse: room:{room_id}:extension:{category}:{extension}:jobs
        parts = key_str.split(":")
        if len(parts) != 6 or parts[2] != "extension" or parts[5] != "jobs":
            continue

        room_id = parts[1]
        category = parts[3]
        ext_name = parts[4]

        # Validate category
        if category not in valid_categories:
            continue

        if category_filter and category != category_filter:
            continue

        if search and search not in ext_name.lower():
            continue

        map_key = f"{category}:{ext_name}"
        if map_key not in extensions_map:
            extensions_map[map_key] = {"rooms": set(), "jobs": []}

        extensions_map[map_key]["rooms"].add(room_id)

        # Get jobs for this extension
        job_ids = redis_client.smembers(key)
        for job_id in job_ids:
            job_id_str = job_id.decode() if isinstance(job_id, bytes) else job_id
            job_data = redis_client.hgetall(f"job:{job_id_str}")
            if job_data:
                extensions_map[map_key]["jobs"].append(
                    {
                        k.decode() if isinstance(k, bytes) else k: v.decode()
                        if isinstance(v, bytes)
                        else v
                        for k, v in job_data.items()
                    }
                )

    # Aggregate global stats
    result = []
    for map_key, data in extensions_map.items():
        category, ext_name = map_key.split(":", 1)
        stats = aggregate_job_stats(data["jobs"])

        # Only include extensions with at least one job
        if stats.total_jobs == 0:
            continue

        server_extensions = get_server_extensions(category)
        provider = "celery" if ext_name in server_extensions else "client"

        result.append(
            {
                "name": ext_name,
                "category": category,
                "provider": provider,
                "rooms": list(data["rooms"]),
                "global_stats": {
                    "total_jobs": stats.total_jobs,
                    "avg_success_rate": stats.success_rate,
                    "avg_wait_time_ms": stats.avg_wait_time_ms,
                    "avg_execution_time_ms": stats.avg_execution_time_ms,
                },
            }
        )

    return {"extensions": result}


@extensions.route(
    "/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/analytics",
    methods=["GET"],
)
def get_extension_detailed_analytics(room_id: str, category: str, extension: str):
    """Get detailed analytics for a specific extension.

    Returns:
        {
            daily_stats: [{date, job_count, success_rate, avg_wait_ms, avg_exec_ms}],
            total_stats: {total_jobs, overall_success_rate, ...},
            error_breakdown: [{error, count}]
        }
    """
    from datetime import datetime

    from dateutil.parser import isoparse

    from .analytics import aggregate_job_stats, get_daily_stats

    redis_client = current_app.extensions["redis"]

    # Get all jobs for this extension
    job_ids = redis_client.smembers(
        f"room:{room_id}:extension:{category}:{extension}:jobs"
    )
    jobs = []
    for job_id in job_ids:
        job_id_str = job_id.decode() if isinstance(job_id, bytes) else job_id
        job_data = redis_client.hgetall(f"job:{job_id_str}")
        if job_data:
            jobs.append(
                {
                    k.decode() if isinstance(k, bytes) else k: v.decode()
                    if isinstance(v, bytes)
                    else v
                    for k, v in job_data.items()
                }
            )

    # Aggregate total stats
    total_stats = aggregate_job_stats(jobs)

    # Calculate days based on earliest job, or default to 7 if no jobs
    days = 7
    if jobs:
        try:
            # Find earliest created_at timestamp
            valid_timestamps = []
            for job in jobs:
                if job.get("created_at"):
                    try:
                        valid_timestamps.append(isoparse(job["created_at"]))
                    except Exception:
                        pass

            if valid_timestamps:
                earliest = min(valid_timestamps)
                now = datetime.utcnow()
                days_diff = (now - earliest).days + 1  # +1 to include today
                days = max(days_diff, 1)  # At least 1 day
        except Exception:
            days = 7

    # Get daily breakdown
    daily_stats = get_daily_stats(jobs, days)

    return {
        "daily_stats": daily_stats,
        "total_stats": {
            "total_jobs": total_stats.total_jobs,
            "overall_success_rate": total_stats.success_rate,
            "overall_avg_wait_ms": total_stats.avg_wait_time_ms,
            "overall_avg_execution_ms": total_stats.avg_execution_time_ms,
        },
        "error_breakdown": total_stats.error_breakdown,
    }
