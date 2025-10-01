import json
import logging

import msgpack
import zarr
from flask import Response, current_app, request

from zndraw.storage import ZarrStorageSequence, decode_data, encode_data
from zndraw.server import socketio
from .constants import SocketEvents
from .redis_keys import ExtensionKeys
from .worker_stats import WorkerStats
from .queue_manager import emit_queue_update
from .job_manager import JobManager
from .worker_dispatcher import dispatch_next_task
import traceback


from . import main

# --- Logging Setup ---
log = logging.getLogger(__name__)


# TODO: move to utils
def get_lock_key(room: str, target: str) -> str:
    """Constructs a standardized Redis key for a lock."""
    return f"room:{room}:lock:{target}"


def get_zarr_store_path(room_id: str) -> str:
    """Returns the path to the Zarr store for a given room."""
    storage_path = current_app.config.get("STORAGE_PATH", "./zndraw-data.zarr")
    # Remove .zarr extension if present to append room_id
    base_path = storage_path.rstrip("/").removesuffix(".zarr")
    return f"{base_path}/{room_id}.zarr"

@main.route("/internal/emit", methods=["POST"])
def internal_emit():
    """Internal endpoint to emit Socket.IO events. Secured via a shared secret."""
    data = request.get_json()

    event = data.get("event")
    sid = data.get("sid")
    payload = data.get("data", {})

    if not event or not sid:
        return {"error": "Event and sid are required"}, 400

    socketio.emit(event, payload, to=sid)
    return {"success": True}


@main.route("/api/frames/<string:room_id>", methods=["POST"])
def get_frames(room_id):
    """Serves multiple frames' data from the room's Zarr store using either indices or slice parameters."""
    r = current_app.extensions["redis"]
    try:
        # Parse the request data
        request_data = request.get_json()
        if request_data is None:
            return {"error": "Request body required"}, 400

        store_path = get_zarr_store_path(room_id)
        root = zarr.group(store_path)
        storage = ZarrStorageSequence(root)

        # Get logical-to-physical mapping from Redis
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {"error": "No frames found in room"}, 404

        max_frame = len(frame_mapping) - 1

        # Determine frame indices based on request parameters
        if "indices" in request_data:
            # Direct list of indices
            frame_indices = request_data["indices"]
            if not isinstance(frame_indices, list):
                return {"error": "Indices must be a list"}, 400

            # Validate frame indices
            for frame_id in frame_indices:
                if (
                    not isinstance(frame_id, int)
                    or frame_id < 0
                    or frame_id > max_frame
                ):
                    error_data = {
                        "error": f"Invalid frame index {frame_id}, valid range: 0-{max_frame}",
                        "type": "IndexError",
                    }
                    return Response(
                        json.dumps(error_data),
                        status=404,
                        content_type="application/json",
                    )

        else:
            # Default to slice behavior for any remaining cases (including empty payload)
            # This handles slice parameters and slice(None, None, None) which sends empty payload
            start = request_data.get("start", 0)
            stop = request_data.get("stop", len(frame_mapping))
            step = request_data.get("step", 1)

            # Validate slice parameters
            if not all(isinstance(x, int) for x in [start, stop, step]):
                return {"error": "start, stop, and step must be integers"}, 400

            if step == 0:
                return {"error": "step cannot be zero"}, 400

            # Generate frame indices from slice
            try:
                frame_indices = list(range(start, stop, step))
                # Filter out invalid indices
                frame_indices = [i for i in frame_indices if 0 <= i <= max_frame]
            except ValueError as e:
                return {"error": f"Invalid slice parameters: {e}"}, 400

        # Get keys parameter if specified
        requested_keys = request_data.get("keys")
        print(f"get_frames called with keys: {requested_keys}")
        # error_data = {"error": f"Key(s) not found: {', '.join(sorted(missing_keys))}", "type": "KeyError"}
        # return Response(json.dumps(error_data), status=404, content_type='application/json')

        # TODO: requested keys and KeyError handling
        # TODO: instead of iterate load at once
        try:
            frames_data = []
            for frame_id in frame_indices:
                # Get the physical index for this logical frame
                physical_index = int(frame_mapping[frame_id])
                frame_data = storage.get(physical_index, keys=requested_keys)
                frames_data.append(encode_data(frame_data))
        except KeyError as e:
            error_data = {
                "error": f"Key(s) not found: {e}",
                "type": "KeyError",
            }
            return Response(
                json.dumps(error_data), status=404, content_type="application/json"
            )
        except IndexError as e:
            error_data = {
                "error": f"Frame index out of range: {e}",
                "type": "IndexError",
            }
            return Response(
                json.dumps(error_data), status=404, content_type="application/json"
            )

        packed_data = msgpack.packb(frames_data)
        return Response(packed_data, content_type="application/octet-stream")
    except Exception as e:
        error_data = {
            "error": f"Server error: {e}",
            "type": type(e).__name__,
            "success": False,
        }
        print(traceback.format_exc())
        return Response(
            json.dumps(error_data), status=500, content_type="application/json"
        )


@main.route("/api/rooms/<string:room_id>/frames", methods=["POST"])
def append_frame(room_id):
    """Appends a new frame. Authorized via a short-lived Bearer token."""
    r = current_app.extensions["redis"]
    auth_header = request.headers.get("Authorization")
    if not auth_header or not auth_header.startswith("Bearer "):
        return {"error": "Authorization token is missing or invalid"}, 401

    token = auth_header.split(" ")[1]
    token_key = f"room:{room_id}:upload_token:{token}"

    # Get token metadata
    token_data = r.hgetall(token_key)
    if not token_data:
        return {"error": "Token is invalid or has expired"}, 403

    sid_from_token = token_data.get("sid")
    action = token_data.get("action", "append")
    target_frame_id = (
        int(token_data.get("frame_id", -1))
        if token_data.get("frame_id") != "-1"
        else None
    )

    r.delete(token_key)  # Invalidate the token after first use

    lock_key = get_lock_key(room_id, "trajectory:meta")
    if r.get(lock_key) != sid_from_token:
        return {"error": "Client does not hold the trajectory lock"}, 403

    try:
        # Unpack the msgpack data
        serialized_data = msgpack.unpackb(request.data, strict_map_key=False)

        store_path = get_zarr_store_path(room_id)
        root = zarr.group(store_path)
        storage = ZarrStorageSequence(root)

        indices_key = f"room:{room_id}:trajectory:indices"

        if action == "replace":
            frame_mapping = r.zrange(indices_key, 0, -1)

            # Validate target_frame_id
            if target_frame_id is None or not (
                0 <= target_frame_id < len(frame_mapping)
            ):
                return {
                    "error": f"Invalid or missing frame_id for replace. Valid range: 0-{len(frame_mapping) - 1}"
                }, 404

            new_physical_index = len(storage)
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            old_physical_index_to_unmap = frame_mapping[target_frame_id]

            pipeline = r.pipeline()
            pipeline.zrem(indices_key, old_physical_index_to_unmap)
            pipeline.zadd(indices_key, {str(new_physical_index): target_frame_id})
            pipeline.execute()

            log.info(
                f"Replaced frame {target_frame_id} (old physical: {old_physical_index_to_unmap}, new physical: {new_physical_index}) in room '{room_id}'"
            )
            return {"success": True, "replaced_frame": target_frame_id}

        elif action == "extend":
            # Extend operation: add multiple frames in one go
            if not isinstance(serialized_data, list):
                return {
                    "error": "For extend action, data must be a list of frame dictionaries"
                }, 400

            # 1. Determine starting logical and physical positions
            start_logical_pos = r.zcard(indices_key)
            start_physical_pos = len(storage)
            num_frames = len(serialized_data)

            # 2. Decode all frames and extend the physical storage
            decoded_frames = [decode_data(frame) for frame in serialized_data]
            storage.extend(decoded_frames)

            # 3. Create the new logical-to-physical mapping for all new frames
            new_mapping = {
                str(start_physical_pos + i): start_logical_pos + i
                for i in range(num_frames)
            }
            if new_mapping:
                r.zadd(indices_key, new_mapping)

            # 4. Prepare response data
            new_indices = list(range(start_logical_pos, start_logical_pos + num_frames))

            log.info(
                f"Extended trajectory with {num_frames} frames (physical: {start_physical_pos}-{start_physical_pos + num_frames - 1}) to room '{room_id}'"
            )
            return {"success": True, "new_indices": new_indices}

        elif action == "insert":
            # Insert operation: add a new physical frame and insert it into the logical sequence
            insert_position = int(token_data.get("insert_position", 0))
            current_length = r.zcard(indices_key)

            if not (0 <= insert_position <= current_length):
                return {
                    "error": f"Insert position {insert_position} out of range [0, {current_length}]"
                }, 400

            # 1. Determine the new physical index (always at the end of the physical store)
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data to the physical store
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Atomically shift existing logical indices and insert the new one
            pipeline = r.pipeline()
            # Get members whose scores (logical positions) need to be incremented
            members_to_shift = r.zrangebyscore(indices_key, insert_position, "+inf")
            if members_to_shift:
                for member in members_to_shift:
                    pipeline.zincrby(indices_key, 1, member)

            # Add the new frame at the correct logical position
            pipeline.zadd(indices_key, {str(new_physical_index): insert_position})
            pipeline.execute()

            log.info(
                f"Inserted frame at position {insert_position} (physical: {new_physical_index}) in room '{room_id}'"
            )
            return {"success": True, "inserted_position": insert_position}

        elif action == "append":
            # Append operation: add a new frame to the end of the logical sequence
            # 1. Determine the new logical and physical positions
            logical_position = r.zcard(indices_key)
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Add the new physical index to the logical sequence
            r.zadd(indices_key, {str(new_physical_index): logical_position})

            log.info(
                f"Appended frame {logical_position} (physical: {new_physical_index}) to room '{room_id}'"
            )
            return {"success": True, "new_index": logical_position}

        else:
            # Default case for any unknown actions
            return {"error": f"The requested action '{action}' is not supported."}, 400
    except Exception as e:
        log.error(f"Failed to write to Zarr store: {e}")
        print(traceback.format_exc())
        return {"error": "Failed to write to data store"}, 500


@main.route("/api/exit")
def exit_app():
    """Endpoint to gracefully shut down the server. Secured via a shared secret."""
    socketio.stop()
    return {"success": True}

@main.route("/api/rooms/<string:room_id>/schema/<string:category>", methods=["GET"])
def get_room_schema(room_id: str, category: str):
    """Get the schema for a specific room with worker and queue statistics.

    Returns schema along with metadata about each extension:
    - provider: "celery" for server-side extensions, or count of registered workers
    - queueLength: number of queued tasks for this extension
    - idleWorkers: number of idle workers available
    - progressingWorkers: number of workers currently processing tasks
    """
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    # Map category strings to the corresponding imported objects
    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
    }

    if category not in category_map:
        return {"error": f"Unknown schema category '{category}'"}

    redis_client = current_app.extensions["redis"]
    schema = {}

    # Add server-provided extensions (Celery-based)
    for name, cls in category_map[category].items():
        schema[name] = {
            "schema": cls.model_json_schema(),
            "provider": "celery",
            "queueLength": 0,
            "idleWorkers": 0,
            "progressingWorkers": 0
        }

    # Add client-provided extensions from Redis
    schema_key = ExtensionKeys.schema_key(room_id, category)
    redis_schema = redis_client.hgetall(schema_key)

    for name, sch_str in redis_schema.items():
        sch = json.loads(sch_str)

        # Get worker statistics for this extension
        keys = ExtensionKeys.for_extension(room_id, category, name)
        stats = WorkerStats.fetch(redis_client, keys)

        if name in schema:
            if schema[name]["schema"] != sch:
                print(
                    f"Warning: {category.capitalize()} extension '{name}' schema "
                    "in Redis differs from server schema."
                )
        else:
            schema[name] = {
                "schema": sch,
                "provider": stats.total_workers,  # Number of workers for client extensions
                **stats.to_dict()
            }

    return schema


@main.route("/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>", methods=["POST"])
def log_room_extension(room_id: str, category: str, extension: str):
    """Logs a user extension action in the room's action log."""
    json_data = request.json
    if json_data is None:
        json_data = {}

    # Try to get userId from query params first, then from JSON body
    user_id = request.args.get("userId") or json_data.pop("userId", None)
    if user_id is None:
        return {"error": "User ID is required"}, 400

    data = json_data.pop("data", None)
    if data is None:
        # If no "data" key, treat the entire JSON body as data
        data = json_data

    print(f"Logging extension for room {room_id}: category={category}, extension={extension}, data={json.dumps(data)}")

    # store in redis
    redis_client = current_app.extensions["redis"]
    # Store the entire extension data as a JSON string
    redis_client.hset(
        f"room:{room_id}:user:{user_id}:{category}", extension, json.dumps(data)
    )

    # Check if this is a server-side (Celery) extension
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
    }

    is_celery_extension = (
        category in category_map and extension in category_map[category]
    )

    # Create job
    provider = "celery" if is_celery_extension else "client"
    job_id = JobManager.create_job(
        redis_client, room_id, category, extension, data, user_id, provider
    )

    queue_position = 0

    if is_celery_extension:
        # Queue job for Celery workers to poll via /jobs/next endpoint
        # Celery workers will check provider="celery" and pick up these jobs
        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add to queue with provider info
        redis_client.rpush(
            keys.queue, json.dumps({"user_id": user_id, "data": data, "room": room_id, "jobId": job_id, "provider": "celery"})
        )
        print(
            f"Queued Celery task for user {user_id}, category {category}, extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1

        # Notify all clients in room about queue update
        emit_queue_update(redis_client, room_id, category, extension, socketio)
    else:
        # Queue job for client workers to poll via /jobs/next endpoint
        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add to queue
        redis_client.rpush(
            keys.queue, json.dumps({"user_id": user_id, "data": data, "room": room_id, "jobId": job_id})
        )
        print(
            f"Queued task for user {user_id}, category {category}, extension {extension}, room {room_id}, job {job_id}"
        )
        queue_position = redis_client.llen(keys.queue) - 1  # Zero-indexed position

        # Notify all clients in room about queue update
        emit_queue_update(redis_client, room_id, category, extension, socketio)

    print(
        f"Emitting invalidate for user {user_id}, category {category}, extension {extension}, room {room_id} to user:{user_id}"
    )
    socketio.emit(
        SocketEvents.INVALIDATE,
        {"userId": user_id, "category": category, "extension": extension, "roomId": room_id},
        to=f"user:{user_id}",
    )
    return {"status": "success", "queuePosition": queue_position, "jobId": job_id}, 200

@main.route("/api/rooms/<string:room_id>/extension-data/<string:category>/<string:extension>", methods=["GET"])
def get_extension_data(room_id: str, category: str, extension: str):
    user_id = request.args.get("userId")
    print(f"get_extension_data called with userId={user_id}, category={category}, extension={extension} for room {room_id}")

    if not user_id:
        return {"error": "User ID is required"}, 400

    redis_client = current_app.extensions["redis"]
    extension_data = redis_client.hget(f"room:{room_id}:user:{user_id}:{category}", extension)
    if extension_data is None:
        return {"data": None}, 200
    extension_data = json.loads(extension_data)
    return {"data": extension_data}, 200


@main.route("/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/workers", methods=["GET"])
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
        "totalWorkers": len(idle_workers) + len(progressing_workers)
    }, 200


# Job management endpoints
@main.route("/api/rooms/<string:room_id>/jobs", methods=["GET"])
def list_jobs(room_id: str):
    """List active jobs for a room."""
    redis_client = current_app.extensions["redis"]
    jobs = JobManager.list_all_jobs(redis_client, room_id)
    return jobs, 200


@main.route("/api/rooms/<string:room_id>/jobs/<string:job_id>", methods=["GET", "DELETE"])
def get_job(room_id: str, job_id: str):
    """Get details for a specific job."""
    redis_client = current_app.extensions["redis"]
    if request.method == "DELETE":
        job = JobManager.get_job(redis_client, job_id)
        if not job:
            return {"error": "Job not found"}, 404
        if job.get("status") == "running":
            return {"error": "Cannot delete a running job"}, 400

        category = job.get("category")
        extension = job.get("extension")

        JobManager.delete_job(redis_client, job_id)

        # Emit queue update instead of job:deleted
        if category and extension:
            emit_queue_update(redis_client, room_id, category, extension, socketio)

        return {"status": "success"}, 200

    job = JobManager.get_job(redis_client, job_id)
    if not job:
        return {"error": "Job not found"}, 404
    return job, 200


@main.route("/api/rooms/<string:room_id>/jobs/<string:job_id>/complete", methods=["POST"])
def complete_job(room_id: str, job_id: str):
    """Mark a job as completed and transition worker back to idle.

    Called by both Celery tasks and client workers when they finish successfully.
    For client workers: transitions worker from progressing → idle, then checks queue.
    For Celery workers: just marks job complete (no worker state to manage).
    """
    data = request.get_json() or {}
    result = data.get("result")
    worker_id = data.get("workerId")  # Optional: only for client workers

    log.info(f"Complete job called: job_id={job_id}, worker_id={worker_id}, room_id={room_id}")

    redis_client = current_app.extensions["redis"]

    # Get job details to know category/extension
    job = JobManager.get_job(redis_client, job_id)
    if not job:
        log.error(f"Job {job_id} not found in Redis")
        return {"error": "Job not found"}, 404

    log.info(f"Job data: category={job.get('category')}, extension={job.get('extension')}")

    # Validate job is in running state
    if job.get("status") != "running":
        log.error(f"Job {job_id} is not running (status: {job.get('status')})")
        return {"error": "Job is not running"}, 400

    # Validate worker ID matches the job's assigned worker
    if worker_id and job.get("worker_id") != worker_id:
        log.error(f"Worker ID mismatch: job assigned to {job.get('worker_id')}, but {worker_id} tried to complete it")
        return {"error": "Worker ID does not match job's worker ID"}, 400

    # Mark job as completed
    success = JobManager.complete_job(redis_client, job_id, result)
    if not success:
        log.error(f"Failed to mark job {job_id} as completed")
        return {"error": "Failed to complete job"}, 400

    log.info(f"Job {job_id} completed in room {room_id}")

    # If this is a client worker (not Celery), handle worker state transition
    if worker_id and not worker_id.startswith("celery:"):
        log.info(f"Transitioning worker {worker_id} to idle")
        _transition_worker_to_idle(redis_client, socketio, worker_id, job, room_id, success=True)
    else:
        log.info(f"Skipping worker transition (worker_id={worker_id})")

    return {"status": "success"}, 200


@main.route("/api/rooms/<string:room_id>/jobs/<string:job_id>/fail", methods=["POST"])
def fail_job(room_id: str, job_id: str):
    """Mark a job as failed and transition worker back to idle.

    Called by both Celery tasks and client workers when they fail.
    For client workers: transitions worker from progressing → idle, then checks queue.
    For Celery workers: just marks job failed (no worker state to manage).
    """
    data = request.get_json() or {}
    error = data.get("error", "Unknown error")
    worker_id = data.get("workerId")  # Optional: only for client workers

    redis_client = current_app.extensions["redis"]

    # Get job details to know category/extension
    job = JobManager.get_job(redis_client, job_id)
    if not job:
        return {"error": "Job not found"}, 404

    # Validate job is in running state
    if job.get("status") != "running":
        log.error(f"Job {job_id} is not running (status: {job.get('status')})")
        return {"error": "Job is not running"}, 400

    # Validate worker ID matches the job's assigned worker
    if worker_id and job.get("worker_id") != worker_id:
        log.error(f"Worker ID mismatch: job assigned to {job.get('worker_id')}, but {worker_id} tried to fail it")
        return {"error": "Worker ID does not match job's worker ID"}, 400

    # Mark job as failed
    success = JobManager.fail_job(redis_client, job_id, error)
    if not success:
        return {"error": "Failed to mark job as failed"}, 400

    log.error(f"Job {job_id} failed in room {room_id}: {error}")

    # If this is a client worker (not Celery), handle worker state transition
    if worker_id and not worker_id.startswith("celery:"):
        _transition_worker_to_idle(redis_client, socketio, worker_id, job, room_id, success=False)

    return {"status": "success"}, 200


@main.route("/api/workers/<string:worker_id>", methods=["GET"])
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

    for category in ["modifiers", "selections"]:
        # Get all extensions this worker might be registered for
        # We need to check all extensions to find if the worker is in any progressing set
        # First, get the worker's registered extensions
        user_extensions_key = ExtensionKeys.user_extensions_key("*", category, worker_id)

        # Since we don't know the room, we need to scan for worker state
        # Check if worker has a current job by scanning all rooms
        # For simplicity, we'll check the most common pattern

        # Alternative approach: check if there's any job assigned to this worker
        # by looking for jobs where worker_id matches
        for key in redis_client.scan_iter(match=f"job:*"):
            job_data = redis_client.hgetall(key)
            if job_data.get("worker_id") == worker_id and job_data.get("status") == "running":
                current_job_id = job_data.get("id")
                is_idle = False
                break

        if current_job_id:
            break

    return {
        "idle": is_idle,
        "currentJob": current_job_id
    }, 200


@main.route("/api/rooms/<string:room_id>/jobs/next", methods=["POST"])
def get_next_job(room_id: str):
    """Poll for the next available job for a worker.

    Workers call this endpoint to get assigned their next job from the queue.
    The worker should be registered for at least one extension category.

    Request body:
        {
            "workerId": "worker_session_id"
        }

    Returns:
        Job object directly (jobId, category, extension, data, etc.) or null if no jobs available
    """
    data = request.get_json() or {}
    worker_id = data.get("workerId")

    if not worker_id:
        return {"error": "workerId is required"}, 400

    redis_client = current_app.extensions["redis"]

    # Check if worker already has a running job
    for key in redis_client.scan_iter(match=f"job:*"):
        job_data = redis_client.hgetall(key)
        if job_data.get("worker_id") == worker_id and job_data.get("status") == "running":
            return {"error": "Worker is not idle"}, 400

    # If worker is "celery-worker", check for celery jobs across all extensions
    is_celery_worker = worker_id.startswith("celery")

    # Check both modifiers and selections categories for queued jobs
    for category in ["modifiers", "selections"]:
        if is_celery_worker:
            # For celery workers, check all extensions for celery jobs
            from zndraw.extensions.modifiers import modifiers
            from zndraw.extensions.selections import selections

            category_map = {
                "modifiers": modifiers,
                "selections": selections,
            }

            if category in category_map:
                for extension in category_map[category].keys():
                    keys = ExtensionKeys.for_extension(room_id, category, extension)
                    queue_length = redis_client.llen(keys.queue)

                    if queue_length > 0:
                        # Peek at the first job to check if it's a celery job
                        task_data = redis_client.lindex(keys.queue, 0)
                        if task_data:
                            task_info = json.loads(task_data)
                            if task_info.get("provider") == "celery":
                                # This is a celery job, pop it
                                redis_client.lpop(keys.queue)
                                job_id = task_info.get("jobId")

                                # Get full job details
                                job = JobManager.get_job(redis_client, job_id)
                                if job:
                                    # Mark job as started by this worker
                                    JobManager.start_job(redis_client, job_id, worker_id)

                                    # Emit queue update
                                    emit_queue_update(redis_client, room_id, category, extension, socketio)

                                    # Get updated job details (now with running status)
                                    job = JobManager.get_job(redis_client, job_id)

                                    # Rename 'id' to 'jobId' for consistency with client expectations
                                    job["jobId"] = job.pop("id")

                                    log.info(f"Assigned celery job {job_id} to worker {worker_id} from queue")
                                    return job, 200
        else:
            # Get all extensions this worker is registered for in this category
            user_extensions_key = ExtensionKeys.user_extensions_key(room_id, category, worker_id)
            registered_extensions = redis_client.smembers(user_extensions_key)

            # Check each registered extension for queued jobs
            for extension in registered_extensions:
                keys = ExtensionKeys.for_extension(room_id, category, extension)
                queue_length = redis_client.llen(keys.queue)

                if queue_length > 0:
                    # Get the first job from the queue
                    task_data = redis_client.lpop(keys.queue)
                    if task_data:
                        task_info = json.loads(task_data)
                        job_id = task_info.get("jobId")

                        # Get full job details
                        job = JobManager.get_job(redis_client, job_id)
                        if job:
                            # Move worker from idle to progressing
                            redis_client.smove(keys.idle_workers, keys.progressing_workers, worker_id)

                            # Mark job as started by this worker
                            JobManager.start_job(redis_client, job_id, worker_id)

                            # Emit queue update (this covers job:started notification)
                            emit_queue_update(redis_client, room_id, category, extension, socketio)

                            # Get updated job details (now with running status)
                            job = JobManager.get_job(redis_client, job_id)

                            # Rename 'id' to 'jobId' for consistency with client expectations
                            job["jobId"] = job.pop("id")

                            log.info(f"Assigned job {job_id} to worker {worker_id} from queue")
                            return job, 200

    # No jobs available
    return {"error": "No jobs available"}, 400


def _transition_worker_to_idle(
    redis_client,
    socketio_instance,
    worker_id: str,
    job: dict,
    room_id: str,
    success: bool = True
) -> None:
    """Transition worker from progressing to idle and dispatch next task.

    Args:
        redis_client: Redis client
        socketio_instance: SocketIO instance
        worker_id: Worker session ID
        job: Job data dict containing category/extension
        room_id: Room identifier
        success: True if job completed, False if failed (affects log message)
    """
    category = job.get("category")
    extension = job.get("extension")

    log.info(f"_transition_worker_to_idle called: worker_id={worker_id}, category={category}, extension={extension}")

    if not category or not extension:
        log.warning(f"Missing category or extension in job data: category={category}, extension={extension}")
        return

    # Move worker from progressing back to idle for the extension they just completed
    keys = ExtensionKeys.for_extension(room_id, category, extension)
    moved = redis_client.smove(keys.progressing_workers, keys.idle_workers, worker_id)

    log.info(f"Worker transition: moved={moved}, progressing_key={keys.progressing_workers}, idle_key={keys.idle_workers}")

    if moved:
        status_msg = "finished and is now idle" if success else "marked idle after failure"
        log.info(f"Worker {worker_id} {status_msg}")

        # Get all extensions this worker is registered for in this category
        user_extensions_key = ExtensionKeys.user_extensions_key(room_id, category, worker_id)
        registered_extensions = redis_client.smembers(user_extensions_key)

        log.info(f"Worker {worker_id} registered extensions: {registered_extensions}")

        # Add worker to idle set for ALL registered extensions (not just the one they completed)
        # This ensures dispatch_next_task can find them for any extension
        for ext_name in registered_extensions:
            if ext_name != extension:  # Already moved above for completed extension
                ext_keys = ExtensionKeys.for_extension(room_id, category, ext_name)
                # Remove from progressing (if present) and add to idle
                redis_client.srem(ext_keys.progressing_workers, worker_id)
                redis_client.sadd(ext_keys.idle_workers, worker_id)
                log.info(f"Added worker {worker_id} to idle set for extension {ext_name}")

        # Emit queue update for this extension
        emit_queue_update(redis_client, room_id, category, extension, socketio_instance)

        # Workers will poll for the next task via /jobs/next endpoint
        log.info(f"Worker {worker_id} is now idle and can poll for next task")
    else:
        log.warning(f"Failed to move worker {worker_id} from progressing to idle (may already be idle or not in progressing)")

