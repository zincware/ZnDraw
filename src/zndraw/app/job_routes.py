"""Job management routes.

Handles job listing, status updates, and worker job polling.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio

from .job_manager import JobManager
from .queue_manager import emit_queue_update
from .redis_keys import ExtensionKeys

log = logging.getLogger(__name__)

jobs = Blueprint("jobs", __name__)


def _transition_worker_to_idle(
    redis_client,
    socketio_instance,
    worker_id: str,
    job: dict,
    room_id: str,
    success: bool = True,
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

    log.info(
        f"_transition_worker_to_idle called: worker_id={worker_id}, category={category}, extension={extension}"
    )

    if not category or not extension:
        log.warning(
            f"Missing category or extension in job data: category={category}, extension={extension}"
        )
        return

    # Move worker from progressing back to idle for the extension they just completed
    keys = ExtensionKeys.for_extension(room_id, category, extension)
    moved = redis_client.smove(keys.progressing_workers, keys.idle_workers, worker_id)

    log.info(
        f"Worker transition: moved={moved}, progressing_key={keys.progressing_workers}, idle_key={keys.idle_workers}"
    )

    if moved:
        status_msg = (
            "finished and is now idle" if success else "marked idle after failure"
        )
        log.info(f"Worker {worker_id} {status_msg}")

        # Get all extensions this worker is registered for in this category
        user_extensions_key = ExtensionKeys.user_extensions_key(
            room_id, category, worker_id
        )
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
                log.info(
                    f"Added worker {worker_id} to idle set for extension {ext_name}"
                )

        # Emit queue update for this extension
        emit_queue_update(redis_client, room_id, category, extension, socketio_instance)

        # Workers will poll for the next task via /jobs/next endpoint
        log.info(f"Worker {worker_id} is now idle and can poll for next task")
    else:
        log.warning(
            f"Failed to move worker {worker_id} from progressing to idle (may already be idle or not in progressing)"
        )


@jobs.route("/api/rooms/<string:room_id>/jobs", methods=["GET"])
def list_jobs(room_id: str):
    """List active jobs for a room."""
    redis_client = current_app.extensions["redis"]
    jobs = JobManager.list_all_jobs(redis_client, room_id)
    return jobs, 200


@jobs.route(
    "/api/rooms/<string:room_id>/jobs/<string:job_id>", methods=["GET", "DELETE"]
)
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


@jobs.route("/api/rooms/<string:room_id>/jobs/<string:job_id>/status", methods=["PUT"])
def update_job_status(room_id: str, job_id: str):
    """Update a job's status (complete or fail) and transition worker back to idle.

    Called by both Celery tasks and client workers when they finish.
    For client workers: transitions worker from progressing → idle, then checks queue.
    For Celery workers: just marks job complete/failed (no worker state to manage).

    Request body:
        {
            "status": "completed" | "failed",
            "result": <result data>,  // for completed status
            "error": <error message>,  // for failed status
            "workerId": <worker_id>  // optional, only for client workers
        }
    """
    data = request.get_json() or {}
    status = data.get("status")
    worker_id = data.get("workerId")  # Optional: only for client workers

    if status not in ["completed", "failed"]:
        return {"error": "Status must be 'completed' or 'failed'"}, 400

    redis_client = current_app.extensions["redis"]

    # Get job details to know category/extension
    job = JobManager.get_job(redis_client, job_id)
    if not job:
        log.error(f"Job {job_id} not found in Redis")
        return {"error": "Job not found"}, 404

    log.info(
        f"Update job status called: job_id={job_id}, status={status}, worker_id={worker_id}, room_id={room_id}"
    )
    log.info(
        f"Job data: category={job.get('category')}, extension={job.get('extension')}"
    )

    # Validate job is in running state
    if job.get("status") != "running":
        log.error(f"Job {job_id} is not running (status: {job.get('status')})")
        return {"error": "Job is not running"}, 400

    # Validate worker ID matches the job's assigned worker
    if worker_id and job.get("worker_id") != worker_id:
        log.error(
            f"Worker ID mismatch: job assigned to {job.get('worker_id')}, but {worker_id} tried to update it"
        )
        return {"error": "Worker ID does not match job's worker ID"}, 400

    # Update job status based on requested status
    if status == "completed":
        result = data.get("result")
        success = JobManager.complete_job(redis_client, job_id, result)
        if not success:
            log.error(f"Failed to mark job {job_id} as completed")
            return {"error": "Failed to complete job"}, 400
        log.info(f"Job {job_id} completed in room {room_id}")
        worker_success = True
    else:  # status == "failed"
        error = data.get("error", "Unknown error")
        success = JobManager.fail_job(redis_client, job_id, error)
        if not success:
            return {"error": "Failed to mark job as failed"}, 400
        log.error(f"Job {job_id} failed in room {room_id}: {error}")

        # Log error to client UI via SocketIO event (similar to Celery workers)
        try:
            from zndraw.app.chat import send_message_to_room

            extension_name = job.get("extension", "unknown")
            error_message = f"❌ Error in {extension_name}: {error}"
            send_message_to_room(
                redis_client, socketio, room_id, error_message, worker_id or "system"
            )
        except Exception as e:
            log.warning(f"Failed to log error to client UI: {e}")

        worker_success = False

    # If this is a client worker (not Celery), handle worker state transition
    if worker_id and not worker_id.startswith("celery:"):
        log.info(f"Transitioning worker {worker_id} to idle")
        _transition_worker_to_idle(
            redis_client, socketio, worker_id, job, room_id, success=worker_success
        )
    else:
        log.info(f"Skipping worker transition (worker_id={worker_id})")

    return {"status": "success"}, 200


@jobs.route("/api/jobs/next", methods=["POST"])
def get_next_job_agnostic():
    """Poll for the next available job for a worker (room-agnostic).

    Workers call this endpoint to get assigned their next job from any queue
    (room-specific or global) where they have registered extensions.

    Request body:
        {
            "workerId": "worker_session_id"
        }

    Returns:
        Job object directly (jobId, room, category, extension, data, etc.) or error if no jobs available
    """
    try:
        data = request.get_json() or {}
        worker_id = data.get("workerId")

        if not worker_id:
            return {"error": "workerId is required"}, 400

        redis_client = current_app.extensions["redis"]

        # Check if worker already has a running job
        for key in redis_client.scan_iter(match="job:*"):
            job_data = redis_client.hgetall(key)
            if (
                job_data.get("worker_id") == worker_id
                and job_data.get("status") == "running"
            ):
                return {"error": "Worker is not idle"}, 400

        # Special handling for Celery workers
        is_celery_worker = worker_id.startswith("celery")

        # Check both modifiers, selections, and analysis categories for queued jobs
        for category in ["modifiers", "selections", "analysis"]:
            # If worker is Celery, check for celery jobs across all rooms
            if is_celery_worker:
                try:
                    from zndraw.extensions.analysis import analysis
                    from zndraw.extensions.modifiers import modifiers
                    from zndraw.extensions.selections import selections
                except ImportError as e:
                    log.error(f"Failed to import extensions: {e}")
                    continue

                category_map = {
                    "modifiers": modifiers,
                    "selections": selections,
                    "analysis": analysis,
                }

                if category in category_map:
                    # Celery extensions use global keys (shared across all rooms)
                    for extension in category_map[category].keys():
                        # Check global queue for this Celery extension
                        keys = ExtensionKeys.for_global_extension(category, extension)
                        queue_length = redis_client.llen(keys.queue)

                        if queue_length > 0:
                            # Peek at first job to check if it's a celery job
                            task_data = redis_client.lindex(keys.queue, 0)
                            if task_data:
                                try:
                                    task_info = json.loads(task_data)
                                except json.JSONDecodeError as e:
                                    log.error(
                                        f"Invalid JSON in queue {keys.queue}: {e}, data: {task_data}"
                                    )
                                    redis_client.lpop(keys.queue)
                                    continue

                                if task_info.get("provider") == "celery":
                                    # This is a celery job, pop it
                                    redis_client.lpop(keys.queue)
                                    job_id = task_info.get("jobId")
                                    room_id = task_info.get("room")  # Room context from job data

                                    # Get full job details
                                    job = JobManager.get_job(
                                        redis_client, job_id
                                    )
                                    if job:
                                        # Mark job as started by this worker
                                        JobManager.start_job(
                                            redis_client, job_id, worker_id
                                        )

                                        # Emit queue update (global for Celery)
                                        emit_queue_update(
                                            redis_client,
                                            None,  # Global notification
                                            category,
                                            extension,
                                            socketio,
                                        )

                                        # Get updated job details
                                        job = JobManager.get_job(
                                            redis_client, job_id
                                        )

                                        # Rename 'id' to 'jobId'
                                        job["jobId"] = job.pop("id")

                                        log.info(
                                            f"Assigned celery job {job_id} to worker {worker_id} from room {room_id}"
                                        )
                                        return job, 200
                continue  # Skip to next category for Celery workers
            # First, check global extensions for this worker
            global_user_extensions_key = ExtensionKeys.global_user_extensions_key(
                category, worker_id
            )
            global_registered_extensions = redis_client.smembers(
                global_user_extensions_key
            )

            # Check each globally registered extension for queued jobs
            for extension in global_registered_extensions:
                keys = ExtensionKeys.for_global_extension(category, extension)
                queue_length = redis_client.llen(keys.queue)

                if queue_length > 0:
                    # Get the first job from the queue
                    task_data = redis_client.lpop(keys.queue)
                    if task_data:
                        try:
                            task_info = json.loads(task_data)
                        except json.JSONDecodeError as e:
                            log.error(
                                f"Invalid JSON in global queue {keys.queue}: {e}, data: {task_data}"
                            )
                            continue

                        job_id = task_info.get("jobId")

                        # Get full job details
                        job = JobManager.get_job(redis_client, job_id)
                        if job:
                            # Move worker from idle to progressing (global keys)
                            redis_client.smove(
                                keys.idle_workers,
                                keys.progressing_workers,
                                worker_id,
                            )

                            # Mark job as started by this worker
                            JobManager.start_job(redis_client, job_id, worker_id)

                            # Emit queue update for global extension (broadcast to all)
                            # Note: For global extensions, emit without room context
                            emit_queue_update(
                                redis_client,
                                None,  # No specific room for global extensions
                                category,
                                extension,
                                socketio,
                            )

                            # Get updated job details (now with running status)
                            job = JobManager.get_job(redis_client, job_id)

                            # Rename 'id' to 'jobId' for consistency with client expectations
                            job["jobId"] = job.pop("id")

                            log.info(
                                f"Assigned global job {job_id} (extension={extension}, room={job.get('room')}) to worker {worker_id}"
                            )
                            return job, 200

            # Second, check room-specific extensions
            # We need to check all rooms where this worker is registered
            # Pattern: room:*:extensions:{category}:{worker_id}
            room_pattern = f"room:*:extensions:{category}:{worker_id}"
            for user_extensions_key in redis_client.scan_iter(match=room_pattern):
                # Extract room_id from key: room:{room_id}:extensions:{category}:{worker_id}
                parts = user_extensions_key.split(":")
                if len(parts) >= 5:
                    room_id = parts[1]

                    registered_extensions = redis_client.smembers(user_extensions_key)

                    # Check each registered extension for queued jobs
                    for extension in registered_extensions:
                        keys = ExtensionKeys.for_extension(room_id, category, extension)
                        queue_length = redis_client.llen(keys.queue)

                        if queue_length > 0:
                            # Get the first job from the queue
                            task_data = redis_client.lpop(keys.queue)
                            if task_data:
                                try:
                                    task_info = json.loads(task_data)
                                except json.JSONDecodeError as e:
                                    log.error(
                                        f"Invalid JSON in queue {keys.queue}: {e}, data: {task_data}"
                                    )
                                    continue

                                job_id = task_info.get("jobId")

                                # Get full job details
                                job = JobManager.get_job(redis_client, job_id)
                                if job:
                                    # Move worker from idle to progressing
                                    redis_client.smove(
                                        keys.idle_workers,
                                        keys.progressing_workers,
                                        worker_id,
                                    )

                                    # Mark job as started by this worker
                                    JobManager.start_job(redis_client, job_id, worker_id)

                                    # Emit queue update
                                    emit_queue_update(
                                        redis_client,
                                        room_id,
                                        category,
                                        extension,
                                        socketio,
                                    )

                                    # Get updated job details (now with running status)
                                    job = JobManager.get_job(redis_client, job_id)

                                    # Rename 'id' to 'jobId' for consistency
                                    job["jobId"] = job.pop("id")

                                    log.info(
                                        f"Assigned room job {job_id} (room={room_id}, extension={extension}) to worker {worker_id}"
                                    )
                                    return job, 200

        # No jobs available
        return {"error": "No jobs available"}, 400

    except Exception as e:
        log.error(
            f"Unexpected error in get_next_job_agnostic for worker {worker_id if 'worker_id' in locals() else 'unknown'}: {e}",
            exc_info=True,
        )
        return {"error": f"Internal server error: {str(e)}"}, 500
