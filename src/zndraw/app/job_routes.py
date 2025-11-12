"""Job management routes.

Handles job listing, status updates, and worker job polling.
"""

import json
import logging

from flask import Blueprint, current_app, jsonify, request

from zndraw.server import socketio

from .job_manager import JobManager, JobStatus
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

        # Assign pending jobs for this extension to the now-idle worker
        from .job_dispatcher import assign_pending_jobs_for_extension

        # Determine if this is a public (global) or room-scoped extension
        public = job.get("public") == "true"
        ext_room_id = None if public else room_id

        assigned_count = assign_pending_jobs_for_extension(
            redis_client, socketio_instance, ext_room_id, category, extension, worker_id
        )

        if assigned_count > 0:
            log.info(f"Assigned {assigned_count} pending job(s) to worker {worker_id}")
        else:
            log.debug(f"No pending jobs to assign to worker {worker_id}")
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
        if job.get("status") in [JobStatus.ASSIGNED, JobStatus.PROCESSING]:
            return {"error": "Cannot delete a job that is assigned or processing"}, 400

        category = job.get("category")
        extension = job.get("extension")

        JobManager.delete_job(redis_client, job_id)

        return {"status": "success"}, 200

    job = JobManager.get_job(redis_client, job_id)
    if not job:
        return {"error": "Job not found"}, 404
    return job, 200


def _update_worker_state_to_running(
    redis_client, worker_id: str, job: dict, room_id: str
) -> None:
    """Update worker state from 'assigned' to 'running'.

    Parameters
    ----------
    redis_client
        Redis client instance
    worker_id : str
        Worker ID
    job : dict
        Job data
    room_id : str
        Room ID
    """
    from datetime import datetime

    category = job.get("category")
    extension = job.get("extension")
    public = job.get("public") == "true"
    job_id = job.get("id")

    # Determine if this is a global or room-scoped extension
    ext_room_id = None if public else room_id

    # Update worker state hash
    worker_state_key = ExtensionKeys.worker_state_key(
        ext_room_id, category, extension, worker_id
    )

    worker_state = {
        "state": "running",
        "started_at": datetime.utcnow().isoformat(),
        "job_id": job_id,
    }

    redis_client.hset(worker_state_key, mapping=worker_state)
    log.debug(f"Updated worker {worker_id} state to running for job {job_id}")


@jobs.route("/api/jobs/<string:job_id>", methods=["GET"])
def get_job_details(job_id: str):
    """Get full job details for a worker to execute.

    This endpoint is called by workers after receiving a job:assigned socket event.
    Returns all information needed to execute the job.

    Returns
    -------
    dict
        Job details including:
        - jobId: Job identifier
        - room: Target room ID
        - category: Extension category
        - extension: Extension name
        - data: Job parameters (parsed JSON)
        - public: Whether this is a global extension
        - status: Current job status
        - created_at: ISO timestamp
    """
    redis_client = current_app.extensions["redis"]

    # Get job from Redis
    job = JobManager.get_job(redis_client, job_id)

    if not job:
        log.error(f"Job {job_id} not found")
        return {"error": "Job not found"}, 404

    # Return job details in camelCase for consistency with other API responses
    response = {
        "jobId": job["id"],
        "room": job["room"],
        "category": job["category"],
        "extension": job["extension"],
        "data": job["data"],  # Already parsed from JSON by JobManager.get_job()
        "public": job.get("public") == "true",
        "status": job["status"],
        "createdAt": job.get("created_at"),
    }

    log.info(f"Worker fetching job {job_id}: {job['category']}/{job['extension']}")

    return jsonify(response), 200


@jobs.route("/api/rooms/<string:room_id>/jobs/<string:job_id>/status", methods=["PUT"])
def update_job_status(room_id: str, job_id: str):
    """Update a job's status (processing, completed, or failed).

    Called by workers during job lifecycle:
    1. Worker receives job:assigned event with jobId
    2. Worker fetches job via GET /api/jobs/{jobId}
    3. Worker calls PUT with status='processing' (ASSIGNED → PROCESSING)
    4. Worker executes job
    5. Worker calls PUT with status='completed' or 'failed' (PROCESSING → COMPLETED/FAILED)

    Request body:
        {
            "status": "processing" | "completed" | "failed",
            "error": <error message>,  // for failed status
            "workerId": <worker_id>
        }
    """
    data = request.get_json() or {}
    status = data.get("status")
    worker_id = data.get("workerId")

    if status not in ["processing", "completed", "failed"]:
        return {"error": "Status must be 'processing', 'completed', or 'failed'"}, 400

    if not worker_id:
        return {"error": "workerId is required"}, 400

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

    # Handle status transitions
    if status == "processing":
        # Transition: ASSIGNED → PROCESSING
        # Validate job is in assigned state
        if job.get("status") != JobStatus.ASSIGNED:
            log.error(
                f"Job {job_id} cannot transition to processing (current status: {job.get('status')}, expected: {JobStatus.ASSIGNED})"
            )
            return {
                "error": f"Job must be in 'assigned' state to start processing (current: {job.get('status')})"
            }, 400

        # Validate worker ID matches
        if job.get("worker_id") != worker_id:
            log.error(
                f"Worker ID mismatch: job assigned to {job.get('worker_id')}, but {worker_id} tried to start it"
            )
            return {"error": "Worker ID does not match job's assigned worker"}, 400

        # Update job status to processing
        success = JobManager.start_processing(redis_client, job_id, worker_id)
        if not success:
            return {"error": "Failed to start job processing"}, 400

        log.info(f"Job {job_id} transitioned to processing by worker {worker_id}")

        # Update worker state to running
        _update_worker_state_to_running(redis_client, worker_id, job, room_id)

        # Emit job state change
        from .job_dispatcher import _emit_job_state_changed
        _emit_job_state_changed(socketio, room_id, job_id, JobStatus.PROCESSING)

        return {"status": "success"}, 200

    elif status == "completed":
        # Transition: PROCESSING → COMPLETED
        # Validate job is in processing state
        if job.get("status") != JobStatus.PROCESSING:
            log.error(
                f"Job {job_id} cannot be completed (current status: {job.get('status')}, expected: {JobStatus.PROCESSING})"
            )
            return {
                "error": f"Job must be in 'processing' state to complete (current: {job.get('status')})"
            }, 400

        # Validate worker ID matches
        if job.get("worker_id") != worker_id:
            log.error(
                f"Worker ID mismatch: job assigned to {job.get('worker_id')}, but {worker_id} tried to complete it"
            )
            return {"error": "Worker ID does not match job's worker"}, 400

        success = JobManager.complete_job(redis_client, job_id)
        if not success:
            log.error(f"Failed to mark job {job_id} as completed")
            return {"error": "Failed to complete job"}, 400
        log.info(f"Job {job_id} completed in room {room_id}")
        worker_success = True

    else:  # status == "failed"
        # Transition: PROCESSING → FAILED
        # Validate job is in processing state
        if job.get("status") != JobStatus.PROCESSING:
            log.error(
                f"Job {job_id} cannot be marked as failed (current status: {job.get('status')}, expected: {JobStatus.PROCESSING})"
            )
            return {
                "error": f"Job must be in 'processing' state to fail (current: {job.get('status')})"
            }, 400

        # Validate worker ID matches
        if job.get("worker_id") != worker_id:
            log.error(
                f"Worker ID mismatch: job assigned to {job.get('worker_id')}, but {worker_id} tried to mark it failed"
            )
            return {"error": "Worker ID does not match job's worker"}, 400

        error = data.get("error", "Unknown error")
        success = JobManager.fail_job(redis_client, job_id, error)
        if not success:
            return {"error": "Failed to mark job as failed"}, 400
        log.error(f"Job {job_id} failed in room {room_id}: {error}")

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


