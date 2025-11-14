"""Job management routes.

Handles job listing, status updates, and worker job polling.
"""

import json
import logging

from flask import Blueprint, current_app, jsonify, request

from zndraw.server import socketio
from zndraw.auth import require_auth

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
    """Release worker capacity and dispatch next task.

    Args:
        redis_client: Redis client
        socketio_instance: SocketIO instance
        worker_id: Worker session ID
        job: Job data dict containing category/extension
        room_id: Room identifier
        success: True if job completed, False if failed (affects log message)
    """
    from .job_dispatcher import assign_pending_jobs_for_extension, release_worker_capacity

    category = job.get("category")
    extension = job.get("extension")
    job_id = job.get("id")

    log.info(
        f"_transition_worker_to_idle called: worker_id={worker_id}, category={category}, extension={extension}, job_id={job_id}"
    )

    if not category or not extension or not job_id:
        log.warning(
            f"Missing required field in job data: category={category}, extension={extension}, job_id={job_id}"
        )
        return

    # Release worker capacity for this job
    release_worker_capacity(redis_client, worker_id, job_id)

    status_msg = (
        "finished and capacity restored" if success else "capacity restored after failure"
    )
    log.info(f"Worker {worker_id} {status_msg}")

    # Get all extensions this worker is registered for in this category
    # Try both room-scoped and global extensions
    registered_extensions_room = []
    registered_extensions_global = []

    user_extensions_key_room = ExtensionKeys.user_extensions_key(room_id, category, worker_id)
    registered_extensions_room = list(redis_client.smembers(user_extensions_key_room))

    user_extensions_key_global = ExtensionKeys.global_user_extensions_key(category, worker_id)
    registered_extensions_global = list(redis_client.smembers(user_extensions_key_global))

    log.info(
        f"Worker {worker_id} registered extensions - room: {registered_extensions_room}, global: {registered_extensions_global}"
    )

    # Try to assign pending jobs from ALL extensions this worker is registered for
    total_assigned = 0

    # Check room-scoped extensions
    for ext_name in registered_extensions_room:
        assigned = assign_pending_jobs_for_extension(
            redis_client, socketio_instance, room_id, category, ext_name, worker_id
        )
        if assigned > 0:
            log.info(f"Assigned {assigned} pending job(s) from room extension {ext_name} to worker {worker_id}")
            total_assigned += assigned
            break  # Worker now busy, stop trying to assign more

    # If no job assigned from room extensions, try global extensions
    if total_assigned == 0:
        for ext_name in registered_extensions_global:
            assigned = assign_pending_jobs_for_extension(
                redis_client, socketio_instance, None, category, ext_name, worker_id
            )
            if assigned > 0:
                log.info(f"Assigned {assigned} pending job(s) from global extension {ext_name} to worker {worker_id}")
                total_assigned += assigned
                break  # Worker now busy, stop trying to assign more

    if total_assigned == 0:
        log.debug(f"No pending jobs to assign to worker {worker_id}")


@jobs.route("/api/rooms/<string:room_id>/jobs", methods=["GET"])
def list_jobs(room_id: str):
    """List active jobs for a room."""
    redis_client = current_app.extensions["redis"]
    jobs = JobManager.list_all_jobs(redis_client, room_id)
    return jobs, 200


@jobs.route("/api/rooms/<string:room_id>/jobs/<string:job_id>", methods=["GET"])
def get_job(room_id: str, job_id: str):
    """Get details for a specific job."""
    redis_client = current_app.extensions["redis"]
    job = JobManager.get_job(redis_client, job_id)
    if not job:
        return {"error": "Job not found"}, 404
    return job, 200


@jobs.route("/api/rooms/<string:room_id>/jobs/<string:job_id>", methods=["DELETE"])
@require_auth
def delete_job(room_id: str, job_id: str):
    """Delete a specific job."""
    redis_client = current_app.extensions["redis"]
    job = JobManager.get_job(redis_client, job_id)
    if not job:
        return {"error": "Job not found"}, 404
    if job.get("status") in [JobStatus.ASSIGNED, JobStatus.PROCESSING]:
        return {"error": "Cannot delete a job that is assigned or processing"}, 400

    JobManager.delete_job(redis_client, job_id)

    return {"status": "success"}, 200


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
@require_auth
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

        # Update job status to processing (automatically emits job:state_changed)
        success = JobManager.start_processing(redis_client, job_id, worker_id, socketio=socketio)
        if not success:
            return {"error": "Failed to start job processing"}, 400

        log.info(f"Job {job_id} transitioned to processing by worker {worker_id}")

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

        # Mark job as completed (automatically emits job:state_changed)
        success = JobManager.complete_job(redis_client, job_id, socketio=socketio)
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
        # Mark job as failed (automatically emits job:state_changed with error metadata)
        success = JobManager.fail_job(redis_client, job_id, error, socketio=socketio)
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


