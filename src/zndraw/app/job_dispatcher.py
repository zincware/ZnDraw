"""Job dispatcher for assigning pending jobs to idle workers.

Handles the logic of matching pending jobs with available workers
and emitting job assignments via Socket.IO.
"""

import logging
from datetime import datetime

from .constants import SocketEvents
from .job_manager import JobManager, JobStatus
from .redis_keys import ExtensionKeys

log = logging.getLogger(__name__)


def get_available_workers(
    redis_client, extension_keys: ExtensionKeys, min_capacity: int = 1
) -> list[str]:
    """Get workers with available capacity for this extension.

    Parameters
    ----------
    redis_client
        Redis client instance
    extension_keys : ExtensionKeys
        ExtensionKeys instance for this extension
    min_capacity : int
        Minimum capacity required (default: 1)

    Returns
    -------
    list[str]
        Worker IDs with available capacity
    """
    # Get all workers registered for this extension
    all_workers = redis_client.hkeys(extension_keys.workers)

    if not all_workers:
        return []

    # Build pipeline to check capacity for all workers
    pipe = redis_client.pipeline()
    for worker_id in all_workers:
        capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
        pipe.get(capacity_key)

    capacities = pipe.execute()

    # Filter workers with sufficient capacity
    available_workers = []
    for worker_id, capacity in zip(all_workers, capacities):
        capacity_val = int(capacity) if capacity else 0
        if capacity_val >= min_capacity:
            available_workers.append(worker_id)

    return available_workers


def assign_job_to_worker(redis_client, worker_id: str, job_id: str) -> bool:
    """Assign a job to a worker, updating capacity.

    Parameters
    ----------
    redis_client
        Redis client instance
    worker_id : str
        Worker ID
    job_id : str
        Job ID

    Returns
    -------
    bool
        True if assignment successful, False if worker has no capacity
    """
    capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
    active_jobs_key = ExtensionKeys.worker_active_jobs_key(worker_id)

    # Check capacity
    capacity = int(redis_client.get(capacity_key) or 0)
    if capacity < 1:
        return False  # No capacity

    # Decrement capacity and track active job
    redis_client.decr(capacity_key)
    redis_client.sadd(active_jobs_key, job_id)
    return True


def release_worker_capacity(redis_client, worker_id: str, job_id: str) -> None:
    """Release worker capacity when job completes.

    Parameters
    ----------
    redis_client
        Redis client instance
    worker_id : str
        Worker ID
    job_id : str
        Job ID
    """
    capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
    active_jobs_key = ExtensionKeys.worker_active_jobs_key(worker_id)

    pipe = redis_client.pipeline()
    pipe.incr(capacity_key)
    pipe.srem(active_jobs_key, job_id)
    pipe.execute()


def assign_pending_jobs_for_extension(
    redis_client,
    socketio,
    room_id: str | None,
    category: str,
    extension: str,
    worker_id: str | None = None,
) -> int:
    """Assign pending jobs to idle workers for a specific extension.

    Parameters
    ----------
    redis_client
        Redis client instance
    socketio
        SocketIO instance for emitting events
    room_id : str | None
        Room ID for room-scoped extensions, None for global extensions
    category : str
        Extension category (modifiers, selections, analysis, settings)
    extension : str
        Extension name
    worker_id : str | None
        Specific worker ID to assign to (if provided), otherwise assigns to any idle worker

    Returns
    -------
    int
        Number of jobs assigned
    """
    # Get extension keys
    if room_id is None:
        keys = ExtensionKeys.for_global_extension(category, extension)
    else:
        keys = ExtensionKeys.for_extension(room_id, category, extension)

    # Get workers with available capacity
    if worker_id is not None:
        # Check if specific worker has capacity
        capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
        capacity = int(redis_client.get(capacity_key) or 0)
        # Also check if worker is registered for this extension
        is_registered = redis_client.hexists(keys.workers, worker_id)
        available_workers = [worker_id] if (is_registered and capacity >= 1) else []
    else:
        # Get all workers with available capacity
        available_workers = get_available_workers(redis_client, keys)

    if not available_workers:
        log.debug(
            f"No available workers for {category}/{extension} (room={room_id})"
        )
        return 0

    # Get pending jobs (sorted by timestamp, oldest first)
    pending_jobs = redis_client.zrange(keys.pending_jobs, 0, -1, withscores=False)

    if not pending_jobs:
        log.debug(
            f"No pending jobs for {category}/{extension} (room={room_id})"
        )
        return 0

    assignments_made = 0

    # Assign pending jobs to available workers
    for idx, job_id in enumerate(pending_jobs):
        if idx >= len(available_workers):
            # No more available workers
            break

        current_worker_id = available_workers[idx]

        # Remove job from pending queue
        removed = redis_client.zrem(keys.pending_jobs, job_id)
        if not removed:
            # Job was already removed (race condition)
            continue

        # Get job data
        job_data = JobManager.get_job(redis_client, job_id)
        if not job_data:
            log.warning(f"Job {job_id} not found in Redis, skipping")
            continue

        # Assign job to worker (updates capacity and active jobs)
        assigned = assign_job_to_worker(redis_client, current_worker_id, job_id)

        if not assigned:
            # Worker has no capacity (may have been assigned another job), re-queue the job
            timestamp = datetime.utcnow().timestamp()
            redis_client.zadd(keys.pending_jobs, {job_id: timestamp})
            log.warning(
                f"Worker {current_worker_id} has no capacity, re-queued job {job_id}"
            )
            continue

        # Update job status to ASSIGNED (automatically emits job:state_changed)
        success = JobManager.assign_job(redis_client, job_id, current_worker_id, socketio=socketio)
        if not success:
            # Job not in PENDING state, roll back worker capacity
            release_worker_capacity(redis_client, current_worker_id, job_id)
            log.warning(
                f"Failed to assign job {job_id} to worker {current_worker_id}, "
                f"job not in PENDING state"
            )
            continue

        # Emit job:assigned to specific worker (just the ID, worker will fetch details)
        _emit_job_assigned(socketio, current_worker_id, job_id)

        assignments_made += 1
        log.info(
            f"Assigned job {job_id} ({category}/{extension}) to worker {current_worker_id}"
        )

    return assignments_made


def _emit_job_assigned(socketio, worker_id: str, job_id: str) -> None:
    """Emit job:assigned event to a specific worker.

    Only sends the job ID - worker must fetch full job details via GET /api/jobs/{job_id}.
    This follows REST principles: Socket for notification, REST API for data.

    Parameters
    ----------
    socketio
        SocketIO instance
    worker_id : str
        Worker socket ID
    job_id : str
        Job ID to assign
    """
    payload = {
        "jobId": job_id,
    }

    try:
        socketio.emit(
            SocketEvents.JOB_ASSIGNED,
            payload,
            to=worker_id,
            namespace="/",
        )
        log.debug(f"Emitted job:assigned to worker {worker_id} for job {job_id}")
    except Exception as e:
        log.error(f"Failed to emit job:assigned to worker {worker_id}: {e}")
