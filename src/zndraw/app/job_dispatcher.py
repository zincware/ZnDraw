"""Job dispatcher for assigning pending jobs to idle workers.

Handles the logic of matching pending jobs with available workers
and emitting job assignments via Socket.IO.
"""

import logging
from datetime import datetime

from flask_socketio import emit

from .constants import SocketEvents
from .job_manager import JobManager, JobStatus
from .redis_keys import ExtensionKeys

log = logging.getLogger(__name__)


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

    # Get idle workers
    if worker_id is not None:
        # Check if specific worker is idle
        is_idle = redis_client.sismember(keys.idle_workers, worker_id)
        idle_workers = [worker_id] if is_idle else []
    else:
        # Get all idle workers
        idle_workers = list(redis_client.smembers(keys.idle_workers))

    if not idle_workers:
        log.debug(
            f"No idle workers for {category}/{extension} (room={room_id})"
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

    # Assign pending jobs to idle workers
    for idx, job_id in enumerate(pending_jobs):
        if idx >= len(idle_workers):
            # No more idle workers available
            break

        current_worker_id = idle_workers[idx]

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

        # Move worker from idle to assigned state
        moved = redis_client.smove(
            keys.idle_workers,
            keys.progressing_workers,
            current_worker_id,
        )

        if not moved:
            # Worker was removed (disconnected), re-queue the job
            timestamp = datetime.utcnow().timestamp()
            redis_client.zadd(keys.pending_jobs, {job_id: timestamp})
            log.warning(
                f"Worker {current_worker_id} no longer idle, re-queued job {job_id}"
            )
            continue

        # Update job status to ASSIGNED
        success = JobManager.assign_job(redis_client, job_id, current_worker_id)
        if not success:
            # Job not in PENDING state, roll back worker state
            redis_client.smove(
                keys.progressing_workers,
                keys.idle_workers,
                current_worker_id,
            )
            log.warning(
                f"Failed to assign job {job_id} to worker {current_worker_id}, "
                f"job not in PENDING state"
            )
            continue

        # Store worker state
        worker_state_key = ExtensionKeys.worker_state_key(
            room_id, category, extension, current_worker_id
        )
        worker_state = {
            "state": "assigned",
            "assigned_at": datetime.utcnow().isoformat(),
            "job_id": job_id,
        }
        redis_client.hset(worker_state_key, mapping=worker_state)

        # Emit job:assigned to specific worker (just the ID, worker will fetch details)
        _emit_job_assigned(socketio, current_worker_id, job_id)

        # Emit job:state_changed to room
        _emit_job_state_changed(socketio, job_data["room"], job_id, JobStatus.ASSIGNED)

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
        emit(
            SocketEvents.JOB_ASSIGNED,
            payload,
            to=worker_id,
            namespace="/",
        )
        log.debug(f"Emitted job:assigned to worker {worker_id} for job {job_id}")
    except Exception as e:
        log.error(f"Failed to emit job:assigned to worker {worker_id}: {e}")


def _emit_job_state_changed(
    socketio,
    room_id: str,
    job_id: str,
    status: str,
    metadata: dict | None = None,
) -> None:
    """Emit job:state_changed event to all clients in a room.

    Parameters
    ----------
    socketio
        SocketIO instance
    room_id : str
        Room ID
    job_id : str
        Job ID
    status : str
        New job status
    metadata : dict | None
        Optional metadata (error message, result, etc.)
    """
    payload = {
        "jobId": job_id,
        "status": status,
        "room": room_id,
    }

    if metadata:
        payload["metadata"] = metadata

    try:
        emit(
            SocketEvents.JOB_STATE_CHANGED,
            payload,
            to=room_id,
            namespace="/",
        )
        log.debug(f"Emitted job:state_changed to room {room_id} for job {job_id}: {status}")
    except Exception as e:
        log.error(f"Failed to emit job:state_changed to room {room_id}: {e}")
