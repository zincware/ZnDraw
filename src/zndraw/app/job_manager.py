"""Job tracking and management for ZnDraw extensions."""

import json
import logging
import uuid
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Any, Optional

from dateutil.parser import isoparse

from .redis_keys import JobKeys, RoomKeys, WorkerKeys

log = logging.getLogger(__name__)


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
        SocketIO instance (can be None if not available)
    room_id : str
        Room ID
    job_id : str
        Job ID
    status : str
        New job status
    metadata : dict | None
        Optional metadata (error message, result, etc.)
    """
    if not socketio:
        return  # Silently skip if socketio not available

    from .constants import SocketEvents

    payload = {
        "jobId": job_id,
        "status": status,
        "room": room_id,
    }

    if metadata:
        payload["metadata"] = metadata

    try:
        socketio.emit(
            SocketEvents.JOB_STATE_CHANGED,
            payload,
            to=f"room:{room_id}",
            namespace="/",
        )
        log.debug(f"Emitted job:state_changed to room {room_id} for job {job_id}: {status}")
    except Exception as e:
        log.error(f"Failed to emit job:state_changed to room {room_id}: {e}")


@dataclass
class Job:
    id: str
    room: str
    category: str
    extension: str
    data: dict
    user_name: str
    status: str
    provider: str  # "celery" or worker count
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    worker_id: Optional[str] = None  # SID or "celery:{task_id}"
    error: Optional[str] = None
    wait_time_ms: Optional[int] = None  # started_at - created_at
    execution_time_ms: Optional[int] = None  # completed_at - started_at


class JobStatus(str, Enum):
    """Job status states."""

    PENDING = "pending"  # Waiting for idle worker
    ASSIGNED = "assigned"  # Emitted to worker, awaiting confirmation
    PROCESSING = "processing"  # Worker confirmed and actively processing
    COMPLETED = "completed"
    FAILED = "failed"


def _calculate_durations(
    created_at: str | None,
    started_at: str | None,
    completed_at: str | None,
) -> tuple[int | None, int | None]:
    """Calculate wait and execution times in milliseconds.

    Args:
        created_at: ISO 8601 timestamp when job was created
        started_at: ISO 8601 timestamp when job started
        completed_at: ISO 8601 timestamp when job completed

    Returns:
        Tuple of (wait_time_ms, execution_time_ms)
    """
    wait_time_ms = None
    execution_time_ms = None

    try:
        if started_at and created_at:
            wait_time_ms = int(
                (isoparse(started_at) - isoparse(created_at)).total_seconds() * 1000
            )

        if completed_at and started_at:
            execution_time_ms = int(
                (isoparse(completed_at) - isoparse(started_at)).total_seconds() * 1000
            )
    except Exception:
        pass  # Invalid timestamps, leave as None

    return wait_time_ms, execution_time_ms


class JobManager:
    """Manages job lifecycle and tracking in Redis."""

    @staticmethod
    def create_job(
        redis_client: Any,
        room: str,
        category: str,
        extension: str,
        data: dict,
        user_name: str,
        provider: str,
        public: bool = False,
        ttl: int = 86400,
        initial_status: str = JobStatus.PENDING,
        socketio=None,
    ) -> str:
        """Create a new job and store in Redis.

        Args:
            redis_client: Redis client instance
            room: Room identifier
            category: Extension category
            extension: Extension name
            data: Job input parameters
            user_name: Username who submitted the job
            provider: "celery" or worker count
            public: Whether this is a public/global extension (default: False)
            ttl: Time-to-live in seconds (default: 86400 = 24 hours)
            initial_status: Initial job status (default: PENDING)
            socketio: SocketIO instance for emitting events (optional)

        Returns:
            job_id: Unique job identifier
        """
        job_id = str(uuid.uuid4())
        now = datetime.utcnow().isoformat()

        job_data = {
            "id": job_id,
            "room": room,
            "category": category,
            "extension": extension,
            "data": json.dumps(data),
            "user_name": user_name,
            "status": initial_status,
            "provider": provider,
            "public": str(public).lower(),  # Store as "true" or "false" string for Redis
            "created_at": now,
            "started_at": "",
            "completed_at": "",
            "worker_id": "",
            "error": "",
            "wait_time_ms": "",
            "execution_time_ms": "",
        }

        # Store job data
        job_keys = JobKeys(job_id)
        room_keys = RoomKeys(room)
        redis_client.hset(job_keys.hash_key(), mapping=job_data)

        # Set TTL
        redis_client.expire(job_keys.hash_key(), ttl)

        # Add to indexes
        redis_client.sadd(room_keys.jobs_active(), job_id)
        redis_client.zadd(
            room_keys.jobs_by_time(), {job_id: datetime.utcnow().timestamp()}
        )
        redis_client.sadd(room_keys.extension_jobs(category, extension), job_id)

        # Update Prometheus metrics - increment counter for submitted jobs
        from zndraw.analytics import (
            running_analysis,
            running_modifiers,
            running_selections,
        )

        gauge_map = {
            'modifiers': running_modifiers,
            'selections': running_selections,
            'analysis': running_analysis,
        }
        gauge = gauge_map.get(category)
        if gauge:
            # Use 'room' as scope since jobs are room-scoped
            gauge.labels(scope='room').inc()

        # Emit job state change event
        _emit_job_state_changed(socketio, room, job_id, initial_status)

        return job_id

    @staticmethod
    def assign_job(redis_client: Any, job_id: str, worker_id: str, socketio=None) -> bool:
        """Mark job as assigned to a worker.

        Transition: PENDING -> ASSIGNED

        Args:
            redis_client: Redis client instance
            job_id: Job identifier
            worker_id: Worker identifier (SID or "celery:{task_id}")
            socketio: SocketIO instance for emitting events (optional)

        Returns:
            True if successful, False if job not in PENDING state
        """
        job_keys = JobKeys(job_id)

        # Check if job exists and is pending
        job_data = redis_client.hgetall(job_keys.hash_key())
        if not job_data or job_data.get("status") != JobStatus.PENDING:
            return False

        # Update status and add to worker's job set
        update_data = {
            "status": JobStatus.ASSIGNED,
            "worker_id": worker_id,
        }

        redis_client.hset(job_keys.hash_key(), mapping=update_data)

        # Remove job from extension's pending_jobs sorted set
        from .redis_keys import ExtensionKeys
        category = job_data.get("category")
        extension = job_data.get("extension")
        room = job_data.get("room")
        is_public = job_data.get("public") == "true"

        if category and extension:
            ext_room = None if is_public else room
            if ext_room:
                keys = ExtensionKeys.for_extension(ext_room, category, extension)
            else:
                keys = ExtensionKeys.for_global_extension(category, extension)
            redis_client.zrem(keys.pending_jobs, job_id)

        # Add job to worker's active jobs set (for disconnect cleanup)
        worker_keys = WorkerKeys(worker_id)
        redis_client.sadd(worker_keys.active_jobs(), job_id)

        # Emit job state change event
        _emit_job_state_changed(socketio, room, job_id, JobStatus.ASSIGNED)

        return True

    @staticmethod
    def start_processing(redis_client: Any, job_id: str, worker_id: str, socketio=None) -> bool:
        """Mark job as processing (worker confirmed and started work).

        Transition: ASSIGNED -> PROCESSING

        Args:
            redis_client: Redis client instance
            job_id: Job identifier
            worker_id: Worker identifier (SID or "celery:{task_id}")
            socketio: SocketIO instance for emitting events (optional)

        Returns:
            True if successful, False if job not in ASSIGNED state
        """
        job_keys = JobKeys(job_id)

        # Check if job exists and is assigned
        status = redis_client.hget(job_keys.hash_key(), "status")
        if not status or status != JobStatus.ASSIGNED:
            return False

        # Verify worker_id matches
        assigned_worker = redis_client.hget(job_keys.hash_key(), "worker_id")
        if assigned_worker != worker_id:
            return False

        # Update status
        now = datetime.utcnow().isoformat()
        created_at = redis_client.hget(job_keys.hash_key(), "created_at")
        room = redis_client.hget(job_keys.hash_key(), "room")

        # Calculate wait time (from creation to processing start)
        wait_time_ms, _ = _calculate_durations(created_at, now, None)

        update_data = {
            "status": JobStatus.PROCESSING,
            "started_at": now,
        }
        if wait_time_ms is not None:
            update_data["wait_time_ms"] = str(wait_time_ms)

        redis_client.hset(job_keys.hash_key(), mapping=update_data)

        # Emit job state change event
        if room:
            _emit_job_state_changed(socketio, room, job_id, JobStatus.PROCESSING)

        return True

    @staticmethod
    def complete_job(redis_client: Any, job_id: str, socketio=None) -> bool:
        """Mark job as completed and calculate execution time.

        Args:
            redis_client: Redis client instance
            job_id: Job identifier
            socketio: SocketIO instance for emitting events (optional)

        Returns:
            True if successful, False if job not found
        """
        job_keys = JobKeys(job_id)

        if not redis_client.exists(job_keys.hash_key()):
            return False

        now = datetime.utcnow().isoformat()
        started_at = redis_client.hget(job_keys.hash_key(), "started_at")

        # Calculate execution time
        _, execution_time_ms = _calculate_durations(None, started_at, now)

        update_data = {
            "status": JobStatus.COMPLETED,
            "completed_at": now,
        }
        if execution_time_ms is not None:
            update_data["execution_time_ms"] = str(execution_time_ms)

        redis_client.hset(job_keys.hash_key(), mapping=update_data)

        # Move from active to inactive set
        job_data = redis_client.hgetall(job_keys.hash_key())
        room = job_data.get("room")
        if room:
            room_keys = RoomKeys(room)
            redis_client.smove(
                room_keys.jobs_active(), room_keys.jobs_inactive(), job_id
            )

        # Remove job from worker's active jobs set
        worker_id = job_data.get("worker_id")
        if worker_id:
            worker_keys = WorkerKeys(worker_id)
            redis_client.srem(worker_keys.active_jobs(), job_id)

        # Emit job state change event
        if room:
            _emit_job_state_changed(socketio, room, job_id, JobStatus.COMPLETED)

        return True

    @staticmethod
    def fail_job(redis_client: Any, job_id: str, error: str, socketio=None) -> bool:
        """Mark job as failed and calculate execution time.

        Args:
            redis_client: Redis client instance
            job_id: Job identifier
            error: Error message
            socketio: SocketIO instance for emitting events (optional)

        Returns:
            True if successful, False if job not found
        """
        job_keys = JobKeys(job_id)

        if not redis_client.exists(job_keys.hash_key()):
            return False

        now = datetime.utcnow().isoformat()
        started_at = redis_client.hget(job_keys.hash_key(), "started_at")

        # Calculate execution time
        _, execution_time_ms = _calculate_durations(None, started_at, now)

        update_data = {
            "status": JobStatus.FAILED,
            "completed_at": now,
            "error": error,
        }
        if execution_time_ms is not None:
            update_data["execution_time_ms"] = str(execution_time_ms)

        redis_client.hset(job_keys.hash_key(), mapping=update_data)

        # Move from active to inactive set
        job_data = redis_client.hgetall(job_keys.hash_key())
        room = job_data.get("room")
        if room:
            room_keys = RoomKeys(room)
            redis_client.smove(
                room_keys.jobs_active(), room_keys.jobs_inactive(), job_id
            )

        # Remove job from worker's active jobs set
        worker_id = job_data.get("worker_id")
        if worker_id:
            worker_keys = WorkerKeys(worker_id)
            redis_client.srem(worker_keys.active_jobs(), job_id)

        # Emit job state change event with error metadata
        if room:
            _emit_job_state_changed(socketio, room, job_id, JobStatus.FAILED, metadata={"error": error})

        return True

    @staticmethod
    def get_job(redis_client: Any, job_id: str) -> Optional[dict]:
        """Get job details.

        Args:
            redis_client: Redis client instance
            job_id: Job identifier

        Returns:
            Job data dict or None if not found
        """
        job_keys = JobKeys(job_id)
        job_data = redis_client.hgetall(job_keys.hash_key())

        if not job_data:
            return None

        # Parse JSON fields
        if job_data.get("data"):
            job_data["data"] = json.loads(job_data["data"])

        return job_data

    @staticmethod
    def list_active_jobs(redis_client: Any, room: str) -> list[dict]:
        """List all active (queued or running) jobs for a room.

        Args:
            redis_client: Redis client instance
            room: Room identifier

        Returns:
            List of job data dicts
        """
        room_keys = RoomKeys(room)
        job_ids = redis_client.smembers(room_keys.jobs_active())

        jobs = []
        for job_id in job_ids:
            job_data = JobManager.get_job(redis_client, job_id)
            if job_data:
                jobs.append(job_data)

        return jobs

    @staticmethod
    def list_extension_jobs(
        redis_client: Any, room: str, category: str, extension: str, limit: int = 50
    ) -> list[dict]:
        """List recent jobs for a specific extension.

        Args:
            redis_client: Redis client instance
            room: Room identifier
            category: Extension category
            extension: Extension name
            limit: Maximum number of jobs to return

        Returns:
            List of job data dicts, newest first
        """
        room_keys = RoomKeys(room)
        job_ids = redis_client.smembers(room_keys.extension_jobs(category, extension))

        jobs = []
        for job_id in job_ids:
            job_data = JobManager.get_job(redis_client, job_id)
            if job_data:
                jobs.append(job_data)

        # Sort by created_at, newest first
        jobs.sort(key=lambda x: x.get("created_at", ""), reverse=True)

        return jobs[:limit]

    @staticmethod
    def list_inactive_jobs(redis_client: Any, room: str) -> list[dict]:
        """List all inactive (completed or failed) jobs for a room.
        Args:
            redis_client: Redis client instance
            room: Room identifier
        Returns:
            List of job data dicts
        """
        room_keys = RoomKeys(room)
        job_ids = redis_client.smembers(room_keys.jobs_inactive())

        jobs = []
        for job_id in job_ids:
            job_data = JobManager.get_job(redis_client, job_id)
            if job_data:
                jobs.append(job_data)

        return jobs

    @staticmethod
    def list_all_jobs(redis_client: Any, room: str) -> list[dict]:
        """List all jobs for a room.
        Args:
            redis_client: Redis client instance
            room: Room identifier
        Returns:
            List of job data dicts
        """
        active_jobs = JobManager.list_active_jobs(redis_client, room)
        inactive_jobs = JobManager.list_inactive_jobs(redis_client, room)
        return active_jobs + inactive_jobs

    @staticmethod
    def delete_job(redis_client: Any, job_id: str) -> bool:
        """Delete a job from Redis.

        Args:
            redis_client: Redis client instance
            job_id: Job identifier

        Returns:
            True if successful, False if job not found
        """
        job_keys = JobKeys(job_id)
        job_data = redis_client.hgetall(job_keys.hash_key())

        if not job_data:
            return False

        room = job_data.get("room")
        category = job_data.get("category")
        extension = job_data.get("extension")

        pipe = redis_client.pipeline()
        pipe.delete(job_keys.hash_key())

        if room:
            room_keys = RoomKeys(room)
            pipe.srem(room_keys.jobs_active(), job_id)
            pipe.srem(room_keys.jobs_inactive(), job_id)
            pipe.zrem(room_keys.jobs_by_time(), job_id)
            if category and extension:
                pipe.srem(room_keys.extension_jobs(category, extension), job_id)

        pipe.execute()
        return True

    @staticmethod
    def get_extension_stats(
        redis_client: Any, room: str, category: str, extension: str
    ) -> dict:
        """Get job statistics for an extension.

        Args:
            redis_client: Redis client instance
            room: Room identifier
            category: Extension category
            extension: Extension name

        Returns:
            Dict with pending, assigned, processing, completed, failed counts
        """
        jobs = JobManager.list_extension_jobs(redis_client, room, category, extension)

        stats = {
            "pending": 0,
            "assigned": 0,
            "processing": 0,
            "completed": 0,
            "failed": 0,
        }

        for job in jobs:
            status = job.get("status")
            if status in stats:
                stats[status] += 1

        return stats
