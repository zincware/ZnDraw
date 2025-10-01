"""Job tracking and management for ZnDraw extensions."""

import json
import uuid
from datetime import datetime
from typing import Any, Optional
from enum import Enum
from dataclasses import dataclass

@dataclass
class Job:
    id: str
    room: str
    category: str
    extension: str
    data: dict
    user_id: str
    status: str
    provider: str  # "celery" or worker count
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    worker_id: Optional[str] = None  # SID or "celery:{task_id}"
    error: Optional[str] = None
    result: Optional[dict] = None


class JobStatus(str, Enum):
    """Job status states."""
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class JobManager:
    """Manages job lifecycle and tracking in Redis."""

    @staticmethod
    def create_job(
        redis_client: Any,
        room: str,
        category: str,
        extension: str,
        data: dict,
        user_id: str,
        provider: str,
    ) -> str:
        """Create a new job and store in Redis.

        Args:
            redis_client: Redis client instance
            room: Room identifier
            category: Extension category
            extension: Extension name
            data: Job input parameters
            user_id: User who submitted the job
            provider: "celery" or worker count

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
            "user_id": user_id,
            "status": JobStatus.QUEUED,
            "provider": provider,
            "created_at": now,
            "started_at": "",
            "completed_at": "",
            "worker_id": "",
            "error": "",
            "result": "",
        }

        # Store job data
        job_key = f"job:{job_id}"
        redis_client.hset(job_key, mapping=job_data)

        # Set TTL to 24 hours
        redis_client.expire(job_key, 86400)

        # Add to indexes
        redis_client.sadd(f"room:{room}:jobs:active", job_id)
        redis_client.zadd(
            f"room:{room}:jobs:by_time", {job_id: datetime.utcnow().timestamp()}
        )
        redis_client.sadd(
            f"room:{room}:extension:{category}:{extension}:jobs", job_id
        )

        return job_id

    @staticmethod
    def start_job(redis_client: Any, job_id: str, worker_id: str) -> bool:
        """Mark job as running.

        Args:
            redis_client: Redis client instance
            job_id: Job identifier
            worker_id: Worker identifier (SID or "celery:{task_id}")

        Returns:
            True if successful, False if job not found or already started
        """
        job_key = f"job:{job_id}"

        # Check if job exists and is queued
        status = redis_client.hget(job_key, "status")
        if not status or status != JobStatus.QUEUED:
            return False

        # Update status
        now = datetime.utcnow().isoformat()
        redis_client.hset(
            job_key,
            mapping={
                "status": JobStatus.RUNNING,
                "started_at": now,
                "worker_id": worker_id,
            },
        )

        return True

    @staticmethod
    def complete_job(
        redis_client: Any, job_id: str, result: Optional[dict] = None
    ) -> bool:
        """Mark job as completed.

        Args:
            redis_client: Redis client instance
            job_id: Job identifier
            result: Optional job result data

        Returns:
            True if successful, False if job not found
        """
        job_key = f"job:{job_id}"

        if not redis_client.exists(job_key):
            return False

        now = datetime.utcnow().isoformat()
        redis_client.hset(
            job_key,
            mapping={
                "status": JobStatus.COMPLETED,
                "completed_at": now,
                "result": json.dumps(result) if result else "",
            },
        )
        # Move from active to inactive set
        job_data = redis_client.hgetall(job_key)
        room = job_data.get("room")
        if room:
            redis_client.smove(
                f"room:{room}:jobs:active", f"room:{room}:jobs:inactive", job_id
            )

        return True

    @staticmethod
    def fail_job(redis_client: Any, job_id: str, error: str) -> bool:
        """Mark job as failed.

        Args:
            redis_client: Redis client instance
            job_id: Job identifier
            error: Error message

        Returns:
            True if successful, False if job not found
        """
        job_key = f"job:{job_id}"

        if not redis_client.exists(job_key):
            return False

        now = datetime.utcnow().isoformat()
        redis_client.hset(
            job_key,
            mapping={
                "status": JobStatus.FAILED,
                "completed_at": now,
                "error": error,
            },
        )

        # Move from active to inactive set
        job_data = redis_client.hgetall(job_key)
        room = job_data.get("room")
        if room:
            redis_client.smove(
                f"room:{room}:jobs:active", f"room:{room}:jobs:inactive", job_id
            )

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
        job_key = f"job:{job_id}"
        job_data = redis_client.hgetall(job_key)

        if not job_data:
            return None

        # Parse JSON fields
        if job_data.get("data"):
            job_data["data"] = json.loads(job_data["data"])
        if job_data.get("result"):
            job_data["result"] = json.loads(job_data["result"])

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
        active_key = f"room:{room}:jobs:active"
        job_ids = redis_client.smembers(active_key)

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
        jobs_key = f"room:{room}:extension:{category}:{extension}:jobs"
        job_ids = redis_client.smembers(jobs_key)

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
        inactive_key = f"room:{room}:jobs:inactive"
        job_ids = redis_client.smembers(inactive_key)

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
        job_key = f"job:{job_id}"
        job_data = redis_client.hgetall(job_key)

        if not job_data:
            return False

        room = job_data.get("room")
        category = job_data.get("category")
        extension = job_data.get("extension")

        pipe = redis_client.pipeline()
        pipe.delete(job_key)

        if room:
            pipe.srem(f"room:{room}:jobs:active", job_id)
            pipe.srem(f"room:{room}:jobs:inactive", job_id)
            pipe.zrem(f"room:{room}:jobs:by_time", job_id)
            if category and extension:
                pipe.srem(
                    f"room:{room}:extension:{category}:{extension}:jobs", job_id
                )

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
            Dict with queued, running, completed, failed counts
        """
        jobs = JobManager.list_extension_jobs(redis_client, room, category, extension)

        stats = {
            "queued": 0,
            "running": 0,
            "completed": 0,
            "failed": 0,
        }

        for job in jobs:
            status = job.get("status")
            if status in stats:
                stats[status] += 1

        return stats
