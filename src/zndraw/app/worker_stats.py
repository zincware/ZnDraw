"""Worker statistics for ZnDraw extensions."""

from dataclasses import dataclass
from typing import Any

from .redis_keys import ExtensionKeys


@dataclass
class WorkerStats:
    """Statistics about workers for an extension."""

    idle_count: int
    busy_count: int
    queue_length: int

    @property
    def total_workers(self) -> int:
        """Total number of workers (idle + busy)."""
        return self.idle_count + self.busy_count

    @classmethod
    def fetch(cls, redis_client: Any, keys: ExtensionKeys) -> "WorkerStats":
        """Fetch current worker statistics from Redis.

        Args:
            redis_client: Redis client instance
            keys: ExtensionKeys with the relevant Redis keys

        Returns:
            WorkerStats with current counts
        """
        # Get all workers registered for this extension
        all_workers = redis_client.hkeys(keys.workers)
        total_workers = len(all_workers)

        # Count idle vs busy workers by checking capacity
        idle_count = 0
        busy_count = 0

        if all_workers:
            pipe = redis_client.pipeline()
            for worker_id in all_workers:
                capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
                pipe.get(capacity_key)
            capacities = pipe.execute()

            for capacity in capacities:
                capacity_val = int(capacity) if capacity else 0
                if capacity_val >= 1:
                    idle_count += 1
                else:
                    busy_count += 1

        pending = int(redis_client.zcard(keys.pending_jobs))
        return cls(idle_count=idle_count, busy_count=busy_count, queue_length=pending)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "idleWorkers": self.idle_count,
            "progressingWorkers": self.busy_count,
            "queueLength": self.queue_length,
        }
