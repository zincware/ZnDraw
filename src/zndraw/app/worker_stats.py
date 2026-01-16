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

        Filters out stale workers (registered before the current cluster heartbeat)
        to avoid counting workers from previous server instances.

        Args:
            redis_client: Redis client instance
            keys: ExtensionKeys with the relevant Redis keys

        Returns:
            WorkerStats with current counts
        """
        from .cluster_heartbeat import ClusterHeartbeat

        # Get heartbeat timestamp for stale worker filtering
        heartbeat_ts = ClusterHeartbeat.get_heartbeat_timestamp(redis_client)

        # Get all workers with their registration timestamps
        all_workers = redis_client.hgetall(keys.workers)

        # Filter out stale workers
        valid_workers = []
        for worker_id, reg_ts in all_workers.items():
            if heartbeat_ts is not None and float(reg_ts) < heartbeat_ts:
                # Worker is stale - skip it
                continue
            valid_workers.append(worker_id)

        # Count idle vs busy workers by checking capacity
        idle_count = 0
        busy_count = 0

        if valid_workers:
            pipe = redis_client.pipeline()
            for worker_id in valid_workers:
                capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
                pipe.get(capacity_key)
            capacities = pipe.execute()

            for capacity in capacities:
                # Skip workers without capacity key (likely stale)
                if capacity is None:
                    continue
                capacity_val = int(capacity)
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
