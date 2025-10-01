"""Worker statistics for ZnDraw extensions."""

from dataclasses import dataclass
from typing import Any

from .redis_keys import ExtensionKeys


@dataclass
class WorkerStats:
    """Statistics about workers for an extension."""

    idle_count: int
    progressing_count: int
    queue_length: int

    @property
    def total_workers(self) -> int:
        """Total number of workers (idle + progressing)."""
        return self.idle_count + self.progressing_count

    @classmethod
    def fetch(cls, redis_client: Any, keys: ExtensionKeys) -> "WorkerStats":
        """Fetch current worker statistics from Redis.

        Args:
            redis_client: Redis client instance
            keys: ExtensionKeys with the relevant Redis keys

        Returns:
            WorkerStats with current counts
        """
        idle = int(redis_client.scard(keys.idle_workers))
        progressing = int(redis_client.scard(keys.progressing_workers))
        queued = int(redis_client.llen(keys.queue))
        return cls(
            idle_count=idle, progressing_count=progressing, queue_length=queued
        )

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "idleWorkers": self.idle_count,
            "progressingWorkers": self.progressing_count,
            "queueLength": self.queue_length,
        }
