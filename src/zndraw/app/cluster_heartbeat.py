"""Cluster heartbeat for detecting stale workers after server restart.

This module provides a heartbeat mechanism to detect when workers were
registered with a previous server instance that has since restarted.

The heartbeat key stores a timestamp indicating when the current "cluster session"
began. Workers registered before this timestamp are considered stale and should
be cleaned up.
"""

import logging
import threading
import time
import typing as t
from dataclasses import dataclass, field

from .redis_keys import ClusterKeys

log = logging.getLogger(__name__)


_MIN_REFRESH_INTERVAL = 5  # Minimum refresh interval to prevent tight loops


@dataclass
class ClusterHeartbeat:
    """Manages cluster heartbeat for stale worker detection.

    The heartbeat key stores the timestamp of when the current "cluster session" began.
    Workers registered before this timestamp are stale (from a previous server lifecycle).

    Parameters
    ----------
    redis_client
        Redis client instance
    ttl_seconds : int
        TTL for the heartbeat key in seconds. Default is 60.
        The refresh interval is automatically set to half this value (minimum 5s).
    """

    redis_client: t.Any
    ttl_seconds: int = 60
    _stop_event: threading.Event = field(default_factory=threading.Event, init=False)
    _thread: threading.Thread | None = field(default=None, init=False)
    _start_lock: threading.Lock = field(default_factory=threading.Lock, init=False)

    @property
    def refresh_interval(self) -> int:
        """Refresh interval is half the TTL, minimum 5 seconds."""
        return max(self.ttl_seconds // 2, _MIN_REFRESH_INTERVAL)

    def start(self) -> float:
        """Start the heartbeat background thread.

        If no heartbeat exists (first startup or all servers were down for > TTL),
        creates a new one with the current timestamp. Otherwise, uses the existing
        heartbeat and just starts refreshing it.

        Thread-safe: multiple calls will return the existing heartbeat timestamp
        if a thread is already running.

        Returns
        -------
        float
            The current heartbeat timestamp.
        """
        with self._start_lock:
            # Check if already running
            if self._thread is not None and self._thread.is_alive():
                timestamp = self.redis_client.get(ClusterKeys.HEARTBEAT)
                return float(timestamp) if timestamp else time.time()

            # Reset stop event for potential restart after stop()
            self._stop_event.clear()

            # Try to get existing heartbeat (another server instance may have set it)
            existing = self.redis_client.get(ClusterKeys.HEARTBEAT)

            if existing is None:
                # No heartbeat exists - this means either:
                # 1. First startup ever
                # 2. All servers were down for longer than TTL
                timestamp = time.time()
                self.redis_client.set(
                    ClusterKeys.HEARTBEAT, timestamp, ex=self.ttl_seconds
                )
                log.debug(
                    f"Cluster heartbeat initialized: {timestamp} "
                    f"(TTL: {self.ttl_seconds}s, refresh: {self.refresh_interval}s)"
                )
            else:
                timestamp = float(existing)
                # Refresh the TTL to keep it alive
                self.redis_client.expire(ClusterKeys.HEARTBEAT, self.ttl_seconds)
                log.debug(
                    f"Using existing cluster heartbeat: {timestamp} "
                    f"(TTL: {self.ttl_seconds}s)"
                )

            # Start refresh thread
            self._thread = threading.Thread(
                target=self._refresh_loop, daemon=True, name="cluster-heartbeat"
            )
            self._thread.start()

            return timestamp

    def stop(self) -> None:
        """Stop the heartbeat refresh thread."""
        self._stop_event.set()
        if self._thread:
            self._thread.join(timeout=5)
            self._thread = None

    def _refresh_loop(self) -> None:
        """Background loop that refreshes the heartbeat TTL."""
        while not self._stop_event.is_set():
            try:
                # Refresh TTL without changing the value
                self.redis_client.expire(ClusterKeys.HEARTBEAT, self.ttl_seconds)
                log.debug(f"Refreshed cluster heartbeat TTL to {self.ttl_seconds}s")
            except Exception as e:
                # Log warning but don't crash - Redis might have a brief hiccup
                log.warning(f"Failed to refresh cluster heartbeat (will retry): {e}")

            # Wait for refresh interval or until stopped
            self._stop_event.wait(self.refresh_interval)

    @classmethod
    def get_heartbeat_timestamp(cls, redis_client: t.Any) -> float | None:
        """Get current heartbeat timestamp.

        Parameters
        ----------
        redis_client
            Redis client instance

        Returns
        -------
        float | None
            The heartbeat timestamp, or None if no heartbeat exists.
            None should not happen during normal operation.
        """
        value = redis_client.get(ClusterKeys.HEARTBEAT)
        return float(value) if value else None

    @classmethod
    def is_worker_stale(
        cls, redis_client: t.Any, worker_registration_timestamp: float | str
    ) -> bool:
        """Check if a worker is stale based on its registration timestamp.

        A worker is stale if it was registered before the current cluster heartbeat.
        This indicates the worker was registered with a previous server instance
        that has since restarted.

        Parameters
        ----------
        redis_client
            Redis client instance
        worker_registration_timestamp : float | str
            The timestamp when the worker was registered

        Returns
        -------
        bool
            True if the worker is stale, False otherwise.
            Returns False if no heartbeat exists (can't determine staleness).
        """
        heartbeat_ts = cls.get_heartbeat_timestamp(redis_client)
        if heartbeat_ts is None:
            # No heartbeat - can't determine staleness, assume not stale
            return False

        reg_ts = float(worker_registration_timestamp)
        return reg_ts < heartbeat_ts
