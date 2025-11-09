import dataclasses
import threading
import warnings
import logging
import typing as t

log = logging.getLogger(__name__)

if t.TYPE_CHECKING:
    from src.zndraw.api_manager import APIManager


@dataclasses.dataclass
class ZnDrawLock:
    """A client-side context manager for a distributed lock via REST API.

    Simplified design: delegates all REST calls to APIManager.
    Server controls TTL and refresh interval.
    No client-side re-entrancy - nested locks are not supported.

    Usage:
        with vis.get_lock(msg="Uploading trajectory data") as lock:
            vis.extend(frames)
            lock.update_msg("Processing chunk 5/10")
    """

    api: "APIManager"
    target: str
    msg: str | None = None

    # Server-provided values (set on acquire)
    _lock_token: str | None = dataclasses.field(default=None, init=False)
    _ttl: int | None = dataclasses.field(default=None, init=False)
    _refresh_interval: int | None = dataclasses.field(default=None, init=False)
    _refresh_thread: threading.Thread | None = dataclasses.field(
        default=None, init=False
    )
    _refresh_stop: threading.Event = dataclasses.field(
        default_factory=threading.Event, init=False
    )
    _is_held: bool = dataclasses.field(default=False, init=False)

    def _refresh_lock_periodically(self):
        """Background thread that refreshes the lock periodically."""
        while not self._refresh_stop.is_set():
            if self._refresh_stop.wait(timeout=self._refresh_interval):
                break

            try:
                # Refresh without updating message (requires lock token)
                if self._lock_token:
                    self.api.lock_refresh(self.target, self._lock_token)
                    log.debug(
                        f"Lock refreshed for target '{self.target}' with TTL {self._ttl}s"
                    )
                else:
                    log.error(f"Cannot refresh lock for target '{self.target}': no lock token")
                    break
            except Exception as e:
                log.error(f"Error refreshing lock for target '{self.target}': {e}")

    def acquire(self, timeout: float = 60) -> bool:
        """Acquire the lock via APIManager.

        Nested locking is not supported - attempting to acquire an already-held
        lock will raise RuntimeError.

        Parameters
        ----------
        timeout : float
            Request timeout (not used, kept for compatibility)

        Returns
        -------
        bool
            True if lock acquired successfully

        Raises
        ------
        RuntimeError
            If attempting to acquire an already-held lock or if acquisition fails
        """
        if self._is_held:
            raise RuntimeError(
                f"Lock for target '{self.target}' is already held by this instance. "
                "Nested locking is not supported."
            )

        try:
            # APIManager handles authentication and request
            response = self.api.lock_acquire(self.target, msg=self.msg)

            success = response.get("success", False)
            if success:
                self._is_held = True
                self._lock_token = response.get("lockToken")
                self._ttl = response.get("ttl", 60)
                self._refresh_interval = response.get("refreshInterval", 30)

                # Start refresh thread
                self._refresh_stop.clear()
                self._refresh_thread = threading.Thread(
                    target=self._refresh_lock_periodically,
                    daemon=True,
                    name=f"lock-refresh-{self.target}",
                )
                self._refresh_thread.start()

            return success
        except RuntimeError:
            raise  # Re-raise lock acquisition failures
        except Exception as e:
            raise RuntimeError(f"Failed to acquire lock: {e}") from e

    def refresh(self, msg: str | None = None) -> bool:
        """Refresh lock TTL and optionally update message.

        Parameters
        ----------
        msg : str | None
            Optional updated message

        Returns
        -------
        bool
            True if successfully refreshed
        """
        if not self._lock_token:
            log.warning(f"Cannot refresh lock for target '{self.target}': no lock token")
            return False

        try:
            response = self.api.lock_refresh(self.target, self._lock_token, msg=msg)
            if msg is not None:
                self.msg = msg  # Update instance message
            return response.get("success", False)
        except Exception as e:
            log.warning(f"Failed to refresh lock: {e}")
            return False

    def update_msg(self, msg: str) -> bool:
        """Update lock message (convenience method).

        Parameters
        ----------
        msg : str
            New message to set

        Returns
        -------
        bool
            True if successfully updated
        """
        return self.refresh(msg=msg)

    def release(self) -> bool:
        """Release the lock and stop refresh thread.

        Returns
        -------
        bool
            True if successfully released
        """
        if not self._is_held:
            warnings.warn(
                f"Attempting to release lock for '{self.target}' but not held"
            )
            return False

        # Stop refresh thread
        self._refresh_stop.set()
        if self._refresh_thread and self._refresh_thread.is_alive():
            self._refresh_thread.join(timeout=2)

        try:
            if self._lock_token:
                response = self.api.lock_release(self.target, self._lock_token)
            else:
                log.warning(f"Cannot release lock for target '{self.target}': no lock token")
                self._is_held = False
                return False

            self._is_held = False
            return response.get("success", False)
        except Exception as e:
            log.warning(f"Failed to release lock: {e}")
            self._is_held = False
            return False

    def __enter__(self):
        if not self.acquire():
            raise RuntimeError(f"Failed to acquire lock for target '{self.target}'")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.release():
            warnings.warn(
                f"Failed to release lock for target '{self.target}'. It may have expired."
            )

