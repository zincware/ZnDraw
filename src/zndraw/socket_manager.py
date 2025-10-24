import dataclasses
import logging
import threading
import time
import typing as t
import warnings
import traceback

import socketio

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw

log = logging.getLogger(__name__)


@dataclasses.dataclass
class SocketIOLock:
    """A client-side context manager for a distributed lock via Socket.IO.

    This lock is re-entrant - the same client can acquire it multiple times
    and must release it the same number of times.

    Includes automatic lock renewal to prevent TTL expiration during long operations.

    Optional metadata can be sent to describe the lock's purpose:
        with vis.lock(msg="Uploading trajectory data"):
            vis.extend(frames)
    """

    sio: socketio.Client
    target: str
    ttl: int = 60  # TTL in seconds (must be <= 300 as validated by server)
    _lock_count: int = dataclasses.field(default=0, init=False)
    _refresh_thread: threading.Thread | None = dataclasses.field(default=None, init=False)
    _refresh_stop: threading.Event = dataclasses.field(default_factory=threading.Event, init=False)
    _pending_metadata: dict | None = dataclasses.field(default=None, init=False)

    def __call__(self, msg: str | None = None, metadata: dict | None = None) -> "SocketIOLock":
        """Set optional metadata for this lock acquisition.

        Returns self to maintain singleton pattern and support re-entrant locking.

        Parameters
        ----------
        msg : str | None
            Human-readable message describing the lock purpose
        metadata : dict | None
            Additional metadata fields

        Returns
        -------
        SocketIOLock
            Returns self for use as context manager

        Examples
        --------
        >>> with vis.lock(msg="Uploading trajectory"):
        ...     vis.extend(frames)
        >>> with vis.lock(metadata={"step": 1, "total": 10}):
        ...     process_batch()
        """
        self._pending_metadata = {}
        if msg is not None:
            self._pending_metadata["msg"] = msg
        if metadata:
            self._pending_metadata.update(metadata)
        return self

    def _refresh_lock_periodically(self):
        """Background thread that refreshes the lock periodically to prevent TTL expiration.

        Refreshes at half the TTL interval to ensure the lock stays active.
        """
        refresh_interval = self.ttl / 2
        while not self._refresh_stop.is_set():
            # Wait for half the TTL or until stop signal
            if self._refresh_stop.wait(timeout=refresh_interval):
                break

            # Refresh the lock by calling lock:refresh
            try:
                payload = {"target": self.target, "ttl": self.ttl}
                response = self.sio.call("lock:refresh", payload, timeout=10)
                if response and response.get("success"):
                    log.debug(f"Lock refreshed for target '{self.target}' with TTL {self.ttl}s")
                else:
                    error_msg = response.get("error", "Unknown error") if response else "No response"
                    log.warning(f"Failed to refresh lock for target '{self.target}': {error_msg}")
            except Exception as e:
                log.error(f"Error refreshing lock for target '{self.target}': {e}")

    def _send_metadata(self):
        """Send lock metadata to server after lock acquisition.

        Logs warnings on failure but does not raise exceptions to avoid
        breaking the lock acquisition flow.
        """
        payload = {
            "target": self.target,
            "metadata": self._pending_metadata
        }
        try:
            response = self.sio.call("lock:msg", payload, timeout=5)
            if not (response and response.get("success")):
                log.warning(
                    f"Failed to send lock metadata for target '{self.target}': {response}"
                )
        except Exception as e:
            log.error(
                f"Error sending lock metadata for target '{self.target}': {e}",
                exc_info=True
            )

    def acquire(self, timeout: float = 60) -> bool:
        """
        Acquire a lock for the specific target.
        If already held by this client, increment the lock count.
        Starts a background thread to refresh the lock periodically.
        
        Args:
            timeout: Socket.IO call timeout (not the lock TTL)
        
        Returns:
            True if lock acquired successfully, False otherwise
        """
        # If we already hold the lock, just increment the count
        if self._lock_count > 0:
            self._lock_count += 1
            return True
        
        payload = {"target": self.target, "ttl": self.ttl}
        # sio.call is inherently blocking, so it waits for the server's response.
        response = self.sio.call("lock:acquire", payload, timeout=int(timeout))
        
        # Check if server returned an error
        if response and response.get("error"):
            raise ValueError(f"Failed to acquire lock: {response['error']}")
        
        success = response and response.get("success", False)
        if success:
            self._lock_count = 1
            # Start the refresh thread
            self._refresh_stop.clear()
            self._refresh_thread = threading.Thread(
                target=self._refresh_lock_periodically,
                daemon=True,
                name=f"lock-refresh-{self.target}"
            )
            self._refresh_thread.start()
        return success

    def release(self) -> bool:
        """Release the lock. Only actually releases when count reaches 0.
        Stops the refresh thread when the lock is fully released."""
        if self._lock_count == 0:
            warnings.warn(f"Attempting to release lock for '{self.target}' but not held")
            return False
        
        self._lock_count -= 1
        
        # Only actually release the lock when count reaches 0
        if self._lock_count == 0:
            # Stop the refresh thread first
            self._refresh_stop.set()
            if self._refresh_thread and self._refresh_thread.is_alive():
                self._refresh_thread.join(timeout=2)
            
            payload = {"target": self.target}
            response = self.sio.call("lock:release", payload, timeout=10)
            return response and response.get("success", False)
        
        return True  # Successfully decremented count

    def __enter__(self):
        if not self.acquire():
            raise RuntimeError(f"Failed to acquire lock for target '{self.target}'")

        # Send metadata if present (works for both initial and re-entrant)
        if self._pending_metadata:
            self._send_metadata()

        # Always clear to prevent stale metadata
        self._pending_metadata = None

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.release():
            warnings.warn(
                f"Failed to release lock for target '{self.target}'. It may have expired."
            )


class SocketManager:
    def __init__(self, zndraw_instance: "ZnDraw", join_token: str):
        self.zndraw = zndraw_instance
        self.join_token = join_token
        self.sio = socketio.Client()
        self._register_handlers()

    def _register_handlers(self):
        self.sio.on("connect", self._on_connect)
        self.sio.on("frame_update", self._on_frame_update)
        self.sio.on("selection:update", self._on_selection_update)
        self.sio.on("room:update", self._on_room_update)
        self.sio.on("invalidate", self._on_invalidate)
        self.sio.on("queue:update", self._on_queue_update)
        self.sio.on("frame_selection:update", self._on_frame_selection_update)
        self.sio.on("bookmarks:invalidate", self._on_bookmarks_invalidate)
        self.sio.on("frames:invalidate", self._on_frames_invalidate)
        self.sio.on("invalidate:geometry", self._on_geometry_invalidate)
        self.sio.on("invalidate:figure", self._on_figure_invalidate)

    def connect(self):
        if self.sio.connected:
            print("Already connected.")
            return
        # Connect with join token for authentication
        self.sio.connect(self.zndraw.url, auth={"token": self.join_token}, wait=True)

    def disconnect(self):
        if self.sio.connected:
            self.sio.disconnect()
            print("Disconnected.")

    @property
    def connected(self) -> bool:
        return self.sio.connected

    def _on_connect(self):
        log.debug(
            f"Joined room: '{self.zndraw.room}' with client ID: {self.zndraw.sid}"
        )
        for name, ext in self.zndraw._extensions.items():
            self.zndraw.api.register_extension(
                name=name,
                category=ext["extension"].category,
                schema=ext["extension"].model_json_schema(),
                client_id=self.zndraw.sid,
            )
        self._on_queue_update({})

    def _on_frame_update(self, data):
        if "frame" in data:
            self.zndraw._step = data["frame"]

    def _on_room_update(self, data):
        """Handle room:update events (consolidated room metadata updates)."""
        if "frameCount" in data:
            self.zndraw._len = data["frameCount"]

    def _on_selection_update(self, data):
        if "indices" in data:
            self.zndraw._selection = frozenset(data["indices"])

    def _on_frame_selection_update(self, data):
        if "indices" in data:
            self.zndraw._frame_selection = frozenset(data["indices"])

    def _on_bookmarks_invalidate(self, data):
        """Handle bookmark invalidation by refetching from server."""
        # Refetch all bookmarks from server to update local cache
        bookmarks = self.zndraw.api.get_all_bookmarks()
        self.zndraw._bookmarks = bookmarks

    def _on_geometry_invalidate(self, data):
        if key := data.get("key"):
            # Refresh geometries from server to get updated state
            response = self.zndraw.api.get_geometries()
            if response is not None:
                self.zndraw._geometries = response

    def _on_figure_invalidate(self, data):
        if key := data.get("key"):
            self.zndraw._figures.pop(key, None)

    def _on_queue_update(self, data: dict):
        print(f"Queue update received: {data}")
        if not self.zndraw.auto_pickup_jobs:
            return

        job_data = self.zndraw.api.get_next_job(self.zndraw.sid)
        if job_data and "jobId" in job_data:
            try:
                self._on_task_run(
                    data=job_data.get("data"),
                    extension=job_data.get("extension"),
                    category=job_data.get("category"),
                )
                self.zndraw.api.update_job_status(
                    job_id=job_data.get("jobId"),
                    status="completed",
                    worker_id=self.zndraw.sid,
                )
            except Exception as e:
                log.error(f"Error processing job {job_data.get('jobId')}: {e}")
                traceback.print_exc()
                self.zndraw.api.update_job_status(
                    job_id=job_data.get("jobId"),
                    status="failed",
                    error=str(e),
                    worker_id=self.zndraw.sid,
                )
            self._on_queue_update({})

    def _on_task_run(self, data: dict, extension: str, category: str):
        ext = self.zndraw._extensions[extension]["extension"]
        instance = ext(**(data))
        instance.run(
            self.zndraw, **(self.zndraw._extensions[extension]["run_kwargs"] or {})
        )

    def _on_invalidate(self, data: dict):
        if data["category"] == "settings":
            self.zndraw._settings.pop(data["extension"], None)

    def _on_frames_invalidate(self, data: dict):
        log.debug(f"Received cache invalidation event: {data}")
        if self.zndraw.cache is None:
            return

        operation = data.get("operation")

        if operation == "replace":
            idx = data.get("affectedIndex")
            if idx is not None:
                self.zndraw.cache.pop(idx, None)
        elif operation in ("insert", "delete", "bulk_replace"):
            from_idx = data.get("affectedFrom")
            if from_idx is not None:
                self.zndraw.cache.invalidate_from(from_idx)
        elif operation == "clear_all":
            self.zndraw.cache.clear()
        else:
            log.warning(
                "Unknown or broad invalidation event received. Clearing entire frame cache."
            )
            self.zndraw.cache.clear()
