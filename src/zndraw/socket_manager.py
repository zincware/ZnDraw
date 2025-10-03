import dataclasses
import logging
import typing as t
import warnings

import socketio

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw

log = logging.getLogger(__name__)


@dataclasses.dataclass
class SocketIOLock:
    """A client-side context manager for a distributed lock via Socket.IO."""

    sio: socketio.Client
    target: str

    def acquire(self, timeout: float = 60) -> bool:
        """
        Acquire a lock for the specific target.
        Waits for the server's confirmation.
        """
        payload = {"target": self.target}
        # sio.call is inherently blocking, so it waits for the server's response.
        response = self.sio.call("lock:acquire", payload, timeout=timeout)
        return response and response.get("success", False)

    def release(self) -> bool:
        """Release the lock."""
        payload = {"target": self.target}
        response = self.sio.call("lock:release", payload, timeout=10)
        return response and response.get("success", False)

    def __enter__(self):
        if not self.acquire():
            raise RuntimeError(f"Failed to acquire lock for target '{self.target}'")
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
        self.sio.on("len_frames", self._on_len_frames_update)
        self.sio.on("invalidate", self._on_invalidate)
        self.sio.on("queue:update", self._on_queue_update)
        self.sio.on("frame_selection:update", self._on_frame_selection_update)
        self.sio.on("bookmarks:update", self._on_bookmarks_update)
        self.sio.on("frames:invalidate", self._on_frames_invalidate)

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

    def _on_len_frames_update(self, data):
        if "count" in data:
            self.zndraw._len = data["count"]

    def _on_selection_update(self, data):
        if "indices" in data:
            self.zndraw._selection = frozenset(data["indices"])

    def _on_frame_selection_update(self, data):
        if "indices" in data:
            self.zndraw._frame_selection = frozenset(data["indices"])

    def _on_bookmarks_update(self, data):
        if "bookmarks" in data:
            self.zndraw._bookmarks = {int(k): v for k, v in data["bookmarks"].items()}

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
