"""Job object for tracking job progress."""

import time
from typing import Any

from zndraw.app.job_manager import JobStatus


class Job:
    """Represents a submitted job and allows tracking its progress.

    This object is returned by vis.run() and provides methods to monitor
    job execution, wait for completion, and retrieve results.

    In Jupyter notebooks, displays an iframe showing live progress.

    Parameters
    ----------
    job_id : str
        Unique job identifier
    url : str
        Server URL
    room : str
        Room ID
    api : APIManager
        API manager instance for making requests

    Examples
    --------
    >>> job = vis.run(MyExtension(param=42))
    >>> job.wait()  # Block until completion
    >>> result = job.get_result()
    """

    def __init__(self, job_id: str, url: str, room: str, api: Any, socket: Any = None):
        self.job_id = job_id
        self.url = url
        self.room = room
        self.api = api
        self.socket = socket
        self._cached_data: dict[str, Any] | None = None

    def refresh(self) -> dict[str, Any]:
        """Fetch latest job status from server.

        Returns
        -------
        dict
            Job details including status, timestamps, and result/error
        """
        self._cached_data = self.api.get_job(self.job_id)
        return self._cached_data

    @property
    def status(self) -> str:
        """Get current job status.

        Returns
        -------
        str
            One of: pending, assigned, processing, completed, failed
        """
        if self._cached_data is None:
            self.refresh()
        return self._cached_data.get("status", "unknown")

    def is_pending(self) -> bool:
        """Check if job is pending (waiting for worker)."""
        return self.status == JobStatus.PENDING

    def is_assigned(self) -> bool:
        """Check if job is assigned to a worker."""
        return self.status == JobStatus.ASSIGNED

    def is_processing(self) -> bool:
        """Check if job is currently being processed."""
        return self.status == JobStatus.PROCESSING

    def is_completed(self) -> bool:
        """Check if job completed successfully."""
        return self.status == JobStatus.COMPLETED

    def is_failed(self) -> bool:
        """Check if job failed."""
        return self.status == JobStatus.FAILED

    def is_done(self) -> bool:
        """Check if job is in a terminal state (completed or failed)."""
        return self.is_completed() or self.is_failed()

    def wait(self, timeout: float | None = None, poll_interval: float = 0.5) -> None:
        """Block until job completes or fails.

        Parameters
        ----------
        timeout : float | None
            Maximum time to wait in seconds. None means wait forever.
        poll_interval : float
            How often to check status in seconds

        Raises
        ------
        TimeoutError
            If timeout is reached before job completes
        """
        start_time = time.time()

        while not self.is_done():
            if timeout is not None and (time.time() - start_time) > timeout:
                raise TimeoutError(
                    f"Job {self.job_id} did not complete within {timeout} seconds"
                )

            # Use socketio's sleep if available to avoid blocking the event loop
            if self.socket is not None:
                self.socket.sio.sleep(poll_interval)
            else:
                time.sleep(poll_interval)
            self.refresh()

    def __repr__(self) -> str:
        """Terminal-friendly representation."""
        status = self.status
        return f"Job(id={self.job_id}, status={status})"

    # def _repr_html_(self) -> str:
    #     """Jupyter notebook representation using iframe.

    #     Displays live progress by embedding the server's job page.
    #     """
    #     try:
    #         from IPython.display import IFrame
    #     except ImportError:
    #         raise ImportError(
    #             "IPython is required for viewer display. Install with: uv add / pip install ipython"
    #         )
    #     iframe_url = f"{self.url}/job/{self.job_id}"
    #     return IFrame(src=iframe_url, width="100%", height=600)._repr_html_()
