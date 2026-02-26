"""Shared fixtures and test extensions for worker E2E tests."""

import threading
from collections.abc import Callable
from typing import Any, ClassVar

import pytest
from zndraw_joblib.client import Category, Extension
from zndraw_joblib.schemas import JobSummary, TaskResponse

from zndraw import ZnDraw

# =============================================================================
# Test Extension Classes
# =============================================================================


class _Echo(Extension):
    """Modifier that writes a bookmark to prove room correctness."""

    category: ClassVar[Category] = Category.MODIFIER
    value: str = "hello"

    def run(self, vis: Any, **kwargs: Any) -> None:
        vis.bookmarks[0] = self.value


_Echo.__name__ = "Echo"


class _EchoAnalysis(Extension):
    """Analysis that writes a figure to prove room correctness."""

    category: ClassVar[Category] = Category.ANALYSIS
    value: float = 0.0

    def run(self, vis: Any, **kwargs: Any) -> None:
        import plotly.graph_objects as go

        fig = go.Figure(layout={"title": str(self.value)})
        vis.figures["echo_analysis"] = fig


_EchoAnalysis.__name__ = "EchoAnalysis"


class _FailingExtension(Extension):
    """Extension that always raises an error."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self, vis: Any, **kwargs: Any) -> None:
        raise RuntimeError("intentional failure")


_FailingExtension.__name__ = "FailingExtension"


# =============================================================================
# Extension class fixtures
# =============================================================================


@pytest.fixture
def Echo() -> type[_Echo]:
    """The Echo test extension class."""
    return _Echo


@pytest.fixture
def EchoAnalysis() -> type[_EchoAnalysis]:
    """The EchoAnalysis test extension class."""
    return _EchoAnalysis


@pytest.fixture
def FailingExtension() -> type[_FailingExtension]:
    """The FailingExtension test extension class."""
    return _FailingExtension


# =============================================================================
# Helper fixtures (return callables)
# =============================================================================


@pytest.fixture
def wait_for_task() -> Callable[..., TaskResponse]:
    """Long-poll GET until task reaches terminal status."""

    def _wait(vis: ZnDraw, task_id: str, *, timeout: int = 30) -> TaskResponse:
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/tasks/{task_id}",
            headers={**vis.api.get_headers(), "Prefer": f"wait={timeout}"},
            timeout=timeout + 5,
        )
        resp.raise_for_status()
        return TaskResponse.model_validate(resp.json())

    return _wait


@pytest.fixture
def get_job_list() -> Callable[..., list[JobSummary]]:
    """GET job list for a room (read-only assertion helper)."""

    def _get(vis: ZnDraw, room_id: str = "@global") -> list[JobSummary]:
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/rooms/{room_id}/jobs",
            headers=vis.api.get_headers(),
        )
        resp.raise_for_status()
        return [JobSummary.model_validate(j) for j in resp.json()["items"]]

    return _get


@pytest.fixture
def get_task() -> Callable[..., TaskResponse]:
    """GET a single task by ID (read-only assertion helper)."""

    def _get(vis: ZnDraw, task_id: str) -> TaskResponse:
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/tasks/{task_id}",
            headers=vis.api.get_headers(),
        )
        resp.raise_for_status()
        return TaskResponse.model_validate(resp.json())

    return _get


@pytest.fixture
def run_worker_loop() -> Callable[..., tuple[threading.Thread, threading.Event]]:
    """Start a real worker listen() loop in a background thread.

    Returns (thread, stop_event). Call stop.set(); thread.join(timeout=5)
    to shut down.
    """

    def _start(
        worker: ZnDraw, server_url: str
    ) -> tuple[threading.Thread, threading.Event]:
        stop = threading.Event()

        def _loop() -> None:
            try:
                for task in worker.jobs.listen(polling_interval=0.5, stop_event=stop):
                    worker.jobs.start(task)
                    vis = ZnDraw(url=server_url, room=task.room_id)
                    try:
                        task.extension.run(vis)
                    except Exception as e:
                        worker.jobs.fail(task, str(e))
                    else:
                        worker.jobs.complete(task)
                    finally:
                        vis.disconnect()
            except Exception:  # noqa: S110
                pass  # Worker disconnected — stop gracefully

        thread = threading.Thread(target=_loop, daemon=True)
        thread.start()
        return thread, stop

    return _start
