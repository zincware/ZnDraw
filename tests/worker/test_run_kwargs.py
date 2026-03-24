"""E2E tests for run_kwargs support in register_job / register_extension.

Covers:
- register_job(cls, run_kwargs={...}) forwards to JobManager
- register_extension(cls, run_kwargs={...}) forwards through deprecated shim
- _execute_task unpacks run_kwargs into extension.run()
"""

import threading
import warnings
from typing import Any, ClassVar

import pytest
from zndraw_joblib.client import Category, Extension
from zndraw_joblib.schemas import TaskResponse

from zndraw import ZnDraw

# Shared mutable state to capture what run() actually received.
_captured_kwargs: dict[str, Any] = {}


class _KwargsCapture(Extension):
    """Extension that records **kwargs passed to run()."""

    category: ClassVar[Category] = Category.MODIFIER
    value: str = ""

    def run(self, vis: Any, **kwargs: Any) -> None:
        _captured_kwargs.clear()
        _captured_kwargs.update(kwargs)
        vis.bookmarks[0] = "captured"


_KwargsCapture.__name__ = "KwargsCapture"


@pytest.fixture(autouse=True)
def _clear_captured():
    """Reset captured kwargs before each test."""
    _captured_kwargs.clear()


@pytest.fixture
def KwargsCapture():  # noqa: N802
    return _KwargsCapture


# =============================================================================
# register_job accepts run_kwargs
# =============================================================================


def test_register_job_with_run_kwargs(server, KwargsCapture, get_job_list):
    """register_job(cls, run_kwargs={...}) should register and store kwargs."""
    model = object()  # non-serializable, worker-local
    worker = ZnDraw(url=server)
    try:
        worker.register_job(KwargsCapture, run_kwargs={"model": model})
        jobs = get_job_list(worker, room_id=worker.room)
        names = {j.name for j in jobs}
        assert "KwargsCapture" in names
        # Verify joblib stored them
        full_name = f"{worker.room}:modifiers:KwargsCapture"
        assert worker.jobs.get_run_kwargs(full_name) == {"model": model}
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# register_extension shim forwards run_kwargs
# =============================================================================


def test_register_extension_with_run_kwargs(server, KwargsCapture, get_job_list):
    """Deprecated register_extension(cls, run_kwargs={...}) should forward."""
    model = object()
    worker = ZnDraw(url=server)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            worker.register_extension(KwargsCapture, run_kwargs={"model": model})
            assert any(issubclass(x.category, DeprecationWarning) for x in w)
        jobs = get_job_list(worker, room_id=worker.room)
        names = {j.name for j in jobs}
        assert "KwargsCapture" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# End-to-end: run_kwargs delivered to run()
# =============================================================================


def _wait_for_task(vis: ZnDraw, task_id: str, *, timeout: int = 30) -> TaskResponse:
    """Long-poll GET until task reaches terminal status."""
    resp = vis.api.http.get(
        f"{vis.api.base_url}/v1/joblib/tasks/{task_id}",
        headers={**vis.api.get_headers(), "Prefer": f"wait={timeout}"},
        timeout=timeout + 5,
    )
    resp.raise_for_status()
    return TaskResponse.model_validate(resp.json())


def test_execute_task_unpacks_run_kwargs(server, KwargsCapture):
    """run_kwargs passed at registration should arrive in run(**kwargs)."""
    sentinel = {"key": "from_registration"}

    worker = ZnDraw(url=server)
    try:
        worker.register_job(KwargsCapture, run_kwargs={"extra": sentinel})

        # Submit a task from a separate client
        submitter = ZnDraw(url=server, room=worker.room)
        try:
            full_name = f"{worker.room}:modifiers:KwargsCapture"
            resp = submitter.api.http.post(
                f"{submitter.api.base_url}/v1/joblib/rooms/{worker.room}"
                f"/tasks/{full_name}",
                headers=submitter.api.get_headers(),
                json={"payload": {"value": "test"}},
            )
            assert resp.status_code == 202
            task_id = resp.json()["id"]

            # Run worker loop manually with run_kwargs support
            stop = threading.Event()

            def _loop() -> None:
                try:
                    for task in worker.jobs.listen(
                        polling_interval=0.5, stop_event=stop
                    ):
                        worker.jobs.start(task)
                        vis = ZnDraw(url=server, room=task.room_id)
                        try:
                            task.extension.run(vis, **task.run_kwargs)
                        except Exception as e:  # noqa: BLE001
                            worker.jobs.fail(task, str(e))
                        else:
                            worker.jobs.complete(task)
                        finally:
                            vis.disconnect()
                except Exception:  # noqa: BLE001, S110
                    pass

            t = threading.Thread(target=_loop, daemon=True)
            t.start()

            result = _wait_for_task(submitter, task_id)
            assert result.status.value == "completed"

            # Verify run() received the run_kwargs
            assert "extra" in _captured_kwargs
            assert _captured_kwargs["extra"] == sentinel
        finally:
            stop.set()
            submitter.disconnect()
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
