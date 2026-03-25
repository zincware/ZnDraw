# tests/test_serve.py
"""Tests for JobManager auto-serve: background claim loop, lifecycle wrapping, wait()."""

import threading
from typing import Any, ClassVar
from unittest.mock import MagicMock

from zndraw_joblib.client import (
    Category,
    ClaimedTask,
    Extension,
    JobManager,
)
from zndraw_joblib.events import ProviderRequest, TaskAvailable
from zndraw_joblib.schemas import TaskResponse


class _Ext(Extension):
    """Concrete extension for serve tests."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self, vis: Any, **kwargs: Any) -> None:
        pass


# ---------------------------------------------------------------------------
# SIO handler registration
# ---------------------------------------------------------------------------


def test_task_available_handler_registered_in_init(mock_client_api, client):
    """TaskAvailable SIO handler is registered in __init__."""
    mock_tsio = MagicMock()
    JobManager(mock_client_api(client), tsio=mock_tsio)

    on_events = [c[0][0] for c in mock_tsio.on.call_args_list]
    assert TaskAvailable in on_events


def test_provider_request_handler_registered_in_init(mock_client_api, client):
    """ProviderRequest SIO handler is registered in __init__."""
    mock_tsio = MagicMock()
    JobManager(mock_client_api(client), tsio=mock_tsio)

    on_events = [c[0][0] for c in mock_tsio.on.call_args_list]
    assert ProviderRequest in on_events


# ---------------------------------------------------------------------------
# Background threads lifecycle
# ---------------------------------------------------------------------------


def test_no_threads_before_registration(mock_client_api, client):
    """No background threads start before first register() call."""
    manager = JobManager(mock_client_api(client), execute=lambda t: None)
    assert not manager._threads_started
    manager.disconnect()


def test_threads_start_on_first_register(mock_client_api, client):
    """Background threads start after first register()."""
    manager = JobManager(mock_client_api(client), execute=lambda t: None)

    @manager.register
    class Job(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    assert manager._threads_started
    manager.disconnect()


def test_threads_start_on_register_provider(mock_client_api, fs_provider, client):
    """Background threads start after register_provider()."""
    manager = JobManager(mock_client_api(client))
    manager.register_provider(fs_provider, name="local", handler=object())
    assert manager._threads_started
    manager.disconnect()


def test_no_claim_loop_without_execute(mock_client_api, client):
    """Without execute callback, no claim loop runs (manual mode)."""
    manager = JobManager(mock_client_api(client))

    @manager.register
    class Job(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    # Only heartbeat thread should exist — no claim loop
    assert len(manager._threads) == 1

    # Submit a task
    client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Job",
        json={"payload": {}},
    )

    # Long-poll for 1s — if a claim loop existed, it would have claimed by now
    resp = client.get("/v1/joblib/rooms/room_1/tasks")
    task_id = resp.json()["items"][0]["id"]
    resp = client.get(
        f"/v1/joblib/tasks/{task_id}",
        headers={"Prefer": "wait=1"},
    )
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "pending"

    # Manual claim still works
    claimed = manager.claim()
    assert claimed is not None
    manager.disconnect()


# ---------------------------------------------------------------------------
# Auto claim-execute lifecycle
# ---------------------------------------------------------------------------


def test_auto_claim_executes_task(mock_client_api, threadsafe_client):
    """Background loop claims task and calls execute callback."""
    executed = threading.Event()
    received: list[ClaimedTask] = []

    def execute(task: ClaimedTask) -> None:
        received.append(task)
        executed.set()

    manager = JobManager(
        mock_client_api(threadsafe_client), execute=execute, polling_interval=0.01
    )

    @manager.register
    class AutoJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER
        value: int = 0

    # Submit a task
    threadsafe_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:AutoJob",
        json={"payload": {"value": 42}},
    )

    # Wait for background thread to process it
    assert executed.wait(timeout=5.0), "Task was not executed within timeout"

    # Verify callback received the right data
    assert len(received) == 1
    assert received[0].extension.value == 42
    assert received[0].room_id == "room_1"

    # disconnect() joins threads, ensuring complete() finishes before we check DB
    task_id = received[0].task_id
    manager.disconnect()

    # Now safe to check task status — no concurrent DB access
    resp = threadsafe_client.get(f"/v1/joblib/tasks/{task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "completed"


def test_auto_execute_failure_marks_task_failed(mock_client_api, threadsafe_client):
    """If execute callback raises, task is marked FAILED with error."""
    executed = threading.Event()
    received: list[ClaimedTask] = []

    def execute(task: ClaimedTask) -> None:
        received.append(task)
        executed.set()
        raise RuntimeError("boom")

    manager = JobManager(
        mock_client_api(threadsafe_client), execute=execute, polling_interval=0.01
    )

    @manager.register
    class FailJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    threadsafe_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:FailJob",
        json={"payload": {}},
    )

    assert executed.wait(timeout=5.0), "Task was not executed within timeout"

    # disconnect() joins threads, ensuring fail() finishes before we check DB
    manager.disconnect()

    # Verify task is marked failed on the server
    resp = threadsafe_client.get(f"/v1/joblib/tasks/{received[0].task_id}")
    task = resp.json()
    assert task["status"] == "failed"
    assert "boom" in task["error"]


# ---------------------------------------------------------------------------
# disconnect() and wait()
# ---------------------------------------------------------------------------


def test_disconnect_stops_background_threads(mock_client_api, client):
    """disconnect() signals background threads to stop."""
    manager = JobManager(
        mock_client_api(client), execute=lambda t: None, polling_interval=0.01
    )

    @manager.register
    class StopJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    assert manager._threads_started
    manager.disconnect()
    assert manager._stop.is_set()


def test_wait_blocks_until_disconnect(mock_client_api, client):
    """wait() blocks until disconnect() is called from another thread."""
    manager = JobManager(mock_client_api(client))

    @manager.register
    class WaitJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    unblocked = threading.Event()

    def wait_thread():
        manager.wait()
        unblocked.set()

    t = threading.Thread(target=wait_thread, daemon=True)
    t.start()

    # wait() should be blocking
    assert not unblocked.wait(timeout=0.2)

    # disconnect() should unblock wait()
    manager.disconnect()
    assert unblocked.wait(timeout=2.0), "wait() did not unblock after disconnect()"


def test_context_manager_stops_threads(mock_client_api, client):
    """Exiting context manager stops background threads."""
    with JobManager(
        mock_client_api(client), execute=lambda t: None, polling_interval=0.01
    ) as manager:

        @manager.register
        class CtxJob(_Ext):
            category: ClassVar[Category] = Category.MODIFIER

        assert manager._threads_started

    # After exit, threads should be stopped
    assert manager._stop.is_set()


# ---------------------------------------------------------------------------
# E2E: full lifecycle with both jobs and providers
# ---------------------------------------------------------------------------


def test_e2e_job_and_provider_lifecycle(
    mock_client_api, fs_provider, threadsafe_client
):
    """Full e2e: register job + provider, submit task → auto-executed,
    read provider → dispatched and cached."""
    executed = threading.Event()
    received: list[ClaimedTask] = []
    mock_tsio = MagicMock()

    def execute(task: ClaimedTask) -> None:
        received.append(task)
        executed.set()

    mock_handler = MagicMock()
    mock_handler.list_dir.return_value = [{"name": "a.xyz"}]

    manager = JobManager(
        mock_client_api(threadsafe_client),
        tsio=mock_tsio,
        execute=execute,
        polling_interval=0.01,
    )

    # 1. Register a job
    @manager.register
    class Rotate(_Ext):
        category: ClassVar[Category] = Category.MODIFIER
        angle: float = 0.0

    # 2. Register a provider
    manager.register_provider(fs_provider, name="local", handler=mock_handler)

    # 3. Submit a task — auto-claimed and executed by background thread
    threadsafe_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {"angle": 90.0}},
    )
    assert executed.wait(timeout=5.0), "Task not auto-executed"
    assert received[0].extension.angle == 90.0

    # 4. Read provider — dispatch (immediate timeout), then simulate SIO
    from zndraw_joblib.dependencies import request_hash as compute_hash

    params = {"path": "/data"}
    rhash = compute_hash(params)

    read_resp = threadsafe_client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
        params=params,
        headers={"Prefer": "wait=0"},
    )
    assert read_resp.status_code == 504

    # Get the ProviderRequest handler from __init__
    pr_calls = [c for c in mock_tsio.on.call_args_list if c[0][0] is ProviderRequest]
    pr_callback = pr_calls[0][0][1]

    # Simulate server dispatching the event
    event = ProviderRequest.from_dict_params(
        request_id=rhash,
        provider_name="@global:filesystem:local",
        params=params,
    )
    pr_callback(event)

    # 5. Verify cached result
    cached_resp = threadsafe_client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
        params=params,
    )
    assert cached_resp.status_code == 200
    assert cached_resp.json() == [{"name": "a.xyz"}]

    # 6. disconnect() joins threads — ensures complete() finishes
    task_id = received[0].task_id
    manager.disconnect()

    # 7. Verify task reached completed (exclusive DB access after thread join)
    resp = threadsafe_client.get(f"/v1/joblib/tasks/{task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "completed"
