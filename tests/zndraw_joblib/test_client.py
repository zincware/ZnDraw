# tests/test_client.py
"""Integration tests for the client SDK against the actual server."""

from typing import ClassVar
from unittest.mock import MagicMock

import pytest

from zndraw_joblib.client import (
    Category,
    Extension,
    JobManager,
)
from zndraw_joblib.events import JoinJobRoom, LeaveJobRoom, ProviderRequest
from zndraw_joblib.schemas import (
    JobResponse,
    JobSummary,
    PaginatedResponse,
    TaskResponse,
)


class ConcreteExtension(Extension):
    """Non-abstract extension base for tests that don't need run()."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self) -> None:
        pass


def test_category_enum():
    """Category enum should have correct values."""
    assert Category.MODIFIER.value == "modifiers"
    assert Category.SELECTION.value == "selections"
    assert Category.ANALYSIS.value == "analysis"


def test_extension_requires_category():
    """Extension subclasses must define category."""

    class ValidExtension(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        value: int = 0

    assert ValidExtension.category == Category.MODIFIER


@pytest.mark.parametrize(
    ("category", "expected_value"),
    [
        (Category.MODIFIER, "modifiers"),
        (Category.SELECTION, "selections"),
        (Category.ANALYSIS, "analysis"),
    ],
)
def test_job_manager_register_category(api, client, category, expected_value):
    """JobManager.register should create job with the given category."""
    manager = JobManager(api)

    # Dynamically create a class with the given category
    ext_cls = type(
        "TestExt",
        (ConcreteExtension,),
        {"__annotations__": {"category": ClassVar[Category]}, "category": category},
    )

    manager.register(ext_cls)

    full_name = f"@global:{expected_value}:TestExt"
    assert full_name in manager

    response = client.get("/v1/joblib/rooms/@global/jobs")
    assert response.status_code == 200
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    job_names = [j.full_name for j in page.items]
    assert full_name in job_names


def test_job_manager_register_with_room(api, client):
    """JobManager.register(room=...) should create room-specific job."""
    manager = JobManager(api)

    @manager.register(room="my_room")
    class PrivateJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        data: str = ""

    assert "my_room:modifiers:PrivateJob" in manager

    response = client.get("/v1/joblib/rooms/my_room/jobs")
    assert response.status_code == 200
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    job_names = [j.full_name for j in page.items]
    assert "my_room:modifiers:PrivateJob" in job_names


def test_job_manager_getitem_returns_class(api, client):
    """JobManager[full_name] should return the registered class."""
    manager = JobManager(api)

    @manager.register
    class MyJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        x: int = 0

    assert manager["@global:modifiers:MyJob"] is MyJob


def test_job_manager_len(api, client):
    """len(JobManager) should return number of registered jobs."""
    manager = JobManager(api)

    assert len(manager) == 0

    @manager.register
    class Job1(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    assert len(manager) == 1

    @manager.register
    class Job2(ConcreteExtension):
        category: ClassVar[Category] = Category.SELECTION

    assert len(manager) == 2


def test_job_manager_iter(api, client):
    """iter(JobManager) should iterate over registered job names."""
    manager = JobManager(api)

    @manager.register
    class JobA(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    @manager.register
    class JobB(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    names = list(manager)
    assert "@global:modifiers:JobA" in names
    assert "@global:modifiers:JobB" in names


def test_job_manager_schema_sent_to_server(api, client):
    """JobManager should send the Pydantic schema to the server."""
    manager = JobManager(api)

    @manager.register
    class Rotate(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        angle: float = 0.0
        axis: str = "z"

    response = client.get("/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate")
    assert response.status_code == 200
    job = JobResponse.model_validate(response.json())

    assert "properties" in job.schema_
    assert "angle" in job.schema_["properties"]
    assert "axis" in job.schema_["properties"]


def test_job_manager_listen_yields_extension_instance(api, client):
    """JobManager.listen() should yield Extension instances with payload data."""
    manager = JobManager(api)

    @manager.register
    class Rotate(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        angle: float = 0.0

    # Submit a task
    submit_resp = client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {"angle": 45.0}},
    )
    assert submit_resp.status_code == 202
    _ = TaskResponse.model_validate(submit_resp.json())

    # Listen should yield the Extension instance
    for claimed in manager.listen():
        assert isinstance(claimed.extension, Rotate)
        assert claimed.extension.angle == 45.0
        assert claimed.task_id is not None
        assert claimed.room_id == "room_1"
        break  # Only get one task


def test_job_manager_listen_returns_none_when_empty(api, client):
    """JobManager.listen() with timeout should return when no tasks available."""
    manager = JobManager(api)

    @manager.register
    class EmptyJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    # No tasks submitted - claim should return None
    claimed = manager.claim()
    assert claimed is None


def test_job_manager_claim_until_empty(api, client):
    """Calling claim repeatedly should return None when no more tasks."""
    manager = JobManager(api)

    @manager.register
    class BatchJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        index: int = 0

    # Submit 3 tasks
    for i in range(3):
        resp = client.post(
            "/v1/joblib/rooms/room_1/tasks/@global:modifiers:BatchJob",
            json={"payload": {"index": i}},
        )
        assert resp.status_code == 202

    # Claim all 3
    claimed_indices = []
    for _ in range(3):
        claimed = manager.claim()
        assert claimed is not None
        assert isinstance(claimed.extension, BatchJob)
        claimed_indices.append(claimed.extension.index)

    # Verify we got all 3 (order may vary due to FIFO)
    assert sorted(claimed_indices) == [0, 1, 2]

    # Fourth claim should return None
    assert manager.claim() is None


def test_job_manager_heartbeat(api, client):
    """JobManager.heartbeat() should update worker timestamp."""

    manager = JobManager(api)

    @manager.register
    class HeartbeatJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    # The worker_id was set during register - use that
    assert manager.worker_id is not None

    # First heartbeat
    manager.heartbeat()

    # Verify heartbeat was recorded by calling again
    response = client.patch(f"/v1/joblib/workers/{manager.worker_id}")
    assert response.status_code == 200


def test_job_manager_complete_workflow(api, client):
    """Test complete workflow: register, submit, claim, complete, verify empty."""
    manager = JobManager(api)

    @manager.register
    class ProcessData(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        input_file: str = ""
        output_file: str = ""

    # 1. Submit two tasks
    for i in range(2):
        resp = client.post(
            "/v1/joblib/rooms/room_1/tasks/@global:modifiers:ProcessData",
            json={
                "payload": {"input_file": f"in{i}.txt", "output_file": f"out{i}.txt"}
            },
        )
        task = TaskResponse.model_validate(resp.json())
        assert task.status.value == "pending"

    # 2. Claim first task
    claimed1 = manager.claim()
    assert claimed1 is not None
    assert isinstance(claimed1.extension, ProcessData)

    # 3. Update to running and complete
    manager.start(claimed1)
    manager.complete(claimed1)

    # 4. Verify completed status
    complete_resp = client.get(f"/v1/joblib/tasks/{claimed1.task_id}")
    completed = TaskResponse.model_validate(complete_resp.json())
    assert completed.status.value == "completed"

    # 5. Claim second task
    claimed2 = manager.claim()
    assert claimed2 is not None
    assert isinstance(claimed2.extension, ProcessData)

    # 6. Complete second task
    manager.start(claimed2)
    manager.complete(claimed2)

    # 6. No more tasks - claim should return None
    assert manager.claim() is None

    # 7. Validate finished tasks have completed_at
    assert completed.completed_at is not None


def test_job_manager_claimed_task_has_metadata(api, client):
    """ClaimedTask should include task_id, room_id, job_name, and extension."""
    manager = JobManager(api)

    @manager.register
    class MetadataJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        value: int = 42

    client.post(
        "/v1/joblib/rooms/test_room/tasks/@global:modifiers:MetadataJob",
        json={"payload": {"value": 99}},
    )

    claimed = manager.claim()
    assert claimed is not None
    assert claimed.task_id is not None
    assert claimed.room_id == "test_room"
    assert claimed.job_name == "@global:modifiers:MetadataJob"
    assert isinstance(claimed.extension, MetadataJob)
    assert claimed.extension.value == 99


def test_job_manager_submit(api, client):
    """JobManager.submit() should create a task via the server."""
    manager = JobManager(api)

    @manager.register
    class SubmitJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        value: int = 0

    # Submit a task
    task_id = manager.submit(SubmitJob(value=42), room="submit_room")
    assert task_id is not None

    # Verify the task exists on the server
    task_resp = client.get(f"/v1/joblib/tasks/{task_id}")
    assert task_resp.status_code == 200
    task = TaskResponse.model_validate(task_resp.json())
    assert task.payload == {"value": 42}
    assert task.room_id == "submit_room"
    assert task.job_name == "@global:modifiers:SubmitJob"


def test_job_manager_claim_raises_without_worker_id(api, client):
    """claim() should raise ValueError if worker_id is not set."""
    manager = JobManager(api)
    # Don't register any jobs (so _worker_id stays None)
    with pytest.raises(ValueError, match="Worker ID not set"):
        manager.claim()


def test_job_manager_heartbeat_raises_without_worker_id(api, client):
    """heartbeat() should raise ValueError if worker_id is not set."""
    manager = JobManager(api)
    with pytest.raises(ValueError, match="Worker ID not set"):
        manager.heartbeat()


def test_job_manager_register_emits_join_job_room(api, client):
    """register() should emit JoinJobRoom when tsio is provided."""
    mock_tsio = MagicMock()
    manager = JobManager(api, tsio=mock_tsio)

    @manager.register
    class Rotate(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        angle: float = 0.0

    mock_tsio.emit.assert_called_once()
    event = mock_tsio.emit.call_args[0][0]
    assert isinstance(event, JoinJobRoom)
    assert event.job_name == "@global:modifiers:Rotate"
    assert event.worker_id == str(manager.worker_id)


def test_job_manager_register_emits_join_for_each_job(api, client):
    """register() should emit JoinJobRoom for each registered job."""
    mock_tsio = MagicMock()
    manager = JobManager(api, tsio=mock_tsio)

    @manager.register
    class Job1(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    @manager.register
    class Job2(ConcreteExtension):
        category: ClassVar[Category] = Category.SELECTION

    assert mock_tsio.emit.call_count == 2
    events = [call[0][0] for call in mock_tsio.emit.call_args_list]
    worker_id = str(manager.worker_id)
    assert events[0] == JoinJobRoom(
        job_name="@global:modifiers:Job1", worker_id=worker_id
    )
    assert events[1] == JoinJobRoom(
        job_name="@global:selections:Job2", worker_id=worker_id
    )


def test_job_manager_register_no_tsio_no_emit(api, client):
    """register() should not break when tsio is None (default)."""
    manager = JobManager(api)

    @manager.register
    class NoTsioJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    assert "@global:modifiers:NoTsioJob" in manager


def test_job_manager_register_room_emits_correct_job_name(api, client):
    """register(room=...) should emit JoinJobRoom with room-scoped job name."""
    mock_tsio = MagicMock()
    manager = JobManager(api, tsio=mock_tsio)

    @manager.register(room="my_room")
    class PrivateJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    event = mock_tsio.emit.call_args[0][0]
    assert isinstance(event, JoinJobRoom)
    assert event.job_name == "my_room:modifiers:PrivateJob"
    assert event.worker_id == str(manager.worker_id)


def test_extension_cannot_be_instantiated_directly():
    """Extension base class cannot be instantiated (ABC)."""
    with pytest.raises(TypeError, match="Can't instantiate abstract class"):
        Extension()


def test_job_manager_disconnect_emits_leave_and_deletes_worker(api, client):
    """disconnect() should emit LeaveJobRoom for each job and DELETE the worker."""
    mock_tsio = MagicMock()
    manager = JobManager(api, tsio=mock_tsio)

    @manager.register
    class Job1(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    @manager.register
    class Job2(ConcreteExtension):
        category: ClassVar[Category] = Category.SELECTION

    worker_id = manager.worker_id
    assert worker_id is not None
    mock_tsio.reset_mock()

    manager.disconnect()

    # Should have emitted LeaveJobRoom for each job
    assert mock_tsio.emit.call_count == 2
    events = [call[0][0] for call in mock_tsio.emit.call_args_list]
    assert all(isinstance(e, LeaveJobRoom) for e in events)
    job_names = {e.job_name for e in events}
    assert job_names == {"@global:modifiers:Job1", "@global:selections:Job2"}

    # Local state should be cleared
    assert len(manager) == 0
    assert manager.worker_id is None

    # Worker should be deleted from server
    resp = client.patch(f"/v1/joblib/workers/{worker_id}")
    assert resp.status_code == 404


def test_job_manager_disconnect_no_tsio(api, client):
    """disconnect() should work without tsio (REST-only cleanup)."""
    manager = JobManager(api)

    @manager.register
    class RestJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER

    worker_id = manager.worker_id
    assert worker_id is not None

    manager.disconnect()

    assert len(manager) == 0
    assert manager.worker_id is None

    # Worker should be deleted
    resp = client.patch(f"/v1/joblib/workers/{worker_id}")
    assert resp.status_code == 404


def test_job_manager_disconnect_no_worker_id():
    """disconnect() should be safe to call when no worker_id is set."""
    api = MagicMock()
    manager = JobManager(api)
    manager.disconnect()  # Should not raise
    assert manager.worker_id is None


def _register_claim(api, client):
    """Helper: register a job, submit a task, claim it. Returns (manager, claimed)."""
    manager = JobManager(api)

    @manager.register
    class TaskJob(ConcreteExtension):
        category: ClassVar[Category] = Category.MODIFIER
        value: int = 0

    client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:TaskJob",
        json={"payload": {"value": 7}},
    )
    claimed = manager.claim()
    assert claimed is not None
    return manager, claimed


def test_job_manager_start_transitions_to_running(api, client):
    """manager.start(task) should transition task from CLAIMED to RUNNING."""
    manager, claimed = _register_claim(api, client)

    manager.start(claimed)

    resp = client.get(f"/v1/joblib/tasks/{claimed.task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "running"
    assert task.started_at is not None


def test_job_manager_complete_transitions_to_completed(api, client):
    """manager.complete(task) should transition task from RUNNING to COMPLETED."""
    manager, claimed = _register_claim(api, client)

    manager.start(claimed)
    manager.complete(claimed)

    resp = client.get(f"/v1/joblib/tasks/{claimed.task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "completed"
    assert task.completed_at is not None


def test_job_manager_fail_transitions_to_failed_with_error(api, client):
    """manager.fail(task, error) should transition to FAILED and store error."""
    manager, claimed = _register_claim(api, client)

    manager.start(claimed)
    manager.fail(claimed, "something broke")

    resp = client.get(f"/v1/joblib/tasks/{claimed.task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "failed"
    assert task.error == "something broke"
    assert task.completed_at is not None


def test_job_manager_fail_by_id(api, client):
    """manager.fail_by_id(task_id, error) should fail a task without a ClaimedTask."""
    manager, claimed = _register_claim(api, client)

    manager.start(claimed)
    manager.fail_by_id(claimed.task_id, "payload error")

    resp = client.get(f"/v1/joblib/tasks/{claimed.task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "failed"
    assert task.error == "payload error"
    assert task.completed_at is not None


def test_job_manager_cancel_from_claimed(api, client):
    """manager.cancel(task) should transition from CLAIMED to CANCELLED."""
    manager, claimed = _register_claim(api, client)

    manager.cancel(claimed)

    resp = client.get(f"/v1/joblib/tasks/{claimed.task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "cancelled"


def test_job_manager_cancel_from_running(api, client):
    """manager.cancel(task) should transition from RUNNING to CANCELLED."""
    manager, claimed = _register_claim(api, client)

    manager.start(claimed)
    manager.cancel(claimed)

    resp = client.get(f"/v1/joblib/tasks/{claimed.task_id}")
    task = TaskResponse.model_validate(resp.json())
    assert task.status.value == "cancelled"


def test_job_manager_context_manager(api, client):
    """JobManager should support with-statement for automatic disconnect."""
    mock_tsio = MagicMock()

    with JobManager(api, tsio=mock_tsio) as manager:

        @manager.register
        class CtxJob(ConcreteExtension):
            category: ClassVar[Category] = Category.MODIFIER

        worker_id = manager.worker_id
        assert worker_id is not None

    # After exiting context, state should be cleared
    assert manager.worker_id is None
    assert len(manager) == 0

    # Worker should be deleted from server
    resp = client.patch(f"/v1/joblib/workers/{worker_id}")
    assert resp.status_code == 404


# --- Provider request handler (ProviderRequest SIO dispatch) ---


def test_sio_handlers_registered_in_init(api, client):
    """SIO handlers (ProviderRequest, TaskAvailable) registered in __init__."""
    mock_tsio = MagicMock()
    JobManager(api, tsio=mock_tsio)

    # tsio.on should have been called for ProviderRequest (and TaskAvailable)
    on_calls = [c[0][0] for c in mock_tsio.on.call_args_list]
    assert ProviderRequest in on_calls


def test_sio_handlers_not_registered_without_tsio(api, client):
    """No error when tsio is None (no SIO handlers)."""
    manager = JobManager(api)
    assert manager.tsio is None  # Just works, no error


def test_provider_request_handler_dispatches_read(api, client, fs_provider):
    """ProviderRequest handler should call read() and POST result back."""
    from zndraw_joblib.dependencies import request_hash as compute_hash

    mock_tsio = MagicMock()
    manager = JobManager(api, tsio=mock_tsio)

    mock_handler = MagicMock()
    mock_handler.list_dir.return_value = [{"name": "file.xyz", "size": 42}]

    manager.register_provider(fs_provider, name="local", handler=mock_handler)

    params = {"path": "/data"}
    rhash = compute_hash(params)

    # Trigger a read to dispatch (times out immediately with wait=0)
    read_resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
        params=params,
        headers={"Prefer": "wait=0"},
    )
    assert read_resp.status_code == 504

    # Capture the registered handler callback (registered in __init__)
    pr_calls = [c for c in mock_tsio.on.call_args_list if c[0][0] is ProviderRequest]
    callback = pr_calls[0][0][1]

    # Simulate incoming ProviderRequest (as the server would dispatch)
    event = ProviderRequest.from_dict_params(
        request_id=rhash,
        provider_name="@global:filesystem:local",
        params=params,
    )
    callback(event)

    # Verify provider.read() was called with the handler
    mock_handler.list_dir.assert_called_once_with("/data")

    # Verify result is now cached on the server
    cached_resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
        params=params,
    )
    assert cached_resp.status_code == 200
    assert cached_resp.json() == [{"name": "file.xyz", "size": 42}]


def test_provider_request_handler_ignores_unknown_provider(api, client, fs_provider):
    """Handler should silently ignore ProviderRequest for unknown providers."""
    mock_tsio = MagicMock()
    manager = JobManager(api, tsio=mock_tsio)

    manager.register_provider(fs_provider, name="local", handler=object())

    # Capture callback (registered in __init__)
    pr_calls = [c for c in mock_tsio.on.call_args_list if c[0][0] is ProviderRequest]
    callback = pr_calls[0][0][1]

    # Invoke with unknown provider name — should not raise
    event = ProviderRequest.from_dict_params(
        request_id="xyz",
        provider_name="@global:filesystem:unknown",
        params={"path": "/"},
    )
    callback(event)  # Should not raise
