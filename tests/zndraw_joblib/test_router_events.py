# tests/test_router_events.py
"""Tests that router endpoints emit Socket.IO events."""

from unittest.mock import AsyncMock

import pytest
from fastapi.testclient import TestClient

from zndraw_joblib.events import JobsInvalidate, TaskAvailable, TaskStatusEvent


@pytest.fixture
def mock_tsio():
    """Create a mock tsio wrapper that records emit calls."""
    tsio = AsyncMock()
    tsio.emit = AsyncMock()
    return tsio


@pytest.fixture
def app_with_tsio(app, mock_tsio):
    """App with tsio set on app.state."""
    app.state.tsio = mock_tsio
    return app


@pytest.fixture
def client_with_tsio(app_with_tsio):
    return TestClient(app_with_tsio)


def test_register_job_emits_jobs_invalidate(client_with_tsio, mock_tsio):
    """PUT /rooms/{room_id}/jobs should emit JobsInvalidate."""
    resp = client_with_tsio.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code in (200, 201)

    calls = mock_tsio.emit.call_args_list
    invalidate_calls = [c for c in calls if isinstance(c[0][0], JobsInvalidate)]
    assert len(invalidate_calls) >= 1
    assert invalidate_calls[0].kwargs["room"] == "room:@global"


def test_submit_task_emits_task_available(client_with_tsio, mock_tsio):
    """POST submit should emit TaskAvailable and TaskStatusEvent."""
    client_with_tsio.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    mock_tsio.emit.reset_mock()

    resp = client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    assert resp.status_code == 202

    calls = mock_tsio.emit.call_args_list
    available_calls = [c for c in calls if isinstance(c[0][0], TaskAvailable)]
    status_calls = [c for c in calls if isinstance(c[0][0], TaskStatusEvent)]

    assert len(available_calls) == 1
    assert available_calls[0].kwargs["room"] == "jobs:@global:modifiers:Rotate"
    assert len(status_calls) == 1
    assert status_calls[0].kwargs["room"] == "room:test-room"


def test_submit_internal_task_no_task_available(client_with_tsio, mock_tsio):
    """@internal task submission should NOT emit TaskAvailable."""
    from zndraw_joblib.registry import InternalRegistry

    client_with_tsio.put(
        "/v1/joblib/rooms/@internal/jobs",
        json={"category": "modifiers", "name": "InternalOp", "schema": {}},
    )

    mock_task = AsyncMock()
    registry = InternalRegistry(tasks={"@internal:modifiers:InternalOp": mock_task})
    client_with_tsio.app.state.internal_registry = registry

    mock_tsio.emit.reset_mock()

    resp = client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@internal:modifiers:InternalOp",
        json={"payload": {}},
    )
    assert resp.status_code == 202

    calls = mock_tsio.emit.call_args_list
    available_calls = [c for c in calls if isinstance(c[0][0], TaskAvailable)]
    assert len(available_calls) == 0


def test_submit_internal_task_emits_task_status_to_submitting_room(
    client_with_tsio, mock_tsio
):
    """@internal task submission should emit TaskStatusEvent to the submitting room."""
    from zndraw_joblib.registry import InternalRegistry

    client_with_tsio.put(
        "/v1/joblib/rooms/@internal/jobs",
        json={"category": "modifiers", "name": "InternalOp", "schema": {}},
    )

    mock_task = AsyncMock()
    registry = InternalRegistry(tasks={"@internal:modifiers:InternalOp": mock_task})
    client_with_tsio.app.state.internal_registry = registry

    mock_tsio.emit.reset_mock()

    resp = client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@internal:modifiers:InternalOp",
        json={"payload": {}},
    )
    assert resp.status_code == 202

    calls = mock_tsio.emit.call_args_list
    status_calls = [c for c in calls if isinstance(c[0][0], TaskStatusEvent)]
    assert len(status_calls) == 1
    assert status_calls[0].args[0].status == "claimed"
    assert status_calls[0].kwargs["room"] == "room:test-room"


def test_claim_task_emits_task_status(client_with_tsio, mock_tsio):
    """POST /tasks/claim should emit TaskStatusEvent with status=claimed."""
    resp = client_with_tsio.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    worker_id = resp.json()["worker_id"]

    client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    mock_tsio.emit.reset_mock()

    resp = client_with_tsio.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker_id},
    )
    assert resp.status_code == 200
    assert resp.json()["task"] is not None

    calls = mock_tsio.emit.call_args_list
    status_calls = [c for c in calls if isinstance(c[0][0], TaskStatusEvent)]
    assert len(status_calls) == 1
    assert status_calls[0].args[0].status == "claimed"
    assert status_calls[0].kwargs["room"] == "room:test-room"


def test_update_task_emits_task_status(client_with_tsio, mock_tsio):
    """PATCH /tasks/{id} should emit TaskStatusEvent."""
    resp = client_with_tsio.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    worker_id = resp.json()["worker_id"]

    resp = client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    task_id = resp.json()["id"]

    client_with_tsio.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker_id},
    )
    mock_tsio.emit.reset_mock()

    resp = client_with_tsio.patch(
        f"/v1/joblib/tasks/{task_id}",
        json={"status": "running"},
    )
    assert resp.status_code == 200

    calls = mock_tsio.emit.call_args_list
    status_calls = [c for c in calls if isinstance(c[0][0], TaskStatusEvent)]
    assert len(status_calls) == 1
    assert status_calls[0].args[0].status == "running"


def test_delete_worker_emits_events(client_with_tsio, mock_tsio):
    """DELETE /workers/{id} should emit JobsInvalidate + TaskStatusEvent."""
    resp = client_with_tsio.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    worker_id = resp.json()["worker_id"]

    client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    client_with_tsio.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker_id},
    )
    mock_tsio.emit.reset_mock()

    resp = client_with_tsio.delete(f"/v1/joblib/workers/{worker_id}")
    assert resp.status_code == 204

    calls = mock_tsio.emit.call_args_list
    invalidate_calls = [c for c in calls if isinstance(c[0][0], JobsInvalidate)]
    status_calls = [c for c in calls if isinstance(c[0][0], TaskStatusEvent)]
    assert len(invalidate_calls) >= 1
    assert len(status_calls) >= 1
    assert status_calls[0].args[0].status == "failed"


def test_delete_worker_emits_jobs_invalidate_when_job_has_other_workers(
    client_factory, mock_tsio
):
    """DELETE /workers/{id} should emit JobsInvalidate even when other workers remain."""
    c1 = client_factory("worker-1")
    c1.app.state.tsio = mock_tsio
    c2 = client_factory("worker-2")
    c2.app.state.tsio = mock_tsio

    # Worker 1 registers a job
    resp = c1.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "SharedJob", "schema": {}},
    )
    assert resp.status_code == 201
    worker1_id = resp.json()["worker_id"]

    # Worker 2 registers the same job
    resp = c2.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "SharedJob", "schema": {}},
    )
    assert resp.status_code == 200
    mock_tsio.emit.reset_mock()

    # Delete worker 1 — job still has worker 2
    resp = c1.delete(f"/v1/joblib/workers/{worker1_id}")
    assert resp.status_code == 204

    calls = mock_tsio.emit.call_args_list
    invalidate_calls = [c for c in calls if isinstance(c[0][0], JobsInvalidate)]
    assert len(invalidate_calls) >= 1, (
        "JobsInvalidate should be emitted when worker count changes"
    )
    assert invalidate_calls[0].kwargs["room"] == "room:@global"


def test_cancel_last_pending_task_emits_orphan_invalidate(client_with_tsio, mock_tsio):
    """Cancelling the last pending task on a workerless job should emit JobsInvalidate."""
    # Register a @global job (creates a worker)
    resp = client_with_tsio.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Orphan", "schema": {}},
    )
    worker_id = resp.json()["worker_id"]

    # Submit 2 pending tasks
    task1 = client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@global:modifiers:Orphan",
        json={"payload": {}},
    )
    task2 = client_with_tsio.post(
        "/v1/joblib/rooms/test-room/tasks/@global:modifiers:Orphan",
        json={"payload": {}},
    )
    task1_id = task1.json()["id"]
    task2_id = task2.json()["id"]

    # Delete the worker — tasks are PENDING so they stay PENDING.
    # Job has no workers but 2 pending tasks → not orphaned yet.
    client_with_tsio.delete(f"/v1/joblib/workers/{worker_id}")

    # Cancel task 1 — still 1 pending task left → no orphan
    client_with_tsio.patch(
        f"/v1/joblib/tasks/{task1_id}",
        json={"status": "cancelled"},
    )

    mock_tsio.emit.reset_mock()

    # Cancel task 2 — last pending task gone, no workers → job orphaned
    resp = client_with_tsio.patch(
        f"/v1/joblib/tasks/{task2_id}",
        json={"status": "cancelled"},
    )
    assert resp.status_code == 200

    calls = mock_tsio.emit.call_args_list
    invalidate_calls = [c for c in calls if isinstance(c[0][0], JobsInvalidate)]
    assert len(invalidate_calls) >= 1, (
        "JobsInvalidate should be emitted when orphan job is soft-deleted"
    )
    assert invalidate_calls[0].kwargs["room"] == "room:@global"


def test_no_tsio_does_not_break(client):
    """Endpoints work fine without tsio (default None)."""
    resp = client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "NoTsio", "schema": {}},
    )
    assert resp.status_code == 201
