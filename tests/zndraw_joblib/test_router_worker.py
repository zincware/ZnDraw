# tests/test_router_worker.py
"""Tests for worker heartbeat and deletion endpoints using shared fixtures."""

import time
from uuid import UUID

from zndraw_joblib.exceptions import ProblemDetail
from zndraw_joblib.schemas import (
    JobResponse,
    PaginatedResponse,
    TaskResponse,
    WorkerResponse,
    WorkerSummary,
)


def test_create_worker(client):
    """POST /workers should create a new worker for the authenticated user."""
    response = client.post("/v1/joblib/workers")
    assert response.status_code == 201
    data = WorkerResponse.model_validate(response.json())
    assert data.id is not None
    assert isinstance(data.id, UUID)
    assert data.last_heartbeat is not None


def test_worker_heartbeat_updates_timestamp(seeded_client):
    """Heartbeat should update timestamp, second call should have later timestamp."""
    worker_id = seeded_client.seeded_worker_id

    response1 = seeded_client.patch(f"/v1/joblib/workers/{worker_id}")
    assert response1.status_code == 200
    data1 = WorkerResponse.model_validate(response1.json())
    assert str(data1.id) == worker_id
    assert data1.last_heartbeat is not None

    time.sleep(0.01)  # Small delay to ensure timestamp differs

    response2 = seeded_client.patch(f"/v1/joblib/workers/{worker_id}")
    assert response2.status_code == 200
    data2 = WorkerResponse.model_validate(response2.json())

    assert data2.last_heartbeat > data1.last_heartbeat


def test_worker_heartbeat_not_found(client):
    response = client.patch("/v1/joblib/workers/00000000-0000-0000-0000-000000000000")
    assert response.status_code == 404
    error = ProblemDetail.model_validate(response.json())
    assert error.status == 404


def test_worker_delete(seeded_client):
    worker_id = seeded_client.seeded_worker_id
    response = seeded_client.delete(f"/v1/joblib/workers/{worker_id}")
    assert response.status_code == 204


def test_worker_delete_not_found(client):
    response = client.delete("/v1/joblib/workers/00000000-0000-0000-0000-000000000000")
    assert response.status_code == 404
    error = ProblemDetail.model_validate(response.json())
    assert error.status == 404


def test_worker_delete_removes_orphan_job(seeded_client):
    """Deleting sole worker of a job should remove the job too."""
    worker_id = seeded_client.seeded_worker_id

    seeded_client.delete(f"/v1/joblib/workers/{worker_id}")
    response = seeded_client.get(
        "/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate"
    )
    assert response.status_code == 404


def test_worker_delete_keeps_job_with_pending_task(seeded_client):
    """Job should remain if there are non-terminal tasks, even without workers."""
    worker_id = seeded_client.seeded_worker_id

    # Submit a task (creates a pending task)
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )

    seeded_client.delete(f"/v1/joblib/workers/{worker_id}")

    # Job should still exist because there's a pending task
    response = seeded_client.get(
        "/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate"
    )
    assert response.status_code == 200
    data = JobResponse.model_validate(response.json())
    assert data.full_name == "@global:modifiers:Rotate"
    assert data.workers == []


def test_worker_delete_removes_job_after_task_completes(seeded_client):
    """Job should be removed when sole worker deleted and all tasks are terminal."""
    worker_id = seeded_client.seeded_worker_id

    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    claim_resp = seeded_client.post(
        "/v1/joblib/tasks/claim", json={"worker_id": seeded_client.seeded_worker_id}
    )
    task_id = claim_resp.json()["task"]["id"]
    seeded_client.patch(f"/v1/joblib/tasks/{task_id}", json={"status": "running"})
    seeded_client.patch(f"/v1/joblib/tasks/{task_id}", json={"status": "completed"})

    seeded_client.delete(f"/v1/joblib/workers/{worker_id}")

    # Job should be removed (no workers, no pending tasks)
    response = seeded_client.get(
        "/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate"
    )
    assert response.status_code == 404


def test_workers_list_changes_with_workers(client_factory):
    """workers list should update when workers register and are removed."""
    client1 = client_factory("worker_1")
    client2 = client_factory("worker_2")
    client3 = client_factory("worker_3")

    # Worker 1 registers
    resp = client1.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code == 201
    data = JobResponse.model_validate(resp.json())
    assert len(data.workers) == 1
    worker1_id = data.worker_id

    # Worker 2 registers same job
    resp = client2.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code == 200  # Idempotent
    data = JobResponse.model_validate(resp.json())
    assert len(data.workers) == 2
    worker2_id = data.worker_id

    # Worker 3 registers same job
    resp = client3.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    data = JobResponse.model_validate(resp.json())
    assert len(data.workers) == 3
    worker3_id = data.worker_id

    # Remove worker 2
    client2.delete(f"/v1/joblib/workers/{worker2_id}")

    # Check workers list now has 2 workers
    resp = client1.get("/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate")
    data = JobResponse.model_validate(resp.json())
    assert len(data.workers) == 2

    # Remove worker 1
    client1.delete(f"/v1/joblib/workers/{worker1_id}")

    # Check workers list now has 1 worker
    resp = client3.get("/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate")
    data = JobResponse.model_validate(resp.json())
    assert len(data.workers) == 1
    assert data.workers[0] == worker3_id


def test_worker_delete_fails_running_tasks(seeded_client):
    """Deleting worker should mark their running/claimed tasks as failed."""
    worker_id = seeded_client.seeded_worker_id

    # Submit a task
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )

    # Claim and start running
    claim_resp = seeded_client.post(
        "/v1/joblib/tasks/claim", json={"worker_id": seeded_client.seeded_worker_id}
    )
    task_id = claim_resp.json()["task"]["id"]
    seeded_client.patch(f"/v1/joblib/tasks/{task_id}", json={"status": "running"})

    # Delete worker - task should be marked as failed
    seeded_client.delete(f"/v1/joblib/workers/{worker_id}")

    # Task should still exist and be marked as failed
    task_resp = seeded_client.get(f"/v1/joblib/tasks/{task_id}")
    assert task_resp.status_code == 200
    task_data = TaskResponse.model_validate(task_resp.json())
    assert task_data.status.value == "failed"
    assert "worker disconnected" in task_data.error.lower()

    # Job should be soft-deleted (returns 404 on GET)
    resp = seeded_client.get("/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate")
    assert resp.status_code == 404


def test_worker_delete_soft_deletes_job_but_keeps_task(seeded_client):
    """Job should be soft-deleted but task remains accessible."""
    worker_id = seeded_client.seeded_worker_id

    # Submit and complete a task
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {"test": "data"}},
    )
    claim_resp = seeded_client.post(
        "/v1/joblib/tasks/claim", json={"worker_id": seeded_client.seeded_worker_id}
    )
    task_id = claim_resp.json()["task"]["id"]
    seeded_client.patch(f"/v1/joblib/tasks/{task_id}", json={"status": "running"})
    seeded_client.patch(f"/v1/joblib/tasks/{task_id}", json={"status": "completed"})

    # Delete sole worker - job becomes orphan
    seeded_client.delete(f"/v1/joblib/workers/{worker_id}")

    # Job should return 404 (soft-deleted)
    resp = seeded_client.get("/v1/joblib/rooms/@global/jobs/@global:modifiers:Rotate")
    assert resp.status_code == 404

    # But task should still be accessible with all its data
    task_resp = seeded_client.get(f"/v1/joblib/tasks/{task_id}")
    assert task_resp.status_code == 200
    task_data = TaskResponse.model_validate(task_resp.json())
    assert task_data.status.value == "completed"
    assert task_data.payload == {"test": "data"}
    # Job name should still be available from the soft-deleted job
    assert task_data.job_name == "@global:modifiers:Rotate"


def test_list_workers_empty(client):
    """List workers returns empty list when no workers exist."""
    response = client.get("/v1/joblib/workers")
    assert response.status_code == 200
    page = PaginatedResponse[WorkerSummary].model_validate(response.json())
    assert page.items == []
    assert page.total == 0


def test_list_workers_returns_all(client_factory):
    """List workers returns all workers with job counts."""
    client1 = client_factory("worker-a")
    client2 = client_factory("worker-b")

    # Worker A registers two jobs
    client1.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "job1", "schema": {}},
    )
    client1.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "job2", "schema": {}},
    )

    # Worker B registers one job
    client2.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "job3", "schema": {}},
    )

    response = client1.get("/v1/joblib/workers")
    assert response.status_code == 200
    page = PaginatedResponse[WorkerSummary].model_validate(response.json())
    # Each registration without worker_id creates a new worker:
    # client1: 2 workers (job1 + job2), client2: 1 worker (job3)
    assert len(page.items) == 3


def test_list_workers_for_room_empty(client):
    """List workers for room returns empty list when no workers."""
    response = client.get("/v1/joblib/rooms/my-room/workers")
    assert response.status_code == 200
    page = PaginatedResponse[WorkerSummary].model_validate(response.json())
    assert page.items == []
    assert page.total == 0


def test_list_workers_for_room_filters_by_room(client_factory):
    """List workers for room only returns workers serving that room."""
    client_a = client_factory("worker-a")
    client_b = client_factory("worker-b")
    client_c = client_factory("worker-c")

    # Worker A serves room1 and @global
    client_a.put(
        "/v1/joblib/rooms/room1/jobs",
        json={"category": "modifiers", "name": "job1", "schema": {}},
    )
    resp = client_a.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "global-job", "schema": {}},
    )
    _ = resp.json()["worker_id"]  # Validate response has worker_id

    # Worker B serves room2
    client_b.put(
        "/v1/joblib/rooms/room2/jobs",
        json={"category": "modifiers", "name": "job2", "schema": {}},
    )

    # Worker C serves room1
    resp = client_c.put(
        "/v1/joblib/rooms/room1/jobs",
        json={"category": "selections", "name": "job3", "schema": {}},
    )
    _ = resp.json()["worker_id"]  # Validate response has worker_id

    # List workers for room1 - should include workers from A and C
    response = client_a.get("/v1/joblib/rooms/room1/workers")
    assert response.status_code == 200
    page = PaginatedResponse[WorkerSummary].model_validate(response.json())

    # Workers from room1 jobs (worker_a1, worker_c1) and @global jobs (worker_a2)
    assert len(page.items) == 3


def test_list_workers_for_global_room(client_factory):
    """List workers for @global room only returns workers serving @global jobs."""
    client_a = client_factory("worker-a")
    client_b = client_factory("worker-b")

    resp = client_a.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "global-job", "schema": {}},
    )
    worker_a_id = resp.json()["worker_id"]

    client_b.put(
        "/v1/joblib/rooms/room1/jobs",
        json={"category": "modifiers", "name": "job1", "schema": {}},
    )

    response = client_a.get("/v1/joblib/rooms/@global/workers")
    assert response.status_code == 200
    page = PaginatedResponse[WorkerSummary].model_validate(response.json())
    worker_ids = {str(w.id) for w in page.items}

    assert worker_a_id in worker_ids
    # worker_b shouldn't be here as it only serves room1


def test_worker_heartbeat_forbidden_for_non_owner(client_factory):
    """Heartbeat should return 403 when called by a different user."""
    client_owner = client_factory("owner")
    client_other = client_factory("other")

    # Owner creates a worker
    resp = client_owner.post("/v1/joblib/workers")
    assert resp.status_code == 201
    worker_id = resp.json()["id"]

    # Other user tries to heartbeat
    response = client_other.patch(f"/v1/joblib/workers/{worker_id}")
    assert response.status_code == 403


def test_worker_delete_forbidden_for_non_owner(client_factory):
    """Delete should return 403 when called by a non-owner, non-superuser."""
    client_owner = client_factory("owner")
    client_other = client_factory("other", is_superuser=False)

    # Owner creates a worker
    resp = client_owner.post("/v1/joblib/workers")
    assert resp.status_code == 201
    worker_id = resp.json()["id"]

    # Non-superuser, non-owner tries to delete
    response = client_other.delete(f"/v1/joblib/workers/{worker_id}")
    assert response.status_code == 403


def test_worker_delete_allowed_for_superuser(client_factory):
    """Superuser should be able to delete any worker."""
    client_owner = client_factory("owner")
    client_admin = client_factory("admin", is_superuser=True)

    # Owner creates a worker
    resp = client_owner.post("/v1/joblib/workers")
    assert resp.status_code == 201
    worker_id = resp.json()["id"]

    # Superuser can delete it
    response = client_admin.delete(f"/v1/joblib/workers/{worker_id}")
    assert response.status_code == 204
