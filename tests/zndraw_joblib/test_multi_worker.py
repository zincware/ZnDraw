# tests/test_multi_worker.py
"""Tests for multi-worker scenarios - fixing the MultipleResultsFound bug."""

import uuid

from zndraw_joblib.schemas import TaskClaimResponse


def test_two_workers_same_user_can_claim_tasks(client):
    """
    Core bug reproduction: A user with two workers registered for the same job
    should be able to claim tasks without MultipleResultsFound error.
    """
    # Create first worker and register job
    resp1 = client.post("/v1/joblib/workers")
    assert resp1.status_code == 201
    worker1_id = resp1.json()["id"]

    resp = client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "SharedJob", "worker_id": worker1_id},
    )
    assert resp.status_code == 201

    # Create second worker and register same job
    resp2 = client.post("/v1/joblib/workers")
    assert resp2.status_code == 201
    worker2_id = resp2.json()["id"]

    resp = client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "SharedJob", "worker_id": worker2_id},
    )
    assert resp.status_code == 200  # Job already exists, worker added

    # Submit a task
    resp = client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:SharedJob",
        json={"payload": {"data": "test"}},
    )
    assert resp.status_code == 202

    # Worker 1 claims the task - should NOT raise MultipleResultsFound
    response = client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker1_id},
    )
    assert response.status_code == 200
    data = TaskClaimResponse.model_validate(response.json())
    assert data.task is not None
    assert str(data.task.worker_id) == worker1_id


def test_claim_requires_worker_id(seeded_client):
    """Claim endpoint should return 422 if worker_id is missing."""
    # Submit a task first
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )

    # Try to claim without worker_id - should fail validation
    response = seeded_client.post("/v1/joblib/tasks/claim", json={})
    assert response.status_code == 422  # Validation error


def test_claim_with_unknown_worker_id_returns_404(seeded_client):
    """Claim with unknown worker_id should return 404."""
    # Submit a task first
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )

    unknown_worker_id = str(uuid.uuid4())
    response = seeded_client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": unknown_worker_id},
    )
    assert response.status_code == 404


def test_claim_with_other_users_worker_returns_403(client_factory):
    """Claim with another user's worker_id should return 403."""
    client1 = client_factory("user_1")
    client2 = client_factory("user_2")

    # User 1 creates a worker
    resp1 = client1.post("/v1/joblib/workers")
    assert resp1.status_code == 201
    worker1_id = resp1.json()["id"]

    # User 1 registers job and submits task
    client1.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "PrivateJob", "worker_id": worker1_id},
    )
    client1.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:PrivateJob",
        json={"payload": {}},
    )

    # User 2 creates their own worker
    resp2 = client2.post("/v1/joblib/workers")
    assert resp2.status_code == 201

    # User 2 tries to claim using User 1's worker_id - should fail with 403
    response = client2.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker1_id},
    )
    assert response.status_code == 403


def test_multiple_workers_claim_different_tasks(client):
    """Two workers registered for same job can each claim different tasks."""
    # Create first worker and register job
    resp1 = client.post("/v1/joblib/workers")
    worker1_id = resp1.json()["id"]

    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "ParallelJob", "worker_id": worker1_id},
    )

    # Create second worker and register same job
    resp2 = client.post("/v1/joblib/workers")
    worker2_id = resp2.json()["id"]

    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "ParallelJob", "worker_id": worker2_id},
    )

    # Submit two tasks
    for i in range(2):
        client.post(
            "/v1/joblib/rooms/room_1/tasks/@global:modifiers:ParallelJob",
            json={"payload": {"index": i}},
        )

    # Worker 1 claims first task
    response1 = client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker1_id},
    )
    assert response1.status_code == 200
    data1 = TaskClaimResponse.model_validate(response1.json())
    assert data1.task is not None
    assert data1.task.payload["index"] == 0

    # Worker 2 claims second task
    response2 = client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker2_id},
    )
    assert response2.status_code == 200
    data2 = TaskClaimResponse.model_validate(response2.json())
    assert data2.task is not None
    assert data2.task.payload["index"] == 1


def test_two_workers_can_complete_tasks(client):
    """
    Bug reproduction: Two workers complete tasks for the same job.
    The orphan cleanup code should not fail with MultipleResultsFound.
    """
    # Create first worker and register job
    resp1 = client.post("/v1/joblib/workers")
    worker1_id = resp1.json()["id"]

    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={
            "category": "modifiers",
            "name": "CompletionJob",
            "worker_id": worker1_id,
        },
    )

    # Create second worker and register same job
    resp2 = client.post("/v1/joblib/workers")
    worker2_id = resp2.json()["id"]

    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={
            "category": "modifiers",
            "name": "CompletionJob",
            "worker_id": worker2_id,
        },
    )

    # Submit two tasks
    task_ids = []
    for i in range(2):
        resp = client.post(
            "/v1/joblib/rooms/room_1/tasks/@global:modifiers:CompletionJob",
            json={"payload": {"index": i}},
        )
        task_ids.append(resp.json()["id"])

    # Worker 1 claims and completes first task
    response1 = client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker1_id},
    )
    assert response1.status_code == 200
    claimed_task_1 = response1.json()["task"]
    assert claimed_task_1 is not None

    # Transition: claimed -> running -> completed
    running_resp = client.patch(
        f"/v1/joblib/tasks/{claimed_task_1['id']}",
        json={"status": "running"},
    )
    assert running_resp.status_code == 200

    # Complete the task - this triggers orphan cleanup check
    complete_resp = client.patch(
        f"/v1/joblib/tasks/{claimed_task_1['id']}",
        json={"status": "completed"},
    )
    assert complete_resp.status_code == 200, (
        f"Failed: {complete_resp.json()}"
    )  # Should NOT fail with MultipleResultsFound

    # Worker 2 claims and completes second task
    response2 = client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker2_id},
    )
    assert response2.status_code == 200
    claimed_task_2 = response2.json()["task"]
    assert claimed_task_2 is not None

    # Transition: claimed -> running -> completed
    running_resp2 = client.patch(
        f"/v1/joblib/tasks/{claimed_task_2['id']}",
        json={"status": "running"},
    )
    assert running_resp2.status_code == 200

    # Complete the second task
    complete_resp2 = client.patch(
        f"/v1/joblib/tasks/{claimed_task_2['id']}",
        json={"status": "completed"},
    )
    assert complete_resp2.status_code == 200, f"Failed: {complete_resp2.json()}"


def test_worker_can_only_claim_registered_jobs(client):
    """Worker can only claim tasks for jobs they are registered for."""
    # Create two workers
    resp1 = client.post("/v1/joblib/workers")
    worker1_id = resp1.json()["id"]

    resp2 = client.post("/v1/joblib/workers")
    worker2_id = resp2.json()["id"]

    # Worker 1 registers job A
    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "JobA", "worker_id": worker1_id},
    )

    # Worker 2 registers job B
    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "JobB", "worker_id": worker2_id},
    )

    # Submit task for job A
    client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:JobA",
        json={"payload": {}},
    )

    # Worker 2 tries to claim - should get nothing (not registered for JobA)
    response = client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker2_id},
    )
    assert response.status_code == 200
    data = TaskClaimResponse.model_validate(response.json())
    assert data.task is None

    # Worker 1 can claim it
    response = client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker1_id},
    )
    assert response.status_code == 200
    data = TaskClaimResponse.model_validate(response.json())
    assert data.task is not None
