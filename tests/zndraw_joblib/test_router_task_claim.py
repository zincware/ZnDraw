# tests/test_router_task_claim.py
"""Tests for task claim endpoint using shared fixtures."""

from zndraw_joblib.schemas import TaskClaimResponse


def test_claim_task_returns_null_when_empty(seeded_client):
    worker_id = seeded_client.seeded_worker_id
    response = seeded_client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker_id},
    )
    assert response.status_code == 200
    data = TaskClaimResponse.model_validate(response.json())
    assert data.task is None


def test_claim_task_returns_oldest_first(seeded_client):
    worker_id = seeded_client.seeded_worker_id
    # Submit two tasks
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {"order": 1}},
    )
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {"order": 2}},
    )

    # Claim should return first (oldest)
    response = seeded_client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker_id},
    )
    assert response.status_code == 200
    data = TaskClaimResponse.model_validate(response.json())
    assert data.task is not None
    assert data.task.payload["order"] == 1
    assert data.task.status.value == "claimed"


def test_claim_task_marks_as_claimed(seeded_client):
    worker_id = seeded_client.seeded_worker_id
    seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )

    response = seeded_client.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker_id},
    )
    data = TaskClaimResponse.model_validate(response.json())
    assert data.task is not None
    assert data.task.status.value == "claimed"


def test_claim_task_only_registered_jobs(client_factory):
    """Worker can only claim tasks for jobs they are registered for."""
    client1 = client_factory("worker_1")
    client2 = client_factory("worker_2")

    # Worker 1 creates worker and registers job
    resp1 = client1.post("/v1/joblib/workers")
    worker1_id = resp1.json()["id"]
    client1.put(
        "/v1/joblib/rooms/@global/jobs",
        json={
            "category": "modifiers",
            "name": "Rotate",
            "schema": {},
            "worker_id": worker1_id,
        },
    )

    # Worker 2 creates worker (but doesn't register the job)
    resp2 = client2.post("/v1/joblib/workers")
    worker2_id = resp2.json()["id"]

    # Submit task
    client1.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )

    # Worker 2 tries to claim (not registered)
    response = client2.post(
        "/v1/joblib/tasks/claim",
        json={"worker_id": worker2_id},
    )
    data = TaskClaimResponse.model_validate(response.json())
    assert data.task is None  # Can't claim - not registered
