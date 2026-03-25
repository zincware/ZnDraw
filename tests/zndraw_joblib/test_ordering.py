# tests/test_ordering.py
"""Tests that all paginated list endpoints return items newest-first (created_at DESC)."""

import time

import pytest

from zndraw_joblib.schemas import (
    JobSummary,
    PaginatedResponse,
    TaskResponse,
    WorkerSummary,
)

# ── Jobs ────────────────────────────────────────────────────────────────


@pytest.fixture
def ordered_job_client(client):
    """Client with 3 global jobs registered in known order."""
    for name in ["First", "Second", "Third"]:
        resp = client.put(
            "/v1/joblib/rooms/@global/jobs",
            json={"category": "modifiers", "name": name, "schema": {}},
        )
        assert resp.status_code == 201
        time.sleep(0.01)
    return client


def test_list_jobs_newest_first(ordered_job_client):
    """GET /rooms/{room_id}/jobs returns newest job first."""
    response = ordered_job_client.get("/v1/joblib/rooms/@global/jobs")
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    names = [j.name for j in page.items]
    assert names == ["Third", "Second", "First"]


def test_list_jobs_pagination_preserves_order(ordered_job_client):
    """Paginating through jobs maintains newest-first order."""
    # Page 1: limit=2 → Third, Second
    resp1 = ordered_job_client.get("/v1/joblib/rooms/@global/jobs?limit=2&offset=0")
    page1 = PaginatedResponse[JobSummary].model_validate(resp1.json())
    # Page 2: limit=2, offset=2 → First
    resp2 = ordered_job_client.get("/v1/joblib/rooms/@global/jobs?limit=2&offset=2")
    page2 = PaginatedResponse[JobSummary].model_validate(resp2.json())

    all_names = [j.name for j in page1.items] + [j.name for j in page2.items]
    assert all_names == ["Third", "Second", "First"]


# ── Tasks for room ─────────────────────────────────────────────────────


@pytest.fixture
def ordered_task_client(seeded_client):
    """Seeded client with 3 tasks submitted in known order."""
    for i in range(3):
        resp = seeded_client.post(
            "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
            json={"payload": {"index": i}},
        )
        assert resp.status_code == 202
        time.sleep(0.01)
    return seeded_client


def test_list_tasks_for_room_newest_first(ordered_task_client):
    """GET /rooms/{room_id}/tasks returns newest task first."""
    response = ordered_task_client.get("/v1/joblib/rooms/room_1/tasks")
    page = PaginatedResponse[TaskResponse].model_validate(response.json())
    indices = [t.payload["index"] for t in page.items]
    assert indices == [2, 1, 0]


def test_list_tasks_for_room_pagination_preserves_order(ordered_task_client):
    """Paginating through tasks maintains newest-first order."""
    resp1 = ordered_task_client.get("/v1/joblib/rooms/room_1/tasks?limit=2&offset=0")
    page1 = PaginatedResponse[TaskResponse].model_validate(resp1.json())
    resp2 = ordered_task_client.get("/v1/joblib/rooms/room_1/tasks?limit=2&offset=2")
    page2 = PaginatedResponse[TaskResponse].model_validate(resp2.json())

    indices = [t.payload["index"] for t in page1.items] + [
        t.payload["index"] for t in page2.items
    ]
    assert indices == [2, 1, 0]


# ── Tasks for job ──────────────────────────────────────────────────────


def test_list_tasks_for_job_newest_first(ordered_task_client):
    """GET /rooms/{room_id}/jobs/{job}/tasks returns newest task first."""
    response = ordered_task_client.get(
        "/v1/joblib/rooms/room_1/jobs/@global:modifiers:Rotate/tasks"
    )
    page = PaginatedResponse[TaskResponse].model_validate(response.json())
    indices = [t.payload["index"] for t in page.items]
    assert indices == [2, 1, 0]


# ── Workers (global list) ──────────────────────────────────────────────


def test_list_workers_newest_first(client_factory):
    """GET /workers returns newest worker first."""
    worker_ids = []
    for name in ["worker-a", "worker-b", "worker-c"]:
        c = client_factory(name)
        resp = c.put(
            "/v1/joblib/rooms/@global/jobs",
            json={"category": "modifiers", "name": f"job-{name}", "schema": {}},
        )
        assert resp.status_code == 201
        worker_ids.append(resp.json()["worker_id"])
        time.sleep(0.01)

    response = c.get("/v1/joblib/workers")
    page = PaginatedResponse[WorkerSummary].model_validate(response.json())
    returned_ids = [str(w.id) for w in page.items]
    # Newest (worker-c) should be first
    assert returned_ids[0] == worker_ids[2]
    assert returned_ids[-1] == worker_ids[0]


# ── Workers for room ───────────────────────────────────────────────────


def test_list_workers_for_room_newest_first(client_factory):
    """GET /rooms/{room_id}/workers returns newest worker first."""
    worker_ids = []
    for name in ["worker-x", "worker-y", "worker-z"]:
        c = client_factory(name)
        resp = c.put(
            "/v1/joblib/rooms/room_1/jobs",
            json={"category": "modifiers", "name": f"job-{name}", "schema": {}},
        )
        assert resp.status_code == 201
        worker_ids.append(resp.json()["worker_id"])
        time.sleep(0.01)

    response = c.get("/v1/joblib/rooms/room_1/workers")
    page = PaginatedResponse[WorkerSummary].model_validate(response.json())
    returned_ids = [str(w.id) for w in page.items]
    # Newest (worker-z) should be first
    assert returned_ids[0] == worker_ids[2]
    assert returned_ids[-1] == worker_ids[0]
