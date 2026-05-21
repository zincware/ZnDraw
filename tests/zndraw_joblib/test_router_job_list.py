# tests/test_router_job_list.py
"""Tests for job listing and details endpoints using shared fixtures."""

import pytest
from conftest import make_room_address

from zndraw_joblib.exceptions import ProblemDetail
from zndraw_joblib.schemas import JobResponse, JobSummary, PaginatedResponse


@pytest.fixture
def multi_job_client(client, test_user_id):
    """Client with multiple jobs registered."""
    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "selections", "name": "All", "schema": {}},
    )
    addr = make_room_address(test_user_id, "room_123")
    client.put(
        f"/v1/joblib/rooms/{addr}/jobs",
        json={"category": "modifiers", "name": "Translate", "schema": {}},
    )
    client.room_123_address = addr  # expose for consumers
    return client


def test_list_jobs_global_only(multi_job_client):
    response = multi_job_client.get("/v1/joblib/rooms/@global/jobs")
    assert response.status_code == 200
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.total == 2
    names = [j.full_name for j in page.items]
    assert "@global:modifiers:Rotate" in names
    assert "@global:selections:All" in names


def test_list_jobs_room_includes_global(multi_job_client):
    response = multi_job_client.get(
        f"/v1/joblib/rooms/{multi_job_client.room_123_address}/jobs"
    )
    assert response.status_code == 200
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.total == 3  # 2 global + 1 room
    names = [j.full_name for j in page.items]
    assert "@global:modifiers:Rotate" in names
    assert f"{multi_job_client.room_123_address}:modifiers:Translate" in names


def test_get_job_details(multi_job_client):
    response = multi_job_client.get(
        f"/v1/joblib/rooms/{multi_job_client.room_123_address}/jobs/@global:modifiers:Rotate"
    )
    assert response.status_code == 200
    data = JobResponse.model_validate(response.json())
    assert data.full_name == "@global:modifiers:Rotate"
    assert data.category == "modifiers"


def test_get_job_not_found(multi_job_client):
    response = multi_job_client.get(
        f"/v1/joblib/rooms/{multi_job_client.room_123_address}/jobs/@global:modifiers:NonExistent"
    )
    assert response.status_code == 404
    error = ProblemDetail.model_validate(response.json())
    assert error.status == 404
