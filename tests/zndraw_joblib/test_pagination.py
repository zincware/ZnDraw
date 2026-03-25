# tests/test_pagination.py
"""Tests for pagination boundary conditions across list endpoints."""

import pytest

from zndraw_joblib.schemas import JobSummary, PaginatedResponse


@pytest.fixture
def five_job_client(client):
    """Client with 5 global jobs registered."""
    for i in range(5):
        resp = client.put(
            "/v1/joblib/rooms/@global/jobs",
            json={"category": "modifiers", "name": f"Job{i}", "schema": {}},
        )
        assert resp.status_code == 201
    return client


def test_pagination_default_returns_all_when_under_limit(five_job_client):
    """Default limit (50) returns all items when total < limit."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs")
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.total == 5
    assert len(page.items) == 5
    assert page.limit == 50
    assert page.offset == 0


def test_pagination_limit_restricts_items(five_job_client):
    """Setting limit=2 returns only 2 items but total stays at 5."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs?limit=2")
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.total == 5
    assert len(page.items) == 2
    assert page.limit == 2
    assert page.offset == 0


def test_pagination_offset_skips_items(five_job_client):
    """Setting offset=3 skips first 3 items, returns remaining 2."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs?offset=3")
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.total == 5
    assert len(page.items) == 2
    assert page.offset == 3


def test_pagination_limit_and_offset_combined(five_job_client):
    """limit=2&offset=1 returns items 2-3 of 5."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs?limit=2&offset=1")
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.total == 5
    assert len(page.items) == 2
    assert page.limit == 2
    assert page.offset == 1


def test_pagination_offset_beyond_total_returns_empty(five_job_client):
    """Offset past total returns empty items but correct total."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs?offset=100")
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.total == 5
    assert len(page.items) == 0
    assert page.offset == 100


def test_pagination_total_stable_across_pages(five_job_client):
    """Total count stays the same regardless of limit/offset."""
    totals = []
    for offset in range(0, 6, 2):
        response = five_job_client.get(
            f"/v1/joblib/rooms/@global/jobs?limit=2&offset={offset}"
        )
        page = PaginatedResponse[JobSummary].model_validate(response.json())
        totals.append(page.total)
    assert all(t == 5 for t in totals)


def test_pagination_all_items_reachable(five_job_client):
    """Walking through pages collects all items without duplicates."""
    all_names = []
    offset = 0
    while True:
        response = five_job_client.get(
            f"/v1/joblib/rooms/@global/jobs?limit=2&offset={offset}"
        )
        page = PaginatedResponse[JobSummary].model_validate(response.json())
        if not page.items:
            break
        all_names.extend(j.full_name for j in page.items)
        offset += 2

    assert len(all_names) == 5
    assert len(set(all_names)) == 5  # No duplicates


def test_pagination_limit_zero_returns_count_only(five_job_client):
    """limit=0 should return count only: items=[], total=5."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs?limit=0")
    assert response.status_code == 200
    page = PaginatedResponse[JobSummary].model_validate(response.json())
    assert page.items == []
    assert page.total == 5
    assert page.limit == 0


def test_pagination_negative_offset_rejected(five_job_client):
    """offset=-1 should be rejected (ge=0 validation)."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs?offset=-1")
    assert response.status_code == 422


def test_pagination_limit_above_max_rejected(five_job_client):
    """limit=501 should be rejected (le=500 validation)."""
    response = five_job_client.get("/v1/joblib/rooms/@global/jobs?limit=501")
    assert response.status_code == 422
