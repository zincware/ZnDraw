"""Unit tests for job endpoints and Job object."""

import pytest
import requests

from zndraw import ZnDraw
from zndraw.extensions import Category, Extension
from zndraw.job import Job


class TestModifier(Extension):
    """Test modifier extension."""

    category = Category.MODIFIER
    param: int = 42

    def run(self, vis: ZnDraw, **kwargs):
        """Simple run method for testing."""
        pass


# =============================================================================
# Test Suite A: Job Submission
# =============================================================================


def test_submit_job_returns_job_object(server):
    """Test that vis.run() returns a Job object."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))

    assert isinstance(job, Job)
    assert job.job_id is not None
    assert isinstance(job.job_id, str)


def test_submit_job_creates_assigned_status(server):
    """Test that submitted job gets assigned when worker available (but not processed with auto_pickup=False)."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=100))

    # Job is assigned to worker immediately, but worker doesn't process it (auto_pickup=False)
    assert job.status == "assigned"
    assert job.is_assigned()
    assert not job.is_processing()
    assert not job.is_completed()
    assert not job.is_failed()


def test_submit_job_without_workers_stays_pending(server):
    """Test that second job stays pending when first job occupies the only worker."""
    vis = ZnDraw(
        url=server, room="test-pending", user="testuser", auto_pickup_jobs=False
    )
    vis.register_extension(TestModifier)

    # Submit two jobs - only one worker available
    job1 = vis.run(TestModifier(param=1))
    job2 = vis.run(TestModifier(param=2))

    # First job is assigned to the worker
    assert job1.status == "assigned"

    # Second job should be pending (worker is busy)
    assert job2.status == "pending"


# =============================================================================
# Test Suite B: GET /api/jobs/{job_id}
# =============================================================================


def test_get_job_details_success(server):
    """Test fetching job details via GET endpoint."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))

    # Fetch via REST API directly
    response = requests.get(f"{server}/api/jobs/{job.job_id}")

    assert response.status_code == 200
    details = response.json()
    assert details["jobId"] == job.job_id
    assert (
        details["status"] == "assigned"
    )  # Job assigned but not processed (auto_pickup=False)
    assert details["category"] == "modifiers"
    assert details["extension"] == "TestModifier"
    assert details["data"]["param"] == 42


def test_get_job_details_not_found(server):
    """Test GET with invalid job ID returns 404."""
    response = requests.get(f"{server}/api/jobs/nonexistent-job-id")

    assert response.status_code == 404


def test_get_job_details_contains_all_fields(server):
    """Test that job details contain all expected fields."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=123))

    details = job.refresh()

    # Verify all expected fields are present
    required_fields = [
        "jobId",
        "room",
        "category",
        "extension",
        "data",
        "public",
        "status",
    ]
    for field in required_fields:
        assert field in details, f"Missing required field: {field}"


# =============================================================================
# Test Suite C: PUT /api/jobs/{job_id}/status
# =============================================================================


def test_update_job_status_to_processing(server, get_jwt_auth_headers):
    """Test updating job status from assigned to processing."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))

    # Job is already assigned, transition to processing
    # Use job.worker_id (the assigned worker) instead of vis.sid
    response = requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={"status": "processing", "workerId": job.worker_id},
        headers=get_jwt_auth_headers(server, "testuser"),
    )

    assert response.status_code == 200
    job.refresh()
    assert job.is_processing()


def test_update_job_status_to_completed(server, get_jwt_auth_headers):
    """Test updating job status to completed."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))
    auth_headers = get_jwt_auth_headers(server, "testuser")
    worker_id = job.worker_id

    # Transition: assigned → processing → completed
    requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={"status": "processing", "workerId": worker_id},
        headers=auth_headers,
    )
    response = requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={"status": "completed", "workerId": worker_id},
        headers=auth_headers,
    )

    assert response.status_code == 200
    job.refresh()
    assert job.is_completed()
    assert job.is_done()


def test_update_job_status_to_failed(server, get_jwt_auth_headers):
    """Test updating job status to failed with error message."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))
    auth_headers = get_jwt_auth_headers(server, "testuser")
    worker_id = job.worker_id

    # Transition to processing first
    requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={"status": "processing", "workerId": worker_id},
        headers=auth_headers,
    )

    # Mark as failed
    response = requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={
            "status": "failed",
            "workerId": worker_id,
            "error": "Something went wrong",
        },
        headers=auth_headers,
    )

    assert response.status_code == 200
    job.refresh()
    assert job.is_failed()
    assert job.is_done()


# =============================================================================
# Test Suite D: Job Object Behavior
# =============================================================================


def test_job_refresh(server, get_jwt_auth_headers):
    """Test that job.refresh() updates status from server."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))

    initial_status = job.status
    assert initial_status == "assigned"  # Assigned immediately with worker available

    # Update via API to processing using the assigned worker ID
    requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={"status": "processing", "workerId": job.worker_id},
        headers=get_jwt_auth_headers(server, "testuser"),
    )

    # Refresh and check
    job.refresh()
    assert job.status != initial_status
    assert job.status == "processing"


def test_job_repr(server):
    """Test Job.__repr__() for terminal output."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))

    repr_str = repr(job)
    assert "Job(" in repr_str
    assert job.job_id in repr_str
    assert "assigned" in repr_str  # Job assigned (worker available)


def test_job_wait_for_completion_with_timeout(server):
    """Test job.wait() with timeout."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))

    # Job will stay pending (no worker), so wait should timeout
    with pytest.raises(TimeoutError):
        job.wait(timeout=0.5, poll_interval=0.1)


def test_job_completion_lifecycle(server, get_jwt_auth_headers):
    """Test complete job lifecycle from submission to completion."""
    vis = ZnDraw(url=server, room="test", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))
    auth_headers = get_jwt_auth_headers(server, "testuser")
    worker_id = job.worker_id

    # Job is assigned but not yet processing
    assert job.is_assigned()
    assert not job.is_done()

    # Transition to processing
    requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={"status": "processing", "workerId": worker_id},
        headers=auth_headers,
    )

    # Mark as completed (extensions work through side effects, no result needed)
    response = requests.put(
        f"{server}/api/rooms/test/jobs/{job.job_id}/status",
        json={"status": "completed", "workerId": worker_id},
        headers=auth_headers,
    )

    assert response.status_code == 200
    job.refresh()
    assert job.is_completed()
    assert job.is_done()


# =============================================================================
# Test Suite E: Integration with Workers
# =============================================================================


def test_multiple_jobs_queue_correctly(server):
    """Test that multiple jobs can be submitted and have unique IDs."""
    vis = ZnDraw(url=server, room="test-queue", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(TestModifier)

    # Submit multiple jobs
    job1 = vis.run(TestModifier(param=1))
    job2 = vis.run(TestModifier(param=2))
    job3 = vis.run(TestModifier(param=3))

    # First job is assigned, others are pending (only one worker)
    assert job1.is_assigned()  # Worker picks up first job
    assert job2.is_pending()  # Second job waits in queue
    assert job3.is_pending()  # Third job waits in queue

    # All should have different IDs
    assert job1.job_id != job2.job_id
    assert job2.job_id != job3.job_id
    assert job1.job_id != job3.job_id

    vis.disconnect()


# =============================================================================
# Test Suite F: ASSIGNED Job Timeout
# =============================================================================


def test_stale_assigned_job_is_failed_on_list(server, redis_client):
    """Test that jobs stuck in ASSIGNED state are failed during list_jobs."""
    from datetime import datetime, timezone

    from zndraw.app.redis_keys import JobKeys

    vis = ZnDraw(
        url=server, room="test-timeout", user="testuser", auto_pickup_jobs=False
    )
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))
    assert job.is_assigned()

    # Artificially age the job by setting created_at to be older than the timeout
    job_keys = JobKeys(job.job_id)
    old_time = (
        datetime.now(timezone.utc).replace(year=2020).isoformat()
    )  # Very old timestamp
    redis_client.hset(job_keys.hash_key(), "created_at", old_time)

    # List jobs - this triggers lazy cleanup
    response = requests.get(f"{server}/api/rooms/test-timeout/jobs")
    assert response.status_code == 200

    # Refresh the job and verify it's now failed
    job.refresh()
    assert job.is_failed()
    assert job.is_done()
    assert job.error is not None
    assert "timed out in ASSIGNED state" in job.error

    vis.disconnect()


def test_assigned_job_within_timeout_is_not_failed(server):
    """Test that jobs in ASSIGNED state within timeout are not affected."""
    vis = ZnDraw(
        url=server, room="test-no-timeout", user="testuser", auto_pickup_jobs=False
    )
    vis.register_extension(TestModifier)

    job = vis.run(TestModifier(param=42))
    assert job.is_assigned()

    # List jobs - this triggers cleanup, but job is fresh so should not be failed
    response = requests.get(f"{server}/api/rooms/test-no-timeout/jobs")
    assert response.status_code == 200

    # Job should still be assigned
    job.refresh()
    assert job.is_assigned()
    assert not job.is_failed()

    vis.disconnect()
