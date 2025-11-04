"""Tests for error handling in the job queue system."""

import json

import pytest
import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw
from zndraw.extensions import Extension, ExtensionType


class TestExtension(Extension):
    category = ExtensionType.MODIFIER
    parameter: int

    def run(self, vis: ZnDraw, **kwargs):
        pass


def test_get_next_job_with_malformed_json_in_queue(server):
    """Test that malformed JSON in queue doesn't crash the server."""
    import redis

    room = "testroom"
    user = "testuser"

    # Create a worker
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(TestExtension)

    # Directly inject malformed JSON into the queue
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    queue_key = f"room:{room}:extensions:modifiers:{TestExtension.__name__}:queue"
    r.rpush(queue_key, "this is not valid JSON {{{")

    # Try to get next job - should handle the error gracefully
    response = requests.post(
        f"{server}/api/jobs/next", json={"workerId": vis.sid}
    )

    # Should return 400 (no jobs available) after skipping the malformed job
    assert response.status_code == 400
    assert response.json() == {"error": "No jobs available"}

    # Verify the malformed job was removed from the queue
    queue_length = r.llen(queue_key)
    assert queue_length == 0

    vis.disconnect()


def test_get_next_job_with_missing_job_id_in_queue(server):
    """Test that missing jobId in queue item doesn't crash the server."""
    import redis

    room = "testroom"
    user = "testuser"

    # Create a worker
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(TestExtension)

    # Create a valid job but with missing jobId
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    queue_key = f"room:{room}:extensions:modifiers:{TestExtension.__name__}:queue"

    # Valid JSON but missing jobId field
    malformed_job = json.dumps({"provider": "client", "data": {"parameter": 42}})
    r.rpush(queue_key, malformed_job)

    # Try to get next job
    response = requests.post(
        f"{server}/api/jobs/next", json={"workerId": vis.sid}
    )

    # Should return 400 (no jobs available) since the job is invalid
    assert response.status_code == 400
    assert response.json() == {"error": "No jobs available"}

    vis.disconnect()


def test_get_next_job_with_nonexistent_job_id(server):
    """Test that nonexistent job ID in queue doesn't crash the server."""
    import redis

    room = "testroom"
    user = "testuser"

    # Create a worker
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(TestExtension)

    # Create a queue item with a job ID that doesn't exist
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)
    queue_key = f"room:{room}:extensions:modifiers:{TestExtension.__name__}:queue"

    fake_job = json.dumps({"jobId": "nonexistent-job-id", "provider": "client"})
    r.rpush(queue_key, fake_job)

    # Try to get next job
    response = requests.post(
        f"{server}/api/jobs/next", json={"workerId": vis.sid}
    )

    # Should return 400 (no jobs available) since the job doesn't exist
    assert response.status_code == 400
    assert response.json() == {"error": "No jobs available"}

    # Verify the queue item was consumed
    queue_length = r.llen(queue_key)
    assert queue_length == 0

    vis.disconnect()


def test_get_next_job_without_worker_id(server):
    """Test that missing workerId returns proper error."""
    room = "testroom"

    response = requests.post(f"{server}/api/jobs/next", json={})

    assert response.status_code == 400
    assert response.json() == {"error": "workerId is required"}


def test_get_next_job_with_running_job(server):
    """Test that worker with running job cannot pick up another job."""
    room = "testroom"
    user = "testuser"

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(TestExtension)

    # Get authentication headers
    auth_headers = get_jwt_auth_headers(server, user)

    # Submit two jobs
    for i in range(2):
        response = requests.post(
            f"{server}/api/rooms/{room}/extensions/modifiers/{TestExtension.__name__}/submit",
            json={"data": {"parameter": i}, "userId": user},
            headers=auth_headers,
        )
        assert response.status_code == 200

    # Pick up first job
    response = requests.post(
        f"{server}/api/jobs/next", json={"workerId": vis.sid}
    )
    assert response.status_code == 200
    job_id = response.json()["jobId"]

    # Try to pick up second job while first is still running
    response = requests.post(
        f"{server}/api/jobs/next", json={"workerId": vis.sid}
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Worker is not idle"}

    # Complete first job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "completed", "workerId": vis.sid},
    )
    assert response.status_code == 200

    # Now should be able to pick up second job
    response = requests.post(
        f"{server}/api/jobs/next", json={"workerId": vis.sid}
    )
    assert response.status_code == 200

    vis.disconnect()


def test_extension_execution_with_auto_pickup(server):
    """Test that extensions execute correctly with auto_pickup_jobs=True."""
    room = "testroom"
    user = "testuser"

    # Create a worker with auto_pickup enabled (default)
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=True)
    vis.register_extension(TestExtension)

    # Get authentication headers
    auth_headers = get_jwt_auth_headers(server, user)

    # Submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{TestExtension.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=auth_headers,
    )
    assert response.status_code == 200
    job_id = response.json()["jobId"]

    # Wait for the job to be picked up and completed
    vis.socket.sio.sleep(1)

    # Check job status
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{job_id}")
    assert response.status_code == 200
    assert response.json()["status"] == "completed"

    vis.disconnect()
