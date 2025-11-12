"""Tests for job execution with auto-pickup.

Note: Tests for the old /api/jobs/next polling endpoint have been removed
since the system now uses Socket.IO push-based job assignment.
"""

import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw
from zndraw.extensions import Extension, Category


class TestExtension(Extension):
    category = Category.MODIFIER
    parameter: int

    def run(self, vis: ZnDraw, **kwargs):
        pass


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
        f"{server}/api/rooms/{room}/extensions/private/modifiers/{TestExtension.__name__}/submit",
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
