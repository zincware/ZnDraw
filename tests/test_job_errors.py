"""Tests for job execution with auto-pickup.

Note: Tests for the old /api/jobs/next polling endpoint have been removed
since the system now uses Socket.IO push-based job assignment.
"""

import requests

from zndraw import ZnDraw
from zndraw.extensions import Category, Extension


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
    vis.socket.sio.sleep(0.5)  # Give registration time to complete

    # Create a second client to submit the job
    client = ZnDraw(url=server, room=room, user="submitter")

    # Submit and wait for job using the Job.wait() pattern
    job = client.run(TestExtension(parameter=42))
    job.wait(timeout=30)

    # Check job completed
    assert job.status == "completed"

    client.disconnect()
    vis.disconnect()
