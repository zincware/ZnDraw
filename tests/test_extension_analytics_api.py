"""Tests for extension analytics API endpoints."""

import pytest
import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw
from zndraw.app.job_manager import JobManager
from zndraw.extensions import Extension, ExtensionType


class DemoModifier(Extension):
    category = ExtensionType.MODIFIER
    parameter: int = 0

    def run(self, vis: ZnDraw, **kwargs):
        pass


class DemoSelection(Extension):
    category = ExtensionType.SELECTION

    def run(self, vis: ZnDraw, **kwargs):
        pass


def test_room_extensions_overview_empty(server, redis_client):
    """Test room extensions overview with no extensions registered."""
    room = "testroom"
    auth_headers = get_jwt_auth_headers(server)

    # Join room
    requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=auth_headers)

    # Get extensions overview
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert "extensions" in data
    assert "summary" in data
    assert data["extensions"] == []
    assert data["summary"]["total_extensions"] == 0
    assert data["summary"]["active_workers"] == 0
    assert data["summary"]["total_jobs_24h"] == 0


def test_room_extensions_overview_with_extension(server, redis_client):
    """Test room extensions overview with a registered extension and job."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Create worker and register extension
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job so the extension appears in overview
    requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )

    # Get extensions overview
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    assert len(data["extensions"]) == 1

    ext = data["extensions"][0]
    assert ext["name"] == "DemoModifier"
    assert ext["category"] == "modifiers"
    assert ext["provider"] == "client"  # Client-side extension
    assert "workers" in ext
    assert ext["workers"]["idle_count"] == 1
    assert "analytics" in ext
    assert ext["analytics"]["total_jobs"] == 1

    vis.disconnect()


def test_room_extensions_overview_with_jobs(server, redis_client):
    """Test room extensions overview with jobs."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Create worker
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job and complete it
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=auth_headers,
    )
    assert response.status_code == 200
    job_id = response.json()["jobId"]

    # Pick up and complete the job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid}
    )
    assert response.status_code == 200

    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "completed", "workerId": vis.sid},
    )
    assert response.status_code == 200

    # Get extensions overview
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )

    assert response.status_code == 200
    data = response.json()
    ext = data["extensions"][0]
    assert ext["analytics"]["total_jobs"] == 1
    assert ext["analytics"]["success_rate"] == 100.0
    assert len(ext["recent_jobs"]) == 1
    assert ext["recent_jobs"][0]["status"] == "completed"

    vis.disconnect()


def test_room_extensions_overview_category_filter(server, redis_client):
    """Test room extensions overview with category filter."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Register extensions from different categories
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)
    vis.register_extension(DemoSelection)

    # Submit jobs for both extensions
    requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    requests.post(
        f"{server}/api/rooms/{room}/extensions/selections/DemoSelection/submit",
        json={"data": {}, "userId": user},
        headers=auth_headers,
    )

    # Get all extensions
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )
    assert response.status_code == 200
    assert len(response.json()["extensions"]) == 2

    # Filter by modifiers
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview?category=modifiers",
        headers=auth_headers,
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["extensions"]) == 1
    assert data["extensions"][0]["category"] == "modifiers"

    # Filter by selections
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview?category=selections",
        headers=auth_headers,
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["extensions"]) == 1
    assert data["extensions"][0]["category"] == "selections"

    vis.disconnect()


def test_room_extensions_overview_search_filter(server, redis_client):
    """Test room extensions overview with search filter."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job so extension appears
    requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )

    # Search for "Modifier"
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview?search=Modifier",
        headers=auth_headers,
    )
    assert response.status_code == 200
    assert len(response.json()["extensions"]) == 1

    # Search for non-existent extension
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview?search=NonExistent",
        headers=auth_headers,
    )
    assert response.status_code == 200
    assert len(response.json()["extensions"]) == 0

    vis.disconnect()


def test_global_extensions_overview_empty(server, redis_client):
    """Test global extensions overview with no extensions."""
    auth_headers = get_jwt_auth_headers(server)

    response = requests.get(f"{server}/api/extensions", headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert "extensions" in data
    assert data["extensions"] == []


def test_global_extensions_overview_with_extensions(server, redis_client):
    """Test global extensions overview with extensions in multiple rooms."""
    room1 = "room1"
    room2 = "room2"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Register extension in room1
    vis1 = ZnDraw(url=server, room=room1, user=user, auto_pickup_jobs=False)
    vis1.register_extension(DemoModifier)

    # Register same extension in room2
    vis2 = ZnDraw(url=server, room=room2, user=user, auto_pickup_jobs=False)
    vis2.register_extension(DemoModifier)

    # Submit jobs in both rooms so extensions appear
    requests.post(
        f"{server}/api/rooms/{room1}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    requests.post(
        f"{server}/api/rooms/{room2}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 2}, "userId": user},
        headers=auth_headers,
    )

    # Get global overview
    response = requests.get(f"{server}/api/extensions", headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert len(data["extensions"]) == 1

    ext = data["extensions"][0]
    assert ext["name"] == "DemoModifier"
    assert ext["category"] == "modifiers"
    assert len(ext["rooms"]) == 2
    assert set(ext["rooms"]) == {room1, room2}

    vis1.disconnect()
    vis2.disconnect()


def test_global_extensions_overview_ignores_non_schema_keys(server, redis_client):
    """Test that global overview scans job keys, not schema keys.

    The global overview now scans for job keys:
    - room:testroom:extension:modifiers:DemoModifier:jobs (SET of job IDs)

    And should ignore schema/worker keys:
    - room:testroom:extensions:modifiers (HASH - schema key)
    - room:testroom:extensions:modifiers:DemoModifier:idle_workers (SET - worker key)
    """
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Register extension which creates multiple Redis keys
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job which creates the job key
    requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )

    # Verify various keys exist
    assert redis_client.exists(f"room:{room}:extensions:modifiers")  # HASH - schema key
    assert redis_client.exists(
        f"room:{room}:extensions:modifiers:DemoModifier:idle_workers"
    )  # SET - worker key
    assert redis_client.exists(
        f"room:{room}:extension:modifiers:DemoModifier:jobs"
    )  # SET - job key

    # Global overview should find the extension via job keys
    response = requests.get(f"{server}/api/extensions", headers=auth_headers)

    assert response.status_code == 200
    data = response.json()
    assert len(data["extensions"]) == 1
    assert data["extensions"][0]["name"] == "DemoModifier"

    vis.disconnect()


def test_extension_detailed_analytics_empty(server, redis_client):
    """Test detailed analytics with no jobs."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/analytics",
        headers=auth_headers,
    )

    assert response.status_code == 200
    data = response.json()
    assert "daily_stats" in data
    assert "total_stats" in data
    assert "error_breakdown" in data
    assert data["total_stats"]["total_jobs"] == 0
    assert len(data["daily_stats"]) == 7  # Default 7 days

    vis.disconnect()


def test_extension_detailed_analytics_with_jobs(server, redis_client):
    """Test detailed analytics with completed jobs."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit and complete multiple jobs
    for i in range(5):
        response = requests.post(
            f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
            json={"data": {"parameter": i}, "userId": user},
            headers=auth_headers,
        )
        assert response.status_code == 200
        job_id = response.json()["jobId"]

        # Pick up and complete
        response = requests.post(
            f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid}
        )
        response = requests.put(
            f"{server}/api/rooms/{room}/jobs/{job_id}/status",
            json={"status": "completed", "workerId": vis.sid},
        )

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/analytics",
        headers=auth_headers,
    )

    assert response.status_code == 200
    data = response.json()
    assert data["total_stats"]["total_jobs"] == 5
    assert data["total_stats"]["overall_success_rate"] == 100.0

    vis.disconnect()


def test_extension_detailed_analytics_auto_days(server, redis_client):
    """Test detailed analytics automatically calculates days from job history."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    assert response.status_code == 200

    # Get analytics - should auto-calculate days based on job history
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/analytics",
        headers=auth_headers,
    )
    assert response.status_code == 200
    data = response.json()

    # Should have at least 1 day of stats (today)
    assert len(data["daily_stats"]) >= 1
    # Should have the job in total stats
    assert data["total_stats"]["total_jobs"] == 1

    vis.disconnect()


def test_extension_detailed_analytics_error_breakdown(server, redis_client):
    """Test error breakdown in detailed analytics."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit and fail a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    job_id = response.json()["jobId"]

    # Pick up and fail
    requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid})
    requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "failed", "workerId": vis.sid, "error": "Test error message"},
    )

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/analytics",
        headers=auth_headers,
    )

    assert response.status_code == 200
    data = response.json()
    assert len(data["error_breakdown"]) == 1
    assert data["error_breakdown"][0]["error"] == "Test error message"
    assert data["error_breakdown"][0]["count"] == 1

    vis.disconnect()


def test_duration_calculation_in_overview(server, redis_client):
    """Test that wait_time_ms and execution_time_ms are calculated correctly."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    job_id = response.json()["jobId"]

    # Pick up job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid}
    )
    assert response.status_code == 200

    # Complete job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "completed", "workerId": vis.sid},
    )
    assert response.status_code == 200

    # Check overview for duration data
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )
    assert response.status_code == 200

    ext = response.json()["extensions"][0]
    # Should have average times calculated
    assert ext["analytics"]["avg_wait_time_ms"] > 0
    assert ext["analytics"]["avg_execution_time_ms"] >= 0

    # Check recent jobs have duration data
    assert len(ext["recent_jobs"]) == 1
    job = ext["recent_jobs"][0]
    assert job["wait_time_ms"] is not None
    assert job["execution_time_ms"] is not None

    vis.disconnect()


def test_room_overview_filters_extensions_with_no_jobs(server, redis_client):
    """Test that extensions with 0 jobs are not shown in room overview."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Register two extensions
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)
    vis.register_extension(DemoSelection)

    # Submit and complete a job only for DemoModifier
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    job_id = response.json()["jobId"]

    requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid})
    requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "completed", "workerId": vis.sid},
    )

    # Get overview - should only show DemoModifier (has 1 job)
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )
    assert response.status_code == 200
    data = response.json()

    # Only DemoModifier should be shown (DemoSelection has 0 jobs)
    assert len(data["extensions"]) == 1
    assert data["extensions"][0]["name"] == "DemoModifier"
    assert data["extensions"][0]["analytics"]["total_jobs"] == 1

    vis.disconnect()


def test_running_jobs_contribute_to_wait_time(server, redis_client):
    """Test that running jobs contribute to avg_wait_time_ms in overview."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit and start a job (but don't complete it - leave it running)
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    job_id = response.json()["jobId"]

    # Pick up job (starts it, calculates wait_time_ms)
    requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid})

    # Get overview - running job should have wait_time_ms > 0
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )
    assert response.status_code == 200
    ext = response.json()["extensions"][0]

    # Running job should contribute to wait time average
    assert ext["analytics"]["avg_wait_time_ms"] > 0
    # But not to execution time (job not completed yet)
    assert ext["analytics"]["avg_execution_time_ms"] == 0.0
    # Success rate should be 0 (no completed jobs)
    assert ext["analytics"]["success_rate"] == 0.0

    vis.disconnect()


def test_global_overview_finds_extensions_via_job_keys(server, redis_client):
    """Test that global overview finds extensions by scanning job keys."""
    room1 = "room1"
    room2 = "room2"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Create jobs in two different rooms
    vis1 = ZnDraw(url=server, room=room1, user=user, auto_pickup_jobs=False)
    vis1.register_extension(DemoModifier)

    vis2 = ZnDraw(url=server, room=room2, user=user, auto_pickup_jobs=False)
    vis2.register_extension(DemoSelection)

    # Submit jobs
    requests.post(
        f"{server}/api/rooms/{room1}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )

    requests.post(
        f"{server}/api/rooms/{room2}/extensions/selections/DemoSelection/submit",
        json={"data": {}, "userId": user},
        headers=auth_headers,
    )

    # Global overview should find both extensions via job keys
    response = requests.get(f"{server}/api/extensions", headers=auth_headers)
    assert response.status_code == 200
    data = response.json()

    # Should find both extensions
    assert len(data["extensions"]) == 2
    extension_names = {ext["name"] for ext in data["extensions"]}
    assert "DemoModifier" in extension_names
    assert "DemoSelection" in extension_names

    # Each extension should show correct room
    modifier = next(e for e in data["extensions"] if e["name"] == "DemoModifier")
    selection = next(e for e in data["extensions"] if e["name"] == "DemoSelection")

    assert room1 in modifier["rooms"]
    assert room2 in selection["rooms"]

    vis1.disconnect()
    vis2.disconnect()


def test_celery_worker_endpoint_format(server, redis_client):
    """Test that job completion uses correct endpoint format."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Submit a job for a server-side extension
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/Empty/submit",
        json={"data": {}, "userId": user},
        headers=auth_headers,
    )
    assert response.status_code == 200
    job_id = response.json()["jobId"]

    # Start the job
    worker_id = "celery:test-worker"
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next",
        json={"workerId": worker_id},
        headers=auth_headers,
    )
    assert response.status_code == 200

    # Complete using the correct endpoint format
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "completed", "workerId": worker_id},
        headers=auth_headers,
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # Verify job is completed
    response = requests.get(
        f"{server}/api/rooms/{room}/jobs/{job_id}",
        headers=auth_headers,
    )
    assert response.status_code == 200
    assert response.json()["status"] == "completed"


def test_celery_worker_failure_endpoint_format(server, redis_client):
    """Test that job failure uses correct endpoint format."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/Empty/submit",
        json={"data": {}, "userId": user},
        headers=auth_headers,
    )
    job_id = response.json()["jobId"]

    # Start the job
    worker_id = "celery:test-worker"
    requests.post(
        f"{server}/api/rooms/{room}/jobs/next",
        json={"workerId": worker_id},
        headers=auth_headers,
    )

    # Fail using the correct endpoint format
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "failed", "error": "Test error", "workerId": worker_id},
        headers=auth_headers,
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # Verify job is failed
    response = requests.get(
        f"{server}/api/rooms/{room}/jobs/{job_id}",
        headers=auth_headers,
    )
    assert response.status_code == 200
    job_data = response.json()
    assert job_data["status"] == "failed"
    assert job_data["error"] == "Test error"


def test_extension_history_persists_after_client_disconnect(server, redis_client):
    """Test that extension job history remains visible after client disconnects.

    When a client extension worker disconnects:
    - The extension schema should be deleted (correct - no workers to handle new jobs)
    - But job history should still be accessible in the room extensions overview
    - The extension should appear as read-only with historical data
    """
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Register extension and submit/complete a job
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=auth_headers,
    )
    assert response.status_code == 200
    job_id = response.json()["jobId"]

    # Pick up and complete the job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid}
    )
    assert response.status_code == 200

    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "completed", "workerId": vis.sid},
    )
    assert response.status_code == 200

    # Verify extension is visible with job history BEFORE disconnect
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["extensions"]) == 1
    assert data["extensions"][0]["name"] == "DemoModifier"
    assert data["extensions"][0]["analytics"]["total_jobs"] == 1
    assert data["extensions"][0]["workers"]["idle_count"] == 1

    # Verify the schema exists in Redis before disconnect
    schema_key = f"room:{room}:extensions:modifiers"
    assert redis_client.hexists(schema_key, "DemoModifier")

    # Disconnect the client (this should delete the schema)
    vis.disconnect()

    # Give a moment for the disconnect handler to complete
    import time

    time.sleep(0.1)

    # Verify the schema is deleted after disconnect
    assert not redis_client.hexists(schema_key, "DemoModifier")

    # Verify job history key still exists
    job_key = f"room:{room}:extension:modifiers:DemoModifier:jobs"
    assert redis_client.exists(job_key)
    assert redis_client.scard(job_key) == 1

    # THE FIX: Extension should STILL be visible in room overview with historical data
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/overview", headers=auth_headers
    )
    assert response.status_code == 200
    data = response.json()

    # Extension should still appear with job history
    assert len(data["extensions"]) == 1
    ext = data["extensions"][0]
    assert ext["name"] == "DemoModifier"
    assert ext["category"] == "modifiers"
    assert ext["analytics"]["total_jobs"] == 1
    assert ext["analytics"]["success_rate"] == 100.0

    # Workers should show as 0 (client disconnected)
    assert ext["workers"]["idle_count"] == 0
    assert ext["workers"]["progressing_count"] == 0

    # Recent jobs should still be accessible
    assert len(ext["recent_jobs"]) == 1
    assert ext["recent_jobs"][0]["status"] == "completed"


def test_client_extension_failure_logs_to_ui(server, redis_client):
    """Test that failed client extension jobs attempt to log errors to UI via vis.log()."""
    room = "testroom"
    user = "testuser"
    auth_headers = get_jwt_auth_headers(server, user)

    # Register extension
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(DemoModifier)

    # Submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/DemoModifier/submit",
        json={"data": {"parameter": 1}, "userId": user},
        headers=auth_headers,
    )
    job_id = response.json()["jobId"]

    # Pick up job
    requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.sid})

    # Fail the job with an error
    error_message = "Test error: Something went wrong"
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{job_id}/status",
        json={"status": "failed", "workerId": vis.sid, "error": error_message},
        headers=auth_headers,
    )
    assert response.status_code == 200

    # Verify the job is marked as failed
    job_data = redis_client.hgetall(f"job:{job_id}")
    assert job_data["status"] == "failed"
    assert job_data["error"] == error_message

    # The vis.log() call in the backend will attempt to create a ZnDraw client
    # In production, this displays the error to users
    # In tests, we just verify the job failure was processed correctly

    vis.disconnect()
