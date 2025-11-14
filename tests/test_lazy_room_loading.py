"""Integration tests for lazy room data loading.

Tests that all REST endpoints for room data work correctly after the join
response was changed to be minimal. This ensures the lazy loading refactoring
maintains data consistency and handles edge cases properly.
"""

import pytest
import requests
from conftest import get_jwt_auth_headers

from zndraw.zndraw import ZnDraw


def test_lazy_loading_empty_room(server):
    """Test lazy loading endpoints with an empty room."""
    room = "test-lazy-empty"
    headers = get_jwt_auth_headers(server)

    # Step 1: Join room (minimal response)
    response = requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room
    assert data["created"] is True
    assert "sessionId" in data
    assert "userName" in data

    # Step 2: Fetch room info
    response = requests.get(f"{server}/api/rooms/{room}", headers=headers)
    assert response.status_code == 200
    room_info = response.json()
    assert room_info["id"] == room
    assert room_info["frameCount"] == 0
    assert room_info["locked"] is False

    # Step 3: Fetch selections
    response = requests.get(f"{server}/api/rooms/{room}/selections", headers=headers)
    assert response.status_code == 200
    selections_data = response.json()
    assert selections_data["selections"] == {}
    assert selections_data["groups"] == {}
    assert selections_data["activeGroup"] is None

    # Step 4: Fetch frame selection
    response = requests.get(f"{server}/api/rooms/{room}/frame-selection", headers=headers)
    assert response.status_code == 200
    frame_selection_data = response.json()
    assert frame_selection_data["frameSelection"] is None

    # Step 5: Fetch current step
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=headers)
    assert response.status_code == 200
    step_data = response.json()
    assert step_data["step"] is None or step_data["step"] == 0
    assert step_data["totalFrames"] == 0

    # Step 6: Fetch bookmarks
    response = requests.get(f"{server}/api/rooms/{room}/bookmarks", headers=headers)
    assert response.status_code == 200
    bookmarks_data = response.json()
    assert bookmarks_data["bookmarks"] == {} or bookmarks_data["bookmarks"] is None

    # Step 7: Fetch geometries
    response = requests.get(f"{server}/api/rooms/{room}/geometries", headers=headers)
    assert response.status_code == 200
    geometries_data = response.json()
    # Geometries might contain default schemas even for empty room
    assert isinstance(geometries_data["geometries"], dict)

    # Step 8: Fetch user settings
    response = requests.get(f"{server}/api/rooms/{room}/settings", headers=headers)
    assert response.status_code == 200
    settings_data = response.json()
    assert "settings" in settings_data
    # Settings might be empty dict or have default categories
    assert isinstance(settings_data["settings"], dict)


def test_lazy_loading_with_data(server, s22):
    """Test lazy loading endpoints with a room that has data."""
    room = "test-lazy-with-data"
    headers = get_jwt_auth_headers(server)

    # Setup: Create room with frames
    vis = ZnDraw(url=server, room=room, user="user1")
    vis.extend(s22[:5])  # Add 5 frames

    # Add a bookmark
    vis.bookmarks[2] = "Frame 2"

    # Step 1: Join room as different user
    response = requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["created"] is False  # Room already exists

    # Step 2: Fetch room info
    response = requests.get(f"{server}/api/rooms/{room}", headers=headers)
    assert response.status_code == 200
    room_info = response.json()
    assert room_info["id"] == room
    assert room_info["frameCount"] == 5

    # Step 3: Fetch current step
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=headers)
    assert response.status_code == 200
    step_data = response.json()
    assert step_data["totalFrames"] == 5
    # Current step might be None or a valid index
    if step_data["step"] is not None:
        assert 0 <= step_data["step"] < 5

    # Step 4: Fetch bookmarks
    response = requests.get(f"{server}/api/rooms/{room}/bookmarks", headers=headers)
    assert response.status_code == 200
    bookmarks_data = response.json()
    assert bookmarks_data["bookmarks"] is not None
    assert "2" in bookmarks_data["bookmarks"] or 2 in bookmarks_data["bookmarks"]
    # Convert key to string if needed for comparison
    bookmark_key = "2" if "2" in bookmarks_data["bookmarks"] else 2
    assert bookmarks_data["bookmarks"][bookmark_key] == "Frame 2"

    # Step 5: Fetch geometries (should be empty initially)
    response = requests.get(f"{server}/api/rooms/{room}/geometries", headers=headers)
    assert response.status_code == 200
    geometries_data = response.json()
    assert isinstance(geometries_data["geometries"], dict)


def test_lazy_loading_auth_required(server):
    """Test that all lazy loading endpoints require authentication."""
    room = "test-lazy-auth"
    headers = get_jwt_auth_headers(server)

    # Create room first
    response = requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=headers)
    assert response.status_code == 200

    # Try to access endpoints without auth - should fail with 401
    endpoints = [
        f"/api/rooms/{room}",
        f"/api/rooms/{room}/selections",
        f"/api/rooms/{room}/frame-selection",
        f"/api/rooms/{room}/step",
        f"/api/rooms/{room}/geometries",
        f"/api/rooms/{room}/settings",
    ]

    for endpoint in endpoints:
        response = requests.get(f"{server}{endpoint}")
        assert response.status_code == 401, f"Endpoint {endpoint} should require auth"


def test_lazy_loading_nonexistent_room(server):
    """Test lazy loading endpoints with a room that doesn't exist."""
    room = "nonexistent-room"
    headers = get_jwt_auth_headers(server)

    # Room info should return 404
    response = requests.get(f"{server}/api/rooms/{room}", headers=headers)
    assert response.status_code == 404

    # Other endpoints might return empty data or 404 depending on implementation
    # We just verify they don't crash
    endpoints = [
        f"/api/rooms/{room}/selections",
        f"/api/rooms/{room}/frame-selection",
        f"/api/rooms/{room}/step",
        f"/api/rooms/{room}/geometries",
    ]

    for endpoint in endpoints:
        response = requests.get(f"{server}{endpoint}", headers=headers)
        # Should return either 200 with empty data or 404, but not 500
        assert response.status_code in [200, 404], f"Endpoint {endpoint} returned {response.status_code}"


def test_step_clamping_with_deleted_frames(server, s22):
    """Test that step clamping updates Redis when frames are deleted."""
    room = "test-step-clamping"
    headers = get_jwt_auth_headers(server)

    # Setup: Create room with 10 frames
    vis = ZnDraw(url=server, room=room, user="user1")
    vis.extend(s22[:10])

    # Set current step to frame 8
    vis.step = 8

    # Verify step is 8
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["step"] == 8
    assert data["totalFrames"] == 10

    # Delete frames so only 5 remain (indices 0-4)
    # This simulates the scenario where current step becomes invalid
    del vis[5:]

    # Fetch step again - should be clamped to 4 (last valid index)
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["step"] == 4  # Clamped from 8 to 4
    assert data["totalFrames"] == 5

    # Verify that subsequent fetches return the same clamped value
    # (meaning Redis was updated, not just the response)
    response = requests.get(f"{server}/api/rooms/{room}/step", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data["step"] == 4
    assert data["totalFrames"] == 5


def test_parallel_lazy_loading(server, s22):
    """Test that all lazy loading endpoints can be called in parallel."""
    import concurrent.futures

    room = "test-parallel-loading"
    headers = get_jwt_auth_headers(server)

    # Setup: Create room with data
    vis = ZnDraw(url=server, room=room, user="user1")
    vis.extend(s22[:3])
    vis.bookmarks[1] = "Parallel test"

    # Join room
    response = requests.post(f"{server}/api/rooms/{room}/join", json={}, headers=headers)
    assert response.status_code == 200

    # Define all endpoints to fetch
    endpoints = [
        f"{server}/api/rooms/{room}",
        f"{server}/api/rooms/{room}/selections",
        f"{server}/api/rooms/{room}/frame-selection",
        f"{server}/api/rooms/{room}/step",
        f"{server}/api/rooms/{room}/bookmarks",
        f"{server}/api/rooms/{room}/geometries",
        f"{server}/api/rooms/{room}/settings",
    ]

    # Fetch all endpoints in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(endpoints)) as executor:
        futures = {executor.submit(requests.get, url, headers=headers): url for url in endpoints}
        results = {}
        for future in concurrent.futures.as_completed(futures):
            url = futures[future]
            try:
                response = future.result(timeout=10)
                results[url] = response
            except Exception as e:
                pytest.fail(f"Failed to fetch {url}: {e}")

    # Verify all requests succeeded
    for url, response in results.items():
        assert response.status_code == 200, f"Failed to fetch {url}: {response.status_code}"
        data = response.json()
        assert data is not None, f"Empty response from {url}"

    # Verify data consistency
    room_info = results[f"{server}/api/rooms/{room}"].json()
    step_data = results[f"{server}/api/rooms/{room}/step"].json()
    bookmarks_data = results[f"{server}/api/rooms/{room}/bookmarks"].json()

    assert room_info["frameCount"] == 3
    assert step_data["totalFrames"] == 3
    assert bookmarks_data["bookmarks"] is not None
