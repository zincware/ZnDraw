# tests/test_router_task_submit.py
"""Tests for task submission endpoint using shared fixtures."""

from unittest.mock import AsyncMock, MagicMock

from zndraw_joblib.exceptions import ProblemDetail
from zndraw_joblib.registry import InternalRegistry
from zndraw_joblib.schemas import TaskResponse


def test_submit_task(seeded_client):
    response = seeded_client.post(
        "/v1/joblib/rooms/room_123/tasks/@global:modifiers:Rotate",
        json={"payload": {"angle": 90}},
    )
    assert response.status_code == 202
    assert "Location" in response.headers
    data = TaskResponse.model_validate(response.json())
    assert data.status.value == "pending"
    assert data.payload == {"angle": 90}
    assert data.room_id == "room_123"


def test_submit_task_job_not_found(seeded_client):
    response = seeded_client.post(
        "/v1/joblib/rooms/room_123/tasks/@global:modifiers:NonExistent",
        json={"payload": {}},
    )
    assert response.status_code == 404
    error = ProblemDetail.model_validate(response.json())
    assert error.status == 404


def test_submit_internal_task_dispatches_to_taskiq(app, client):
    """Submitting a task for an @internal job dispatches via taskiq."""
    # Register an @internal job
    resp = client.put(
        "/v1/joblib/rooms/@internal/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code in (200, 201)

    # Set up mock internal registry on app.state
    mock_task_handle = MagicMock()
    mock_task_handle.kiq = AsyncMock()
    registry = InternalRegistry(
        tasks={"@internal:modifiers:Rotate": mock_task_handle},
        extensions={},
    )
    app.state.internal_registry = registry

    # Submit a task
    resp = client.post(
        "/v1/joblib/rooms/test-room/tasks/@internal:modifiers:Rotate",
        json={"payload": {"angle": 90}},
    )
    assert resp.status_code == 202

    # Verify taskiq dispatch was called
    mock_task_handle.kiq.assert_called_once()
    call_kwargs = mock_task_handle.kiq.call_args.kwargs
    assert call_kwargs["room_id"] == "test-room"
    assert "task_id" in call_kwargs


def test_submit_internal_task_no_registry_returns_503(app, client):
    """Submitting to @internal job without registry returns 503."""
    # Register an @internal job
    resp = client.put(
        "/v1/joblib/rooms/@internal/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code in (200, 201)

    # Ensure no registry is set
    if hasattr(app.state, "internal_registry"):
        delattr(app.state, "internal_registry")

    resp = client.post(
        "/v1/joblib/rooms/test-room/tasks/@internal:modifiers:Rotate",
        json={"payload": {}},
    )
    assert resp.status_code == 503


def test_submit_external_task_unchanged(seeded_client):
    """External task submission still works as before (no taskiq dispatch)."""
    resp = seeded_client.post(
        "/v1/joblib/rooms/test-room/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    assert resp.status_code == 202
    assert resp.json()["status"] == "pending"


def test_submit_task_malformed_job_name_no_colons(seeded_client):
    """Submitting with a job name that has no colons returns 404."""
    response = seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/noColons",
        json={"payload": {}},
    )
    assert response.status_code == 404
    error = ProblemDetail.model_validate(response.json())
    assert "Invalid job name format" in error.detail


def test_submit_task_malformed_job_name_one_colon(seeded_client):
    """Submitting with a job name that has only one colon returns 404."""
    response = seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/one:part",
        json={"payload": {}},
    )
    assert response.status_code == 404
    error = ProblemDetail.model_validate(response.json())
    assert "Invalid job name format" in error.detail


def test_submit_task_cross_room_allowed(client):
    """Job registered in room_A can be submitted to target room_B."""
    # Register a room-scoped job in room_A
    client.put(
        "/v1/joblib/rooms/room_A/jobs",
        json={"category": "modifiers", "name": "Private", "schema": {}},
    )

    # Submit a task targeting room_B using room_A's job
    response = client.post(
        "/v1/joblib/rooms/room_B/tasks/room_A:modifiers:Private",
        json={"payload": {}},
    )
    assert response.status_code == 202
    data = TaskResponse.model_validate(response.json())
    assert data.room_id == "room_B"


def test_submit_task_invalid_room_id(seeded_client):
    """Submitting to a room with @ in the name returns 400."""
    response = seeded_client.post(
        "/v1/joblib/rooms/bad@room/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    assert response.status_code == 400
    error = ProblemDetail.model_validate(response.json())
    assert "invalid characters" in error.detail


def test_submit_task_empty_payload(seeded_client):
    """Submitting with no payload key defaults to empty dict."""
    response = seeded_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Rotate",
        json={},
    )
    assert response.status_code == 202
    data = TaskResponse.model_validate(response.json())
    assert data.payload == {}
