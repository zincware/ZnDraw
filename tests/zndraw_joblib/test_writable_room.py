# tests/test_writable_room.py
"""Tests for WritableRoomDep — verifies the DI override point and default behavior."""

from fastapi import HTTPException, Path
from starlette.testclient import TestClient

from zndraw_joblib.dependencies import verify_writable_room
from zndraw_joblib.exceptions import ProblemDetail


def _locked_room_override(room_id: str = Path()) -> str:
    """Dummy DI that rejects room_id 'a' with 423 Locked."""
    if room_id == "a":
        raise HTTPException(status_code=423, detail="Room is locked")
    return room_id


def test_register_job_blocked_by_writable_room(app):
    """register_job respects WritableRoomDep override — locked room returns 423."""
    app.dependency_overrides[verify_writable_room] = _locked_room_override
    locked_client = TestClient(app)

    resp = locked_client.put(
        "/v1/joblib/rooms/a/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code == 423


def test_register_job_allowed_by_writable_room(app):
    """register_job proceeds when WritableRoomDep override allows the room."""
    app.dependency_overrides[verify_writable_room] = _locked_room_override
    allowed_client = TestClient(app)

    resp = allowed_client.put(
        "/v1/joblib/rooms/b/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code == 201


def test_submit_task_blocked_by_writable_room(app, client):
    """submit_task respects WritableRoomDep override — locked room returns 423."""
    # First register a job via an unlocked room so we have something to submit to
    resp = client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code == 201

    # Now lock and try to submit a task to room "a"
    app.dependency_overrides[verify_writable_room] = _locked_room_override
    locked_client = TestClient(app)

    resp = locked_client.post(
        "/v1/joblib/rooms/a/tasks/@global:modifiers:Rotate",
        json={"payload": {"angle": 90}},
    )
    assert resp.status_code == 423


def test_submit_task_allowed_by_writable_room(app, client):
    """submit_task proceeds when WritableRoomDep override allows the room."""
    resp = client.put(
        "/v1/joblib/rooms/@global/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code == 201

    app.dependency_overrides[verify_writable_room] = _locked_room_override
    allowed_client = TestClient(app)

    resp = allowed_client.post(
        "/v1/joblib/rooms/b/tasks/@global:modifiers:Rotate",
        json={"payload": {"angle": 90}},
    )
    assert resp.status_code == 202


def test_read_endpoints_not_affected_by_writable_room(app):
    """Read endpoints (list_jobs, list_tasks) don't use WritableRoomDep."""
    app.dependency_overrides[verify_writable_room] = _locked_room_override
    locked_client = TestClient(app)

    # list_jobs for room "a" should still work (200, not 423)
    resp = locked_client.get("/v1/joblib/rooms/a/jobs")
    assert resp.status_code == 200

    # list_tasks for room "a" should still work
    resp = locked_client.get("/v1/joblib/rooms/a/tasks")
    assert resp.status_code == 200


# --- Default behavior (no override) ---


def test_default_rejects_invalid_room_id_on_register(client):
    """Default verify_writable_room rejects room IDs with @ or : characters."""
    resp = client.put(
        "/v1/joblib/rooms/bad@room/jobs",
        json={"category": "modifiers", "name": "Rotate", "schema": {}},
    )
    assert resp.status_code == 400
    error = ProblemDetail.model_validate(resp.json())
    assert "invalid characters" in error.detail.lower()


def test_default_rejects_invalid_room_id_on_submit(seeded_client):
    """Default verify_writable_room rejects room IDs with @ or : on submit_task."""
    resp = seeded_client.post(
        "/v1/joblib/rooms/bad:room/tasks/@global:modifiers:Rotate",
        json={"payload": {}},
    )
    assert resp.status_code == 400
