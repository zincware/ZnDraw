# tests/test_schemas.py
from datetime import datetime, timezone
from uuid import UUID

from zndraw_joblib.models import TaskStatus
from zndraw_joblib.schemas import (
    JobRegisterRequest,
    JobResponse,
    TaskClaimResponse,
    TaskResponse,
    TaskSubmitRequest,
    TaskUpdateRequest,
)


def test_job_register_request():
    req = JobRegisterRequest(category="modifiers", name="Rotate", schema={"angle": 0})
    assert req.category == "modifiers"
    assert req.name == "Rotate"
    assert req.schema_ == {"angle": 0}


def test_job_response():
    worker1_id = UUID("11111111-1111-1111-1111-111111111111")
    worker2_id = UUID("22222222-2222-2222-2222-222222222222")
    resp = JobResponse(
        id=UUID("12345678-1234-5678-1234-567812345678"),
        room_id="@global",
        category="modifiers",
        name="Rotate",
        full_name="@global:modifiers:Rotate",
        schema={"angle": 0},
        workers=[worker1_id, worker2_id],
    )
    assert resp.full_name == "@global:modifiers:Rotate"
    assert len(resp.workers) == 2
    assert set(resp.workers) == {worker1_id, worker2_id}


def test_task_submit_request():
    req = TaskSubmitRequest(payload={"angle": 90})
    assert req.payload == {"angle": 90}


def test_task_response():
    resp = TaskResponse(
        id=UUID("12345678-1234-5678-1234-567812345678"),
        job_name="@global:modifiers:Rotate",
        room_id="room_1",
        status=TaskStatus.PENDING,
        created_at=datetime.now(timezone.utc),
    )
    assert resp.status == TaskStatus.PENDING


def test_task_update_request_valid_status():
    req = TaskUpdateRequest(status=TaskStatus.RUNNING)
    assert req.status == TaskStatus.RUNNING


def test_task_claim_response_with_task():
    resp = TaskClaimResponse(
        task=TaskResponse(
            id=UUID("12345678-1234-5678-1234-567812345678"),
            job_name="@global:modifiers:Rotate",
            room_id="room_1",
            status=TaskStatus.CLAIMED,
            created_at=datetime.now(timezone.utc),
        )
    )
    assert resp.task is not None


def test_task_claim_response_empty():
    resp = TaskClaimResponse(task=None)
    assert resp.task is None
