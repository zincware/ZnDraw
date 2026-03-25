# tests/test_models.py
from datetime import datetime, timedelta, timezone
from uuid import UUID

from zndraw_joblib.models import (
    Job,
    Task,
    TaskStatus,
    Worker,
)


def test_task_status_enum():
    assert TaskStatus.PENDING.value == "pending"
    assert TaskStatus.CLAIMED.value == "claimed"
    assert TaskStatus.RUNNING.value == "running"
    assert TaskStatus.COMPLETED.value == "completed"
    assert TaskStatus.FAILED.value == "failed"
    assert TaskStatus.CANCELLED.value == "cancelled"


def test_job_full_name():
    job = Job(room_id="@global", category="modifiers", name="Rotate", schema_={})
    assert job.full_name == "@global:modifiers:Rotate"


def test_job_full_name_private():
    job = Job(room_id="room_123", category="selections", name="All", schema_={})
    assert job.full_name == "room_123:selections:All"


def test_worker_is_alive():
    user_id = UUID("12345678-1234-5678-1234-567812345678")
    worker = Worker(user_id=user_id, last_heartbeat=datetime.now(timezone.utc))
    assert worker.is_alive(timedelta(seconds=60)) is True


def test_worker_is_dead():
    user_id = UUID("12345678-1234-5678-1234-567812345678")
    old_time = datetime.now(timezone.utc) - timedelta(seconds=120)
    worker = Worker(user_id=user_id, last_heartbeat=old_time)
    assert worker.is_alive(timedelta(seconds=60)) is False


def test_task_has_uuid_id():
    # SQLAlchemy mapped columns with default=uuid4 only generate ID when added to session
    # We test that the column accepts UUID values
    task_id = UUID("12345678-1234-5678-1234-567812345678")
    job_id = UUID("12345678-1234-5678-1234-567812345679")
    task = Task(id=task_id, job_id=job_id, room_id="room_1")
    assert isinstance(task.id, UUID)
    assert task.id == task_id


def test_task_status_accepts_pending():
    # Task can be created with PENDING status
    # Note: SQLAlchemy mapped column defaults only apply when added to session
    # Default behavior is tested via integration tests in test_router_*.py
    job_id = UUID("12345678-1234-5678-1234-567812345678")
    task = Task(job_id=job_id, room_id="room_1", status=TaskStatus.PENDING)
    assert task.status == TaskStatus.PENDING
