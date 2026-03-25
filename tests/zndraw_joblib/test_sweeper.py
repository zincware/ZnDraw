# tests/test_sweeper.py
"""Tests for the sweeper background task."""

import asyncio
import uuid
from datetime import datetime, timedelta, timezone

import pytest
from sqlalchemy import select

from zndraw_joblib.events import JobsInvalidate, TaskStatusEvent
from zndraw_joblib.models import Job, Task, TaskStatus, Worker, WorkerJobLink
from zndraw_joblib.settings import JobLibSettings
from zndraw_joblib.sweeper import (
    cleanup_stale_workers,
    cleanup_stuck_internal_tasks,
    cleanup_worker,
    run_sweeper,
)


@pytest.mark.asyncio
async def test_cleanup_stale_workers_finds_stale(async_session_factory, test_user_id):
    """Sweeper should clean up workers with old heartbeats."""
    stale_worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        # Create a stale worker (heartbeat 2 minutes ago)
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        worker = Worker(
            id=stale_worker_id, user_id=test_user_id, last_heartbeat=stale_time
        )
        session.add(worker)
        await session.commit()

    # Run cleanup with 60 second timeout
    async with async_session_factory() as session:
        count, _, _ = await cleanup_stale_workers(session, timedelta(seconds=60))
        assert count == 1

    # Verify worker was deleted
    async with async_session_factory() as session:
        result = await session.execute(
            select(Worker).where(Worker.id == stale_worker_id)
        )
        worker = result.scalar_one_or_none()
        assert worker is None


@pytest.mark.asyncio
async def test_cleanup_stale_workers_ignores_alive(async_session_factory, test_user_id):
    """Sweeper should not touch workers with recent heartbeats."""
    alive_worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        # Create an alive worker (heartbeat just now)
        worker = Worker(
            id=alive_worker_id,
            user_id=test_user_id,
            last_heartbeat=datetime.now(timezone.utc),
        )
        session.add(worker)
        await session.commit()

    # Run cleanup with 60 second timeout
    async with async_session_factory() as session:
        count, _, _ = await cleanup_stale_workers(session, timedelta(seconds=60))
        assert count == 0

    # Verify worker still exists
    async with async_session_factory() as session:
        result = await session.execute(
            select(Worker).where(Worker.id == alive_worker_id)
        )
        worker = result.scalar_one_or_none()
        assert worker is not None


@pytest.mark.asyncio
async def test_cleanup_fails_running_tasks(async_session_factory, test_user_id):
    """Sweeper should fail CLAIMED/RUNNING tasks when cleaning up a worker."""
    stale_worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        # Create a stale worker
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        worker = Worker(
            id=stale_worker_id, user_id=test_user_id, last_heartbeat=stale_time
        )
        session.add(worker)

        # Create a job
        job = Job(room_id="@global", category="modifiers", name="TestJob", schema_={})
        session.add(job)
        await session.flush()

        # Link worker to job
        link = WorkerJobLink(worker_id=stale_worker_id, job_id=job.id)
        session.add(link)

        # Create tasks in various states
        task_claimed = Task(
            job_id=job.id,
            room_id="room1",
            status=TaskStatus.CLAIMED,
            worker_id=stale_worker_id,
        )
        task_running = Task(
            job_id=job.id,
            room_id="room1",
            status=TaskStatus.RUNNING,
            worker_id=stale_worker_id,
        )
        task_pending = Task(
            job_id=job.id,
            room_id="room1",
            status=TaskStatus.PENDING,
            worker_id=None,
        )
        session.add_all([task_claimed, task_running, task_pending])
        await session.commit()

        task_claimed_id = task_claimed.id
        task_running_id = task_running.id
        task_pending_id = task_pending.id

    # Run cleanup
    async with async_session_factory() as session:
        count, _, _ = await cleanup_stale_workers(session, timedelta(seconds=60))
        assert count == 1

    # Verify CLAIMED and RUNNING tasks are failed, PENDING is unchanged
    async with async_session_factory() as session:
        result = await session.execute(select(Task).where(Task.id == task_claimed_id))
        claimed = result.scalar_one_or_none()
        result = await session.execute(select(Task).where(Task.id == task_running_id))
        running = result.scalar_one_or_none()
        result = await session.execute(select(Task).where(Task.id == task_pending_id))
        pending = result.scalar_one_or_none()

        assert claimed.status == TaskStatus.FAILED
        assert claimed.error == "Worker disconnected"
        assert claimed.completed_at is not None

        assert running.status == TaskStatus.FAILED
        assert running.error == "Worker disconnected"
        assert running.completed_at is not None

        assert pending.status == TaskStatus.PENDING
        assert pending.error is None


@pytest.mark.asyncio
async def test_cleanup_soft_deletes_orphan_jobs(async_session_factory, test_user_id):
    """Sweeper should soft-delete jobs with no workers and no non-terminal tasks."""
    stale_worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        # Create a stale worker
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        worker = Worker(
            id=stale_worker_id, user_id=test_user_id, last_heartbeat=stale_time
        )
        session.add(worker)

        # Create a job
        job = Job(room_id="@global", category="modifiers", name="OrphanJob", schema_={})
        session.add(job)
        await session.flush()
        job_id = job.id

        # Link worker to job
        link = WorkerJobLink(worker_id=stale_worker_id, job_id=job.id)
        session.add(link)

        # Create a completed task (terminal state)
        task = Task(
            job_id=job.id,
            room_id="room1",
            status=TaskStatus.COMPLETED,
            worker_id=stale_worker_id,
        )
        session.add(task)
        await session.commit()

    # Run cleanup
    async with async_session_factory() as session:
        count, _, _ = await cleanup_stale_workers(session, timedelta(seconds=60))
        assert count == 1

    # Verify job is soft-deleted
    async with async_session_factory() as session:
        result = await session.execute(select(Job).where(Job.id == job_id))
        job = result.scalar_one_or_none()
        assert job is not None
        assert job.deleted is True


@pytest.mark.asyncio
async def test_cleanup_keeps_job_with_pending_tasks(
    async_session_factory, test_user_id
):
    """Sweeper should not soft-delete jobs that have pending tasks."""
    stale_worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        # Create a stale worker
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        worker = Worker(
            id=stale_worker_id, user_id=test_user_id, last_heartbeat=stale_time
        )
        session.add(worker)

        # Create a job
        job = Job(room_id="@global", category="modifiers", name="ActiveJob", schema_={})
        session.add(job)
        await session.flush()
        job_id = job.id

        # Link worker to job
        link = WorkerJobLink(worker_id=stale_worker_id, job_id=job.id)
        session.add(link)

        # Create a pending task (non-terminal state)
        task = Task(
            job_id=job.id,
            room_id="room1",
            status=TaskStatus.PENDING,
            worker_id=None,
        )
        session.add(task)
        await session.commit()

    # Run cleanup
    async with async_session_factory() as session:
        count, _, _ = await cleanup_stale_workers(session, timedelta(seconds=60))
        assert count == 1

    # Verify job is NOT soft-deleted because it has a pending task
    async with async_session_factory() as session:
        result = await session.execute(select(Job).where(Job.id == job_id))
        job = result.scalar_one_or_none()
        assert job is not None
        assert job.deleted is False


@pytest.mark.asyncio
async def test_cleanup_multiple_stale_workers(async_session_factory, test_user_id):
    """Sweeper should clean up multiple stale workers in one pass."""
    stale_ids = [uuid.uuid4() for _ in range(3)]
    alive_id = uuid.uuid4()

    async with async_session_factory() as session:
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        alive_time = datetime.now(timezone.utc)

        # Create 3 stale workers and 1 alive
        for worker_id in stale_ids:
            session.add(
                Worker(id=worker_id, user_id=test_user_id, last_heartbeat=stale_time)
            )
        session.add(
            Worker(id=alive_id, user_id=test_user_id, last_heartbeat=alive_time)
        )
        await session.commit()

    # Run cleanup
    async with async_session_factory() as session:
        count, _, _ = await cleanup_stale_workers(session, timedelta(seconds=60))
        assert count == 3

    # Verify only alive worker remains
    async with async_session_factory() as session:
        result = await session.execute(select(Worker))
        workers = result.scalars().all()
        assert len(workers) == 1
        assert workers[0].id == alive_id


@pytest.mark.asyncio
async def test_cleanup_worker_removes_links(async_session_factory, test_user_id):
    """cleanup_worker should remove worker-job links."""
    worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        # Create worker with a job
        worker = Worker(
            id=worker_id,
            user_id=test_user_id,
            last_heartbeat=datetime.now(timezone.utc),
        )
        session.add(worker)

        job = Job(room_id="@global", category="modifiers", name="TestJob", schema_={})
        session.add(job)
        await session.flush()

        link = WorkerJobLink(worker_id=worker_id, job_id=job.id)
        session.add(link)
        await session.commit()

        # Now cleanup
        result = await session.execute(select(Worker).where(Worker.id == worker_id))
        worker = result.scalar_one_or_none()
        await cleanup_worker(session, worker)
        await session.commit()

    # Verify link is removed
    async with async_session_factory() as session:
        result = await session.execute(
            select(WorkerJobLink).where(WorkerJobLink.worker_id == worker_id)
        )
        links = result.scalars().all()
        assert len(links) == 0


@pytest.mark.asyncio
async def test_cleanup_job_keeps_other_workers(async_session_factory, test_user_id):
    """If a job has other workers, it should not be soft-deleted."""
    stale_worker_id = uuid.uuid4()
    alive_worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        # Create two workers
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        worker1 = Worker(
            id=stale_worker_id, user_id=test_user_id, last_heartbeat=stale_time
        )
        worker2 = Worker(
            id=alive_worker_id,
            user_id=test_user_id,
            last_heartbeat=datetime.now(timezone.utc),
        )
        session.add_all([worker1, worker2])

        # Create a job linked to both
        job = Job(room_id="@global", category="modifiers", name="SharedJob", schema_={})
        session.add(job)
        await session.flush()
        job_id = job.id

        link1 = WorkerJobLink(worker_id=stale_worker_id, job_id=job.id)
        link2 = WorkerJobLink(worker_id=alive_worker_id, job_id=job.id)
        session.add_all([link1, link2])
        await session.commit()

    # Run cleanup
    async with async_session_factory() as session:
        count, _, _ = await cleanup_stale_workers(session, timedelta(seconds=60))
        assert count == 1

    # Verify job is NOT soft-deleted because alive-worker is still linked
    async with async_session_factory() as session:
        result = await session.execute(select(Job).where(Job.id == job_id))
        job = result.scalar_one_or_none()
        assert job is not None
        assert job.deleted is False

        # Verify alive worker's link still exists
        result = await session.execute(
            select(WorkerJobLink).where(WorkerJobLink.worker_id == alive_worker_id)
        )
        link = result.scalar_one_or_none()
        assert link is not None


@pytest.mark.asyncio
async def test_cleanup_stuck_internal_tasks(async_session_factory):
    """Stuck @internal tasks in RUNNING are marked FAILED after timeout."""
    async with async_session_factory() as session:
        job = Job(room_id="@internal", category="modifiers", name="Rotate", schema_={})
        session.add(job)
        await session.flush()

        task = Task(
            job_id=job.id,
            room_id="test-room",
            status=TaskStatus.RUNNING,
            started_at=datetime.now(timezone.utc) - timedelta(hours=2),
        )
        session.add(task)
        await session.commit()
        task_id = task.id

    async with async_session_factory() as session:
        count, _ = await cleanup_stuck_internal_tasks(
            session, timeout=timedelta(hours=1)
        )
        assert count == 1

    async with async_session_factory() as session:
        result = await session.execute(select(Task).where(Task.id == task_id))
        task = result.scalar_one()
        assert task.status == TaskStatus.FAILED
        assert task.error == "Internal worker timeout"
        assert task.completed_at is not None


@pytest.mark.asyncio
async def test_cleanup_stuck_internal_tasks_skips_recent(async_session_factory):
    """Recently started @internal tasks are not cleaned up."""
    async with async_session_factory() as session:
        job = Job(room_id="@internal", category="modifiers", name="Scale", schema_={})
        session.add(job)
        await session.flush()

        task = Task(
            job_id=job.id,
            room_id="test-room",
            status=TaskStatus.RUNNING,
            started_at=datetime.now(timezone.utc) - timedelta(minutes=10),
        )
        session.add(task)
        await session.commit()

    async with async_session_factory() as session:
        count, _ = await cleanup_stuck_internal_tasks(
            session, timeout=timedelta(hours=1)
        )
        assert count == 0


@pytest.mark.asyncio
async def test_cleanup_stuck_skips_external_tasks(async_session_factory):
    """External (@global) RUNNING tasks are NOT cleaned up by this function."""
    async with async_session_factory() as session:
        job = Job(room_id="@global", category="modifiers", name="Rotate", schema_={})
        session.add(job)
        await session.flush()

        task = Task(
            job_id=job.id,
            room_id="test-room",
            status=TaskStatus.RUNNING,
            started_at=datetime.now(timezone.utc) - timedelta(hours=2),
        )
        session.add(task)
        await session.commit()

    async with async_session_factory() as session:
        count, _ = await cleanup_stuck_internal_tasks(
            session, timeout=timedelta(hours=1)
        )
        assert count == 0


@pytest.mark.asyncio
async def test_cleanup_worker_returns_task_status_emissions(
    async_session_factory, test_user_id
):
    """cleanup_worker should return TaskStatusEvent emissions for failed tasks."""
    worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        worker = Worker(
            id=worker_id,
            user_id=test_user_id,
            last_heartbeat=datetime.now(timezone.utc),
        )
        session.add(worker)
        job = Job(room_id="@global", category="modifiers", name="EmitTest", schema_={})
        session.add(job)
        await session.flush()
        link = WorkerJobLink(worker_id=worker_id, job_id=job.id)
        session.add(link)
        task = Task(
            job_id=job.id,
            room_id="room1",
            status=TaskStatus.RUNNING,
            worker_id=worker_id,
        )
        session.add(task)
        await session.commit()
        task_id = task.id

    async with async_session_factory() as session:
        result = await session.execute(select(Worker).where(Worker.id == worker_id))
        worker = result.scalar_one()
        emissions, _frame_rooms = await cleanup_worker(session, worker)
        await session.commit()

    # Should have TaskStatusEvent for the failed task + JobsInvalidate for orphan job
    task_events = [e for e in emissions if isinstance(e.event, TaskStatusEvent)]
    job_events = [e for e in emissions if isinstance(e.event, JobsInvalidate)]
    assert len(task_events) == 1
    assert task_events[0].event.id == str(task_id)
    assert task_events[0].event.status == "failed"
    assert task_events[0].room == "room:room1"
    assert len(job_events) == 1
    assert job_events[0].room == "room:@global"


@pytest.mark.asyncio
async def test_cleanup_stale_workers_returns_emissions(
    async_session_factory, test_user_id
):
    """cleanup_stale_workers should return count and emissions."""
    worker_id = uuid.uuid4()

    async with async_session_factory() as session:
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        worker = Worker(id=worker_id, user_id=test_user_id, last_heartbeat=stale_time)
        session.add(worker)
        await session.commit()

    async with async_session_factory() as session:
        count, emissions, frame_rooms = await cleanup_stale_workers(
            session, timedelta(seconds=60)
        )
        assert count == 1
        assert isinstance(emissions, set)
        assert isinstance(frame_rooms, set)


@pytest.mark.asyncio
async def test_cleanup_stuck_internal_returns_emissions(async_session_factory):
    """cleanup_stuck_internal_tasks should return count and emissions."""
    async with async_session_factory() as session:
        job = Job(
            room_id="@internal", category="modifiers", name="EmitInternal", schema_={}
        )
        session.add(job)
        await session.flush()
        task = Task(
            job_id=job.id,
            room_id="test-room",
            status=TaskStatus.RUNNING,
            started_at=datetime.now(timezone.utc) - timedelta(hours=2),
        )
        session.add(task)
        await session.commit()
        task_id = task.id

    async with async_session_factory() as session:
        count, emissions = await cleanup_stuck_internal_tasks(
            session, timeout=timedelta(hours=1)
        )
        assert count == 1
        task_events = [e for e in emissions if isinstance(e.event, TaskStatusEvent)]
        assert len(task_events) == 1
        assert task_events[0].event.id == str(task_id)
        assert task_events[0].event.status == "failed"
        assert task_events[0].room == "room:test-room"


@pytest.mark.asyncio
async def test_run_sweeper_cleans_up_stale_workers(async_session_factory, test_user_id):
    """run_sweeper should periodically clean up stale workers in a background task."""
    stale_worker_id = uuid.uuid4()

    # Create a stale worker (heartbeat 2 minutes ago)
    async with async_session_factory() as session:
        stale_time = datetime.now(timezone.utc) - timedelta(minutes=2)
        worker = Worker(
            id=stale_worker_id, user_id=test_user_id, last_heartbeat=stale_time
        )
        session.add(worker)
        await session.commit()

    # Create custom settings with very short intervals
    settings = JobLibSettings(
        sweeper_interval_seconds=1,  # 1 second interval
        worker_timeout_seconds=1,  # 1 second timeout
    )

    # Create session generator
    async def get_session():
        async with async_session_factory() as session:
            yield session

    # Start the sweeper as a background task
    sweeper_task = asyncio.create_task(run_sweeper(get_session, settings, None))

    try:
        # Wait for one sweep cycle (give it time to run at least once)
        await asyncio.sleep(1.5)

        # Verify the stale worker was cleaned up
        async with async_session_factory() as session:
            result = await session.execute(
                select(Worker).where(Worker.id == stale_worker_id)
            )
            worker = result.scalar_one_or_none()
            assert worker is None, "Stale worker should have been cleaned up"

    finally:
        # Cancel the sweeper task
        sweeper_task.cancel()
        try:
            await sweeper_task
        except asyncio.CancelledError:
            pass  # Expected when cancelling


@pytest.mark.asyncio
async def test_cleanup_stuck_internal_tasks_includes_claimed(async_session_factory):
    """CLAIMED @internal tasks with old created_at are marked FAILED after timeout."""
    async with async_session_factory() as session:
        job = Job(room_id="@internal", category="selections", name="Stuck", schema_={})
        session.add(job)
        await session.flush()

        task = Task(
            job_id=job.id,
            room_id="test-room",
            status=TaskStatus.CLAIMED,
            started_at=None,
            created_at=datetime.now(timezone.utc) - timedelta(hours=2),
        )
        session.add(task)
        await session.commit()
        task_id = task.id

    async with async_session_factory() as session:
        count, _ = await cleanup_stuck_internal_tasks(
            session, timeout=timedelta(hours=1)
        )
        assert count == 1

    async with async_session_factory() as session:
        result = await session.execute(select(Task).where(Task.id == task_id))
        task = result.scalar_one()
        assert task.status == TaskStatus.FAILED
        assert task.error == "Internal worker timeout"
        assert task.completed_at is not None
