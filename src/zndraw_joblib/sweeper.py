# src/zndraw_joblib/sweeper.py
"""Background sweeper for cleaning up stale workers."""

import asyncio
import logging
import uuid
from collections.abc import AsyncGenerator, Callable
from datetime import UTC, datetime, timedelta
from typing import Any

from sqlalchemy import func as sa_func
from sqlalchemy.orm import selectinload
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession
from zndraw_socketio import AsyncServerWrapper

from zndraw_joblib.events import (
    Emission,
    JobsInvalidate,
    ProvidersInvalidate,
    build_task_status_emission,
    emit,
)
from zndraw_joblib.models import (
    Job,
    ProviderRecord,
    Task,
    TaskStatus,
    Worker,
    WorkerJobLink,
)
from zndraw_joblib.settings import JobLibSettings

logger = logging.getLogger(__name__)


async def _soft_delete_orphan_job(
    session: AsyncSession, job_id: uuid.UUID
) -> set[Emission]:
    """Soft-delete a job with no workers and no pending tasks.

    Only soft-deletes if:
    - Job has no workers
    - Job has no pending tasks (new workers could register and pick them up)

    Skips @internal jobs since they are server-managed and always available.
    Note: Does NOT commit the transaction - caller must commit.
    """
    # @internal jobs are never orphaned — they're registered at startup
    result = await session.exec(select(Job).where(Job.id == job_id))
    job = result.one_or_none()
    if not job or job.room_id == "@internal":
        return set()

    # Check if job has any remaining workers
    result = await session.exec(
        select(WorkerJobLink).where(WorkerJobLink.job_id == job_id).limit(1)
    )
    if result.one_or_none():
        return set()  # Job still has workers

    # Check if job has any pending tasks (new workers could register and pick them up)
    result = await session.exec(
        select(Task)
        .where(
            Task.job_id == job_id,
            Task.status == TaskStatus.PENDING,
        )
        .limit(1)
    )
    if result.one_or_none():
        return set()  # Job has pending tasks, keep it alive

    emissions: set[Emission] = set()

    # Soft-delete the orphaned job (no workers, no pending tasks)
    job.deleted = True
    session.add(job)
    emissions.add(Emission(JobsInvalidate(), f"room:{job.room_id}"))
    return emissions


async def cleanup_worker(
    session: AsyncSession, worker: Worker
) -> tuple[set[Emission], set[str]]:
    """Clean up a worker: fail tasks, remove links, soft-delete orphan jobs.

    This is the shared cleanup logic used by both delete_worker endpoint and sweeper.
    Note: Does NOT commit the transaction - caller must commit.

    Returns
    -------
    tuple[set[Emission], set[str]]
        A tuple of (emissions, room IDs that had frame providers removed).
        Frame provider room IDs let callers clean up external state (e.g. Redis keys).
    """
    emissions: set[Emission] = set()
    now = datetime.now(UTC)

    # Fail any claimed/running tasks owned by this worker
    result = await session.exec(
        select(Task)
        .options(selectinload(Task.job))
        .where(
            Task.worker_id == worker.id,
            Task.status.in_({TaskStatus.CLAIMED, TaskStatus.RUNNING}),
        )
    )
    worker_tasks = result.all()
    for task in worker_tasks:
        task.status = TaskStatus.FAILED
        task.completed_at = now
        task.error = "Worker disconnected"
        session.add(task)
        emissions.add(
            build_task_status_emission(task, task.job.full_name if task.job else "")
        )

    # Get links this worker has (need both job_ids and the link objects)
    result = await session.exec(
        select(WorkerJobLink).where(WorkerJobLink.worker_id == worker.id)
    )
    links = result.all()
    job_ids = [link.job_id for link in links]

    # Fetch room_ids for affected jobs (worker count is changing)
    if job_ids:
        result = await session.exec(
            select(Job.id, Job.room_id).where(Job.id.in_(job_ids))
        )
        job_rooms = {row.id: row.room_id for row in result.all()}
    else:
        job_rooms = {}

    # Emit JobsInvalidate for all affected rooms (worker count changed)
    for room_id in set(job_rooms.values()):
        emissions.add(Emission(JobsInvalidate(), f"room:{room_id}"))

    # Delete providers owned by this worker
    result = await session.exec(
        select(ProviderRecord).where(ProviderRecord.worker_id == worker.id)
    )
    providers = result.all()
    provider_rooms: set[str] = set()
    frame_provider_rooms: set[str] = set()
    for provider in providers:
        if provider.category == "frames":
            frame_provider_rooms.add(provider.room_id)
        provider_rooms.add(provider.room_id)
        await session.delete(provider)
    for room_id in provider_rooms:
        emissions.add(Emission(ProvidersInvalidate(), f"room:{room_id}"))

    # Delete all links
    for link in links:
        await session.delete(link)

    # Delete worker
    await session.delete(worker)
    await session.flush()

    # Clean up orphan jobs (no workers and no non-terminal tasks)
    for job_id in job_ids:
        emissions |= await _soft_delete_orphan_job(session, job_id)

    return emissions, frame_provider_rooms


async def cleanup_stale_workers(
    session: AsyncSession, stale_after: timedelta
) -> tuple[int, set[Emission], set[str]]:
    """Find and clean up workers with stale heartbeats.

    Parameters
    ----------
    session : AsyncSession
        Async database session.
    stale_after : timedelta
        How long since last heartbeat before a worker is considered stale.

    Returns
    -------
    tuple[int, set[Emission], set[str]]
        Tuple of (count of workers cleaned up, emissions, room IDs that had
        frame providers removed).
    """
    cutoff = datetime.now(UTC) - stale_after
    all_emissions: set[Emission] = set()
    all_frame_rooms: set[str] = set()

    # Find all stale workers. @internal providers mirror @internal jobs —
    # jobs have no WorkerJobLink entries, providers have worker_id=None —
    # because both are server-owned (dispatched by the in-process taskiq
    # broker, not a remote client). NULL sits outside cleanup_worker's
    # concrete-UUID delete query, so no sweeper special case is needed.
    result = await session.exec(select(Worker).where(Worker.last_heartbeat < cutoff))
    stale_workers = result.all()

    count = 0
    for worker in stale_workers:
        logger.info("Cleaning up stale worker: %s", worker.id)
        emissions, frame_rooms = await cleanup_worker(session, worker)
        all_emissions |= emissions
        all_frame_rooms |= frame_rooms
        count += 1

    if count > 0:
        await session.commit()

    return count, all_emissions, all_frame_rooms


async def cleanup_stuck_internal_tasks(
    session: AsyncSession, stuck_after: timedelta
) -> tuple[int, set[Emission]]:
    """Find and fail internal tasks stuck in RUNNING or CLAIMED beyond the time limit.

    Parameters
    ----------
    session : AsyncSession
        Async database session.
    stuck_after : timedelta
        How long a task can be in RUNNING/CLAIMED before being considered stuck.

    Returns
    -------
    tuple[int, set[Emission]]
        Tuple of (count of tasks failed, set of emissions).
    """
    cutoff = datetime.now(UTC) - stuck_after
    now = datetime.now(UTC)
    emissions: set[Emission] = set()

    result = await session.exec(
        select(Task)
        .join(Job)
        .options(selectinload(Task.job))
        .where(
            Job.room_id == "@internal",
            Task.status.in_({TaskStatus.RUNNING, TaskStatus.CLAIMED}),
            sa_func.coalesce(Task.started_at, Task.created_at) < cutoff,
        )
    )
    stuck_tasks = result.all()

    count = 0
    for task in stuck_tasks:
        task.status = TaskStatus.FAILED
        task.completed_at = now
        task.error = "Internal worker timeout"
        session.add(task)
        emissions.add(
            build_task_status_emission(task, task.job.full_name if task.job else "")
        )
        count += 1

    if count > 0:
        await session.commit()
        logger.info("Failed %d stuck internal task(s)", count)

    return count, emissions


async def run_sweeper(
    get_session: Callable[[], AsyncGenerator[AsyncSession, None]],
    settings: JobLibSettings,
    tsio: AsyncServerWrapper | None = None,
    on_frame_rooms: Callable[[set[str]], Any] | None = None,
) -> None:
    """Background task that runs cleanup periodically.

    Parameters
    ----------
    get_session
        Async generator yielding database sessions.
    settings
        JobLib configuration.
    tsio
        Optional Socket.IO server for broadcasting events.
    on_frame_rooms
        Optional async callback invoked with room IDs that had frame providers
        removed. Host apps use this to clean up external state (e.g. Redis keys).
    """
    timeout = timedelta(seconds=settings.worker_timeout_seconds)
    internal_timeout = timedelta(seconds=settings.internal_task_timeout_seconds)
    interval = settings.sweeper_interval_seconds

    logger.info(
        "Starting sweeper: interval=%ss, worker_timeout=%ss, internal_task_timeout=%ss",
        interval,
        settings.worker_timeout_seconds,
        settings.internal_task_timeout_seconds,
    )

    while True:
        await asyncio.sleep(interval)
        try:
            async for session in get_session():
                count, emissions, frame_rooms = await cleanup_stale_workers(
                    session, timeout
                )
                if count > 0:
                    logger.info("Cleaned up %s stale worker(s)", count)
                await emit(tsio, emissions)
                if frame_rooms and on_frame_rooms is not None:
                    await on_frame_rooms(frame_rooms)

            async for session in get_session():
                count, emissions = await cleanup_stuck_internal_tasks(
                    session, internal_timeout
                )
                if count > 0:
                    logger.info("Failed %s stuck internal task(s)", count)
                await emit(tsio, emissions)
        except Exception:
            logger.exception("Error in sweeper")
