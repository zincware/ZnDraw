# src/zndraw_joblib/router.py
import asyncio
import json
import logging
import random
import re
from datetime import UTC, datetime
from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends, Header, Query, Request, Response, status
from sqlalchemy import and_, func, update
from sqlalchemy.exc import OperationalError
from sqlalchemy.ext.asyncio import async_sessionmaker
from sqlalchemy.orm import selectinload
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession
from zndraw_socketio import AsyncServerWrapper

from zndraw_auth import (
    User,
    current_active_user,
    current_superuser,
    current_user_scoped_session,
)
from zndraw_auth.db import SessionDep, get_session_maker
from zndraw_joblib.dependencies import (
    FrameRoomCleanupDep,
    JobLibSettingsDep,
    ResultBackendDep,
    WorkerTokenDep,
    WritableRoomDep,
    get_internal_provider_registry,
    get_internal_registry,
    get_tsio,
    mint_internal_worker_token,
    request_hash,
    validate_room_id,
)
from zndraw_joblib.events import (
    Emission,
    JobsInvalidate,
    ProviderRequest,
    ProviderResultReady,
    ProvidersInvalidate,
    TaskAvailable,
    build_task_status_emission,
    emit,
)
from zndraw_joblib.exceptions import (
    Forbidden,
    InternalJobNotConfigured,
    InvalidCategory,
    InvalidTaskTransition,
    JobNotFound,
    NoWorkersAvailable,
    ProblemDetail,
    ProviderNotFound,
    ProviderTimeout,
    SchemaConflict,
    TaskNotFound,
    WorkerNotFound,
)
from zndraw_joblib.models import (
    TERMINAL_STATUSES,
    Job,
    ProviderRecord,
    Task,
    TaskStatus,
    Worker,
    WorkerJobLink,
)
from zndraw_joblib.registry import InternalProviderRegistry, InternalRegistry
from zndraw_joblib.schemas import (
    JobRegisterRequest,
    JobResponse,
    JobSummary,
    PaginatedResponse,
    ProviderRegisterRequest,
    ProviderResponse,
    TaskClaimRequest,
    TaskClaimResponse,
    TaskResponse,
    TaskSubmitRequest,
    TaskUpdateRequest,
    WorkerResponse,
    WorkerSummary,
)
from zndraw_joblib.sweeper import _soft_delete_orphan_job, cleanup_worker

logger = logging.getLogger(__name__)
_rng = random.SystemRandom()

# Type aliases for dependency injection
CurrentUserDep = Annotated[User, Depends(current_active_user)]
CurrentUserFactoryDep = Annotated[User, Depends(current_user_scoped_session)]
SuperUserDep = Annotated[User, Depends(current_superuser)]
SettingsDep = JobLibSettingsDep
SessionMakerDep = Annotated[
    async_sessionmaker[AsyncSession], Depends(get_session_maker)
]
InternalRegistryDep = Annotated[InternalRegistry | None, Depends(get_internal_registry)]
InternalProviderRegistryDep = Annotated[
    InternalProviderRegistry | None, Depends(get_internal_provider_registry)
]
TsioDep = Annotated[AsyncServerWrapper | None, Depends(get_tsio)]

# Valid status transitions
VALID_TRANSITIONS: dict[TaskStatus, set[TaskStatus]] = {
    TaskStatus.PENDING: {TaskStatus.CLAIMED, TaskStatus.CANCELLED},
    TaskStatus.CLAIMED: {TaskStatus.RUNNING, TaskStatus.FAILED, TaskStatus.CANCELLED},
    TaskStatus.RUNNING: {TaskStatus.COMPLETED, TaskStatus.FAILED, TaskStatus.CANCELLED},
    TaskStatus.COMPLETED: set(),
    TaskStatus.FAILED: set(),
    TaskStatus.CANCELLED: set(),
}


def parse_prefer_wait(prefer_header: str | None) -> int | None:
    """
    Parse RFC 7240 Prefer header for wait directive.
    Returns seconds to wait, or None if not specified.
    """
    if not prefer_header:
        return None
    match = re.search(r"\bwait=(\d+)\b", prefer_header)
    if match:
        return int(match.group(1))
    return None


async def _queue_position(session: AsyncSession, task: Task) -> int | None:
    """Get queue position for a single task via the bulk helper."""
    if task.status != TaskStatus.PENDING:
        return None
    positions = await _bulk_queue_positions(session, [task.id])
    return positions.get(task.id)


async def _resolve_job(session: AsyncSession, job_name: str) -> Job:
    parts = job_name.split(":", 2)
    if len(parts) != 3:
        raise JobNotFound.exception(detail=f"Invalid job name format: {job_name}")
    job_room_id, category, name = parts

    result = await session.exec(
        select(Job).where(
            Job.room_id == job_room_id,
            Job.category == category,
            Job.name == name,
        )
    )
    job = result.one_or_none()
    if not job or job.deleted:
        raise JobNotFound.exception(detail=f"Job '{job_name}' not found")
    return job


async def _task_response(session: AsyncSession, task: Task) -> TaskResponse:
    result = await session.exec(select(Job).where(Job.id == task.job_id))
    job = result.one_or_none()
    return TaskResponse(
        id=task.id,
        job_name=job.full_name if job else "",
        room_id=task.room_id,
        status=task.status,
        created_at=task.created_at,
        started_at=task.started_at,
        completed_at=task.completed_at,
        worker_id=task.worker_id,
        error=task.error,
        payload=task.payload,
        queue_position=await _queue_position(session, task),
    )


async def _bulk_queue_positions(
    session: AsyncSession, task_ids: list[UUID]
) -> dict[UUID, int]:
    """Compute queue positions for pending tasks in a single query.

    Uses a window function for efficiency.
    """
    if not task_ids:
        return {}
    # Window function: rank each pending task within its job by created_at
    row_num = (
        func.row_number()
        .over(partition_by=Task.job_id, order_by=Task.created_at.asc())
        .label("pos")
    )
    subq = select(Task.id, row_num).where(Task.status == TaskStatus.PENDING).subquery()
    result = await session.exec(
        select(subq.c.id, subq.c.pos).where(subq.c.id.in_(task_ids))
    )
    return {row.id: row.pos for row in result}


async def _bulk_task_responses(
    session: AsyncSession, tasks: list[Task]
) -> list[TaskResponse]:
    """Build TaskResponse list efficiently by batching job lookups and queue positions.

    Expects tasks to have been loaded with selectinload(Task.job).
    """
    pending_ids = [t.id for t in tasks if t.status == TaskStatus.PENDING]
    queue_positions = await _bulk_queue_positions(session, pending_ids)

    return [
        TaskResponse(
            id=t.id,
            job_name=t.job.full_name if t.job else "",
            room_id=t.room_id,
            status=t.status,
            created_at=t.created_at,
            started_at=t.started_at,
            completed_at=t.completed_at,
            worker_id=t.worker_id,
            error=t.error,
            payload=t.payload,
            queue_position=queue_positions.get(t.id),
        )
        for t in tasks
    ]


async def _task_status_emission(session: AsyncSession, task: Task) -> Emission:
    """Build a TaskStatusEvent emission from a task.

    Queries job name and queue position.
    """
    result = await session.exec(select(Job).where(Job.id == task.job_id))
    job = result.one_or_none()
    return build_task_status_emission(
        task,
        job_full_name=job.full_name if job else "",
        queue_position=await _queue_position(session, task),
    )


def _room_job_filter(room_id: str):
    """Build a SQLAlchemy filter for jobs visible from a given room."""
    if room_id == "@global":
        return Job.room_id == "@global"
    if room_id == "@internal":
        return Job.room_id == "@internal"
    return (Job.room_id.in_(["@global", "@internal"])) | (Job.room_id == room_id)


router = APIRouter(prefix="/v1/joblib", tags=["joblib"])


@router.post(
    "/workers", response_model=WorkerResponse, status_code=status.HTTP_201_CREATED
)
async def create_worker(
    session: SessionDep,
    user: CurrentUserDep,
):
    """Create a new worker for the authenticated user.

    Returns the worker_id which the client should use for heartbeats
    and store locally for future requests.
    """
    worker = Worker(user_id=user.id)
    session.add(worker)
    await session.commit()
    await session.refresh(worker)
    return WorkerResponse(id=worker.id, last_heartbeat=worker.last_heartbeat)


@router.put(
    "/rooms/{room_id}/jobs",
    response_model=JobResponse,
    status_code=status.HTTP_201_CREATED,
)
async def register_job(
    room_id: WritableRoomDep,
    request: JobRegisterRequest,
    response: Response,
    session: SessionDep,
    user: CurrentUserDep,
    settings: SettingsDep,
    tsio: TsioDep,
):
    """Register a job for a room. Creates worker and link if not exists."""
    # Check admin for @global and @internal
    if room_id in ("@global", "@internal") and not user.is_superuser:
        raise Forbidden.exception(
            detail="Admin required for @global/@internal job registration"
        )

    # Validate category
    if request.category not in settings.allowed_categories:
        raise InvalidCategory.exception(
            detail=(
                f"Category '{request.category}' not in allowed list: "
                f"{settings.allowed_categories}"
            )
        )

    # Check if job exists
    result = await session.exec(
        select(Job).where(
            Job.room_id == room_id,
            Job.category == request.category,
            Job.name == request.name,
        )
    )
    existing_job = result.one_or_none()

    if existing_job and existing_job.deleted:
        # Re-activate soft-deleted job with new schema
        existing_job.deleted = False
        existing_job.schema_ = request.schema_
        job = existing_job
    elif existing_job:
        # Validate schema match
        if existing_job.schema_ != request.schema_:
            raise SchemaConflict.exception(
                detail=f"Schema mismatch for job '{existing_job.full_name}'"
            )
        # Idempotent - ensure worker link exists
        job = existing_job
        response.status_code = status.HTTP_200_OK
    else:
        # Create new job
        job = Job(
            room_id=room_id,
            category=request.category,
            name=request.name,
            schema_=request.schema_,
        )
        session.add(job)
        await session.flush()

    # Handle worker: use provided worker_id or auto-create
    if request.worker_id:
        # Verify ownership
        result = await session.exec(
            select(Worker).where(
                Worker.id == request.worker_id,
                Worker.user_id == user.id,
            )
        )
        worker = result.one_or_none()
        if not worker:
            raise WorkerNotFound.exception(
                detail=f"Worker '{request.worker_id}' not found or not owned by user"
            )
    else:
        # Auto-create worker for this user
        worker = Worker(user_id=user.id)
        session.add(worker)
        await session.flush()

    # Ensure worker-job link exists
    result = await session.exec(
        select(WorkerJobLink).where(
            WorkerJobLink.worker_id == worker.id,
            WorkerJobLink.job_id == job.id,
        )
    )
    link = result.one_or_none()
    if not link:
        link = WorkerJobLink(worker_id=worker.id, job_id=job.id)
        session.add(link)

    await session.commit()
    await emit(tsio, {Emission(JobsInvalidate(), f"room:{room_id}")})
    await session.refresh(job)

    # Get worker IDs for this job
    result = await session.exec(
        select(WorkerJobLink).where(WorkerJobLink.job_id == job.id)
    )
    worker_links = result.all()
    worker_ids = [link.worker_id for link in worker_links]

    return JobResponse(
        id=job.id,
        room_id=job.room_id,
        category=job.category,
        name=job.name,
        full_name=job.full_name,
        schema=job.schema_,
        workers=worker_ids,
        worker_id=worker.id,
    )


@router.get("/rooms/{room_id}/jobs", response_model=PaginatedResponse[JobSummary])
async def list_jobs(
    room_id: str,
    session: SessionDep,
    limit: Annotated[int, Query(ge=0, le=500)] = 50,
    offset: Annotated[int, Query(ge=0)] = 0,
):
    """List jobs for a room. Includes @global jobs unless room_id is @global."""
    validate_room_id(room_id)

    base_query = select(Job).where(_room_job_filter(room_id), Job.deleted.is_(False))

    # Total count
    total_result = await session.exec(
        select(func.count()).select_from(base_query.subquery())
    )
    total = total_result.one()

    # Paginated + eager-load workers
    result = await session.exec(
        base_query.options(selectinload(Job.workers))
        .order_by(Job.created_at.desc())
        .offset(offset)
        .limit(limit)
    )
    jobs = result.all()

    items = [
        JobSummary(
            full_name=job.full_name,
            category=job.category,
            name=job.name,
            workers=[w.id for w in job.workers],
        )
        for job in jobs
    ]
    return PaginatedResponse(items=items, total=total, limit=limit, offset=offset)


@router.get("/rooms/{room_id}/workers", response_model=PaginatedResponse[WorkerSummary])
async def list_workers_for_room(
    room_id: str,
    session: SessionDep,
    limit: Annotated[int, Query(ge=0, le=500)] = 50,
    offset: Annotated[int, Query(ge=0)] = 0,
):
    """List workers for a room. Includes @global workers unless room_id is @global."""
    validate_room_id(room_id)

    result = await session.exec(
        select(Job.id).where(_room_job_filter(room_id), Job.deleted.is_(False))
    )
    job_ids = result.all()

    if not job_ids:
        return PaginatedResponse(items=[], total=0, limit=limit, offset=offset)

    # Find distinct worker IDs linked to these jobs
    worker_id_query = (
        select(WorkerJobLink.worker_id)
        .where(WorkerJobLink.job_id.in_(job_ids))
        .distinct()
    )

    # Total count
    total_result = await session.exec(
        select(func.count()).select_from(worker_id_query.subquery())
    )
    total = total_result.one()

    # Paginated workers with eager-loaded jobs
    result = await session.exec(
        select(Worker)
        .where(Worker.id.in_(worker_id_query))
        .options(selectinload(Worker.jobs))
        .order_by(Worker.created_at.desc())
        .offset(offset)
        .limit(limit)
    )
    workers = result.all()

    items = [
        WorkerSummary(
            id=worker.id,
            last_heartbeat=worker.last_heartbeat,
            job_count=len(worker.jobs),
        )
        for worker in workers
    ]
    return PaginatedResponse(items=items, total=total, limit=limit, offset=offset)


@router.get("/rooms/{room_id}/tasks", response_model=PaginatedResponse[TaskResponse])
async def list_tasks_for_room(
    room_id: str,
    session: SessionDep,
    task_status: Annotated[TaskStatus | None, Query(alias="status")] = None,
    limit: Annotated[int, Query(ge=0, le=500)] = 50,
    offset: Annotated[int, Query(ge=0)] = 0,
):
    """List tasks for a room, optionally filtered by status.

    Includes queue position for pending tasks.
    """
    validate_room_id(room_id)

    base_query = select(Task).where(Task.room_id == room_id)
    if task_status:
        base_query = base_query.where(Task.status == task_status)

    # Total count
    total_result = await session.exec(
        select(func.count()).select_from(base_query.subquery())
    )
    total = total_result.one()

    # Paginated + eager-load job relationship
    result = await session.exec(
        base_query.options(selectinload(Task.job))
        .order_by(Task.created_at.desc())
        .offset(offset)
        .limit(limit)
    )
    tasks = result.all()

    items = await _bulk_task_responses(session, tasks)
    return PaginatedResponse(items=items, total=total, limit=limit, offset=offset)


@router.get(
    "/rooms/{room_id}/jobs/{job_name:path}/tasks",
    response_model=PaginatedResponse[TaskResponse],
)
async def list_tasks_for_job(
    room_id: str,
    job_name: str,
    session: SessionDep,
    task_status: Annotated[TaskStatus | None, Query(alias="status")] = None,
    limit: Annotated[int, Query(ge=0, le=500)] = 50,
    offset: Annotated[int, Query(ge=0)] = 0,
):
    """List tasks for a specific job. Includes queue position for pending tasks."""
    validate_room_id(room_id)

    job = await _resolve_job(session, job_name)

    base_query = select(Task).where(Task.job_id == job.id, Task.room_id == room_id)
    if task_status:
        base_query = base_query.where(Task.status == task_status)

    # Total count
    total_result = await session.exec(
        select(func.count()).select_from(base_query.subquery())
    )
    total = total_result.one()

    # Paginated + eager-load job relationship
    result = await session.exec(
        base_query.options(selectinload(Task.job))
        .order_by(Task.created_at.desc())
        .offset(offset)
        .limit(limit)
    )
    tasks = result.all()

    items = await _bulk_task_responses(session, tasks)
    return PaginatedResponse(items=items, total=total, limit=limit, offset=offset)


@router.get("/rooms/{room_id}/jobs/{job_name:path}", response_model=JobResponse)
async def get_job(
    room_id: str,
    job_name: str,
    session: SessionDep,
):
    """Get job details by full name."""
    validate_room_id(room_id)

    job = await _resolve_job(session, job_name)

    result = await session.exec(
        select(WorkerJobLink).where(WorkerJobLink.job_id == job.id)
    )
    worker_links = result.all()
    worker_ids = [link.worker_id for link in worker_links]

    return JobResponse(
        id=job.id,
        room_id=job.room_id,
        category=job.category,
        name=job.name,
        full_name=job.full_name,
        schema=job.schema_,
        workers=worker_ids,
    )


@router.post(
    "/rooms/{room_id}/tasks/{job_name:path}",
    response_model=TaskResponse,
    status_code=status.HTTP_202_ACCEPTED,
)
async def submit_task(
    room_id: WritableRoomDep,
    job_name: str,
    request: TaskSubmitRequest,
    response: Response,
    session: SessionDep,
    user: CurrentUserDep,
    internal_registry: InternalRegistryDep,
    tsio: TsioDep,
    worker_token: WorkerTokenDep,
):
    """Submit a task for processing."""
    job = await _resolve_job(session, job_name)

    # Validate internal registry BEFORE creating the task to avoid orphans
    if job.room_id == "@internal" and (
        internal_registry is None or job.full_name not in internal_registry.tasks
    ):
        raise InternalJobNotConfigured.exception(
            detail=(
                f"Internal job '{job.full_name}' is registered"
                " but no executor is available"
            )
        )

    # Reject submission for non-@internal jobs with no connected workers
    if job.room_id != "@internal":
        worker_count = (
            await session.exec(
                select(func.count()).where(WorkerJobLink.job_id == job.id)
            )
        ).one()
        if worker_count == 0:
            raise NoWorkersAvailable.exception(
                detail=f"Job '{job.full_name}' has no connected workers"
            )

    # Create task
    task = Task(
        job_id=job.id,
        room_id=room_id,
        created_by_id=user.id,
        payload=request.payload,
        status=TaskStatus.PENDING,
    )
    session.add(task)
    await session.commit()
    await session.refresh(task)

    # Dispatch to taskiq for @internal jobs (fail task if dispatch fails)
    if job.room_id == "@internal":
        try:
            await internal_registry.tasks[job.full_name].kiq(
                task_id=str(task.id),
                room_id=room_id,
                payload=request.payload,
                token=worker_token,
            )
        except Exception:  # noqa: BLE001
            task.status = TaskStatus.FAILED
            task.completed_at = datetime.now(UTC)
            task.error = "Failed to dispatch to internal executor"
            await session.commit()
            await session.refresh(task)
            return await _task_response(session, task)
        # Server implicitly claims the task for internal execution
        task.status = TaskStatus.CLAIMED
        await session.commit()
        await session.refresh(task)

    # Emit events (after successful dispatch for @internal)
    emissions: set[Emission] = {await _task_status_emission(session, task)}
    if job.room_id != "@internal":
        emissions.add(
            Emission(
                TaskAvailable(
                    job_name=job.full_name,
                    room_id=room_id,
                    task_id=str(task.id),
                ),
                f"jobs:{job.full_name}",
            )
        )
    await emit(tsio, emissions)

    # Set Location header
    response.headers["Location"] = f"/v1/joblib/tasks/{task.id}"
    response.headers["Retry-After"] = "1"

    return await _task_response(session, task)


@router.post("/tasks/claim", response_model=TaskClaimResponse)
async def claim_task(
    request: TaskClaimRequest,
    session: SessionDep,
    user: CurrentUserDep,
    settings: SettingsDep,
    tsio: TsioDep,
):
    """Claim the oldest pending task for jobs the specified worker is registered for."""
    # Validate that worker_id exists and belongs to the authenticated user
    result = await session.exec(select(Worker).where(Worker.id == request.worker_id))
    worker = result.one_or_none()

    if not worker:
        raise WorkerNotFound.exception(f"Worker {request.worker_id} not found")

    if worker.user_id != user.id:
        raise Forbidden.exception("Worker belongs to a different user")

    # Find job IDs for this specific worker
    result = await session.exec(
        select(WorkerJobLink.job_id).where(WorkerJobLink.worker_id == request.worker_id)
    )
    worker_job_ids = result.all()

    if not worker_job_ids:
        return TaskClaimResponse(task=None)

    # Use atomic UPDATE with WHERE clause for optimistic locking
    # This handles the race condition where multiple workers try to claim the same task
    max_attempts = settings.claim_max_attempts
    base_delay = settings.claim_base_delay_seconds
    claimed_task_id = None

    for attempt in range(max_attempts):
        try:
            # Find oldest pending task for jobs this worker is registered for
            result = await session.exec(
                select(Task.id)
                .where(
                    Task.job_id.in_(worker_job_ids), Task.status == TaskStatus.PENDING
                )
                .order_by(Task.created_at.asc())
                .limit(1)
            )
            task_id = result.one_or_none()

            if not task_id:
                return TaskClaimResponse(task=None)

            # Atomically update only if still PENDING (optimistic locking)
            stmt = (
                update(Task)
                .where(Task.id == task_id, Task.status == TaskStatus.PENDING)
                .values(status=TaskStatus.CLAIMED, worker_id=request.worker_id)
            )
            cursor_result = await session.exec(stmt)
            await session.commit()

            # Check if we actually claimed it (rowcount == 1 means success)
            if cursor_result.rowcount == 1:
                claimed_task_id = task_id
                break

            # Another worker claimed it first - exponential backoff with jitter
            delay = base_delay * (2**attempt) * (0.5 + _rng.random())
            await asyncio.sleep(delay)

        except OperationalError:
            # Database locked/timeout - rollback and retry with backoff
            logger.warning(
                "OperationalError during claim (attempt %d/%d)",
                attempt + 1,
                max_attempts,
            )
            await session.rollback()
            delay = base_delay * (2**attempt) * (0.5 + _rng.random())
            await asyncio.sleep(delay)

    if not claimed_task_id:
        logger.warning(
            "Failed to claim task after %d attempts for worker %s",
            max_attempts,
            request.worker_id,
        )
        return TaskClaimResponse(task=None)

    # Fetch the claimed task
    result = await session.exec(select(Task).where(Task.id == claimed_task_id))
    task = result.one()

    await emit(tsio, {await _task_status_emission(session, task)})

    return TaskClaimResponse(task=await _task_response(session, task))


@router.get("/tasks/{task_id}", response_model=TaskResponse)
async def get_task_status(
    task_id: UUID,
    response: Response,
    session_maker: SessionMakerDep,
    settings: SettingsDep,
    prefer: Annotated[str | None, Header()] = None,
):
    """Get task status. Supports long-polling via Prefer: wait=N header."""
    # Initial lookup
    async with session_maker() as session:
        result = await session.exec(select(Task).where(Task.id == task_id))
        task = result.one_or_none()
        if not task:
            raise TaskNotFound.exception(detail=f"Task '{task_id}' not found")

    requested_wait = parse_prefer_wait(prefer)

    # Only long-poll if wait requested AND task not in terminal state
    if requested_wait and requested_wait > 0 and task.status not in TERMINAL_STATUSES:
        effective_wait = min(requested_wait, settings.long_poll_max_wait_seconds)
        elapsed = 0.0
        poll_interval = 1.0

        while elapsed < effective_wait and task.status not in TERMINAL_STATUSES:
            await asyncio.sleep(poll_interval)
            elapsed += poll_interval
            poll_interval = min(
                poll_interval * 1.5, 5.0
            )  # Exponential backoff, cap at 5s

            async with session_maker() as session:
                result = await session.exec(select(Task).where(Task.id == task_id))
                task = result.one_or_none()
                if not task:
                    raise TaskNotFound.exception(detail=f"Task '{task_id}' not found")

        response.headers["Preference-Applied"] = f"wait={int(effective_wait)}"

    # Build final response — re-fetch to avoid stale detached object
    async with session_maker() as session:
        result = await session.exec(select(Task).where(Task.id == task_id))
        task = result.one()
        return await _task_response(session, task)


@router.patch("/tasks/{task_id}", response_model=TaskResponse)
async def update_task_status(
    task_id: UUID,
    request: TaskUpdateRequest,
    session: SessionDep,
    user: CurrentUserDep,
    tsio: TsioDep,
):
    """Update task status. Requires the task's worker owner or superuser."""
    result = await session.exec(select(Task).where(Task.id == task_id))
    task = result.one_or_none()
    if not task:
        raise TaskNotFound.exception(detail=f"Task '{task_id}' not found")

    # Authorization: worker owner or superuser
    if not user.is_superuser:
        if task.worker_id is None:
            raise Forbidden.exception(detail="Task not claimed by any worker")
        result = await session.exec(select(Worker).where(Worker.id == task.worker_id))
        worker = result.one_or_none()
        if not worker or worker.user_id != user.id:
            raise Forbidden.exception(detail="Not authorized to update this task")

    # Validate transition
    if request.status not in VALID_TRANSITIONS.get(task.status, set()):
        raise InvalidTaskTransition.exception(
            detail=(
                f"Cannot transition from '{task.status.value}'"
                f" to '{request.status.value}'"
            )
        )

    # Update status
    task.status = request.status
    now = datetime.now(UTC)

    if request.status == TaskStatus.RUNNING:
        task.started_at = now
    elif request.status in TERMINAL_STATUSES:
        task.completed_at = now
        if request.error:
            task.error = request.error

    session.add(task)

    # Handle orphan job cleanup in the same transaction
    orphan_emissions: set[Emission] = set()
    if request.status in TERMINAL_STATUSES:
        await session.flush()
        orphan_emissions = await _soft_delete_orphan_job(session, task.job_id)

    await session.commit()
    await session.refresh(task)

    # Emit events after commit
    emissions: set[Emission] = {await _task_status_emission(session, task)}
    emissions |= orphan_emissions
    await emit(tsio, emissions)

    return await _task_response(session, task)


@router.get("/workers", response_model=PaginatedResponse[WorkerSummary])
async def list_workers(
    session: SessionDep,
    limit: Annotated[int, Query(ge=0, le=500)] = 50,
    offset: Annotated[int, Query(ge=0)] = 0,
):
    """List all workers with their job counts."""
    # Total count
    total_result = await session.exec(select(func.count()).select_from(Worker))
    total = total_result.one()

    # Paginated + eager-load jobs
    result = await session.exec(
        select(Worker)
        .options(selectinload(Worker.jobs))
        .order_by(Worker.created_at.desc())
        .offset(offset)
        .limit(limit)
    )
    workers = result.all()

    items = [
        WorkerSummary(
            id=worker.id,
            last_heartbeat=worker.last_heartbeat,
            job_count=len(worker.jobs),
        )
        for worker in workers
    ]
    return PaginatedResponse(items=items, total=total, limit=limit, offset=offset)


@router.patch("/workers/{worker_id}", response_model=WorkerResponse)
async def worker_heartbeat(
    worker_id: UUID,
    session: SessionDep,
    user: CurrentUserDep,
):
    """Update worker heartbeat. Worker must belong to authenticated user."""
    result = await session.exec(select(Worker).where(Worker.id == worker_id))
    worker = result.one_or_none()
    if not worker:
        raise WorkerNotFound.exception(detail=f"Worker '{worker_id}' not found")

    if worker.user_id != user.id:
        raise Forbidden.exception(detail="Worker belongs to different user")

    worker.last_heartbeat = datetime.now(UTC)
    session.add(worker)
    await session.commit()
    await session.refresh(worker)

    return WorkerResponse(id=worker.id, last_heartbeat=worker.last_heartbeat)


@router.delete("/workers/{worker_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_worker(
    worker_id: UUID,
    session: SessionDep,
    user: CurrentUserDep,
    tsio: TsioDep,
    frame_cleanup: FrameRoomCleanupDep,
):
    """Delete worker, fail their tasks, remove job links, and clean up orphan jobs.

    Worker must belong to authenticated user or user must be superuser.
    """
    result = await session.exec(select(Worker).where(Worker.id == worker_id))
    worker = result.one_or_none()
    if not worker:
        raise WorkerNotFound.exception(detail=f"Worker '{worker_id}' not found")

    if worker.user_id != user.id and not user.is_superuser:
        raise Forbidden.exception(detail="Worker belongs to different user")

    emissions, frame_rooms = await cleanup_worker(session, worker)
    await session.commit()
    await emit(tsio, emissions)
    if frame_rooms:
        await frame_cleanup(frame_rooms)


# ---------------------------------------------------------------------------
# Provider endpoints
# ---------------------------------------------------------------------------


def _room_provider_filter(room_id: str):
    """Build a SQLAlchemy filter for providers visible from a given room."""
    if room_id == "@global":
        return ProviderRecord.room_id == "@global"
    if room_id == "@internal":
        return ProviderRecord.room_id == "@internal"
    return ProviderRecord.room_id.in_(["@global", "@internal", room_id])


def _require_internal_filesystem_access(
    provider: ProviderRecord, user, settings
) -> None:
    """Gate @internal:filesystem:* access on superuser status.

    Parameters
    ----------
    provider : ProviderRecord
        The provider record being accessed.
    user : User
        The authenticated user making the request.
    settings : JobLibSettings
        Current joblib settings (reads ``filebrowser_require_superuser``).

    Raises
    ------
    ProblemError
        403 Forbidden when the caller is non-superuser and the
        ``filebrowser_require_superuser`` setting is True.
    """
    if (
        provider.room_id == "@internal"
        and provider.category == "filesystem"
        and settings.filebrowser_require_superuser
        and not user.is_superuser
    ):
        raise Forbidden.exception(
            detail="@internal filesystem access requires superuser"
        )


async def _resolve_provider(
    session: AsyncSession, provider_name: str, room_id: str
) -> ProviderRecord:
    """Resolve a provider full_name into a ProviderRecord, checking room visibility."""
    parts = provider_name.split(":", 2)
    if len(parts) != 3:
        raise ProviderNotFound.exception(
            detail=f"Invalid provider name format: {provider_name}"
        )
    provider_room_id, category, name = parts

    # Visibility check — mirror _room_provider_filter semantics.
    allowed = {"@global"} if room_id == "@global" else {"@global", "@internal", room_id}
    if provider_room_id not in allowed:
        raise ProviderNotFound.exception(
            detail=f"Provider '{provider_name}' not accessible from room '{room_id}'"
        )

    result = await session.exec(
        select(ProviderRecord).where(
            ProviderRecord.room_id == provider_room_id,
            ProviderRecord.category == category,
            ProviderRecord.name == name,
        )
    )
    provider = result.one_or_none()
    if not provider:
        raise ProviderNotFound.exception(detail=f"Provider '{provider_name}' not found")
    return provider


@router.put(
    "/rooms/{room_id}/providers",
    response_model=ProviderResponse,
    status_code=status.HTTP_201_CREATED,
)
async def register_provider(
    room_id: WritableRoomDep,
    request: ProviderRegisterRequest,
    response: Response,
    session: SessionDep,
    user: CurrentUserDep,
    settings: SettingsDep,
    tsio: TsioDep,
):
    """Register or update a provider. Idempotent on (room_id, category, name)."""
    if room_id == "@internal":
        raise Forbidden.exception(
            detail="@internal providers cannot be registered via HTTP"
        )
    if room_id == "@global" and not user.is_superuser:
        raise Forbidden.exception(
            detail="Admin required for @global provider registration"
        )

    if (
        settings.allowed_provider_categories is not None
        and request.category not in settings.allowed_provider_categories
    ):
        raise InvalidCategory.exception(
            detail=(
                f"Provider category '{request.category}' not in allowed list:"
                f" {settings.allowed_provider_categories}"
            )
        )

    # Handle worker: use provided worker_id or auto-create
    if request.worker_id:
        result = await session.exec(
            select(Worker).where(
                Worker.id == request.worker_id,
                Worker.user_id == user.id,
            )
        )
        worker = result.one_or_none()
        if not worker:
            raise WorkerNotFound.exception(
                detail=f"Worker '{request.worker_id}' not found or not owned by user"
            )
    else:
        worker = Worker(user_id=user.id)
        session.add(worker)
        await session.flush()

    # Check if provider already exists (upsert)
    result = await session.exec(
        select(ProviderRecord).where(
            ProviderRecord.category == request.category,
            ProviderRecord.name == request.name,
            ProviderRecord.room_id == room_id,
        )
    )
    existing = result.one_or_none()

    if existing:
        if existing.user_id != user.id and not user.is_superuser:
            raise Forbidden.exception(
                detail="Provider already registered by another user"
            )
        existing.schema_ = request.schema_
        existing.content_type = request.content_type
        existing.worker_id = worker.id
        existing.user_id = user.id
        provider = existing
        response.status_code = status.HTTP_200_OK
    else:
        provider = ProviderRecord(
            category=request.category,
            name=request.name,
            room_id=room_id,
            schema_=request.schema_,
            content_type=request.content_type,
            user_id=user.id,
            worker_id=worker.id,
        )
        session.add(provider)

    await session.commit()
    await session.refresh(provider)
    await emit(tsio, {Emission(ProvidersInvalidate(), f"room:{provider.room_id}")})

    return ProviderResponse.from_record(provider)


@router.get(
    "/rooms/{room_id}/providers",
    response_model=PaginatedResponse[ProviderResponse],
)
async def list_providers(
    room_id: str,
    session: SessionDep,
    _current_user: CurrentUserDep,
    settings: SettingsDep,
    limit: Annotated[int, Query(ge=0, le=500)] = 50,
    offset: Annotated[int, Query(ge=0)] = 0,
):
    """List providers visible from a room (room-scoped + @global)."""
    validate_room_id(room_id)

    filter_ = _room_provider_filter(room_id)
    # Hide @internal:filesystem:* from non-superusers when the gate is on.
    if settings.filebrowser_require_superuser and not _current_user.is_superuser:
        gate_filter = ~(
            (ProviderRecord.room_id == "@internal")
            & (ProviderRecord.category == "filesystem")
        )
        filter_ = and_(filter_, gate_filter)

    base_query = select(ProviderRecord).where(filter_)

    total_result = await session.exec(
        select(func.count()).select_from(base_query.subquery())
    )
    total = total_result.one()

    result = await session.exec(
        base_query.order_by(ProviderRecord.created_at.desc())
        .offset(offset)
        .limit(limit)
    )
    providers = result.all()

    items = [ProviderResponse.from_record(p) for p in providers]
    return PaginatedResponse(items=items, total=total, limit=limit, offset=offset)


@router.get(
    "/rooms/{room_id}/providers/{provider_name:path}/info",
    response_model=ProviderResponse,
)
async def get_provider_info(
    room_id: str,
    provider_name: str,
    session: SessionDep,
    _current_user: CurrentUserDep,
    settings: SettingsDep,
):
    """Get provider details and JSON Schema."""
    validate_room_id(room_id)
    provider = await _resolve_provider(session, provider_name, room_id)
    _require_internal_filesystem_access(provider, _current_user, settings)

    return ProviderResponse.from_record(provider)


@router.get("/rooms/{room_id}/providers/{provider_name:path}")
async def read_provider(
    room_id: str,
    provider_name: str,
    request: Request,
    session_maker: SessionMakerDep,
    _current_user: CurrentUserFactoryDep,
    result_backend: ResultBackendDep,
    settings: SettingsDep,
    tsio: TsioDep,
    internal_provider_registry: InternalProviderRegistryDep,
    prefer: Annotated[str | None, Header()] = None,
):
    """Read data from a provider. Long-polls until result is available."""
    validate_room_id(room_id)

    # Short-lived session — closed before long-poll
    async with session_maker() as session:
        provider = await _resolve_provider(session, provider_name, room_id)

    _require_internal_filesystem_access(provider, _current_user, settings)

    params = dict(request.query_params)
    rhash = request_hash(params)
    cache_key = f"provider-result:{provider.full_name}:{rhash}"
    status_key = f"{cache_key}:status"
    inflight_key = f"provider-inflight:{provider.full_name}:{rhash}"

    # Fast path: cache hit
    cached = await result_backend.get(cache_key)
    if cached is not None:
        cached_status = await result_backend.get(status_key)
        if cached_status == b"error":
            try:
                problem = json.loads(cached)
                status_code = int(problem.get("status", 500))
            except (ValueError, TypeError, KeyError):
                status_code = 500
            return Response(
                content=cached,
                media_type=ProblemDetail.MEDIA_TYPE,
                status_code=status_code,
            )
        return Response(content=cached, media_type=provider.content_type)

    # Dispatch if not already inflight
    acquired = await result_backend.acquire_inflight(
        inflight_key, settings.provider_inflight_ttl_seconds
    )
    if acquired:
        try:
            if provider.room_id == "@internal":
                if (
                    internal_provider_registry is None
                    or provider.full_name not in internal_provider_registry.tasks
                ):
                    raise InternalJobNotConfigured.exception(  # noqa: TRY301
                        detail=(
                            f"Internal provider '{provider.full_name}' is registered"
                            " in the DB but no executor task is available"
                        )
                    )
                params_json = json.dumps(params, sort_keys=True, separators=(",", ":"))
                worker_token = await mint_internal_worker_token(request.app)
                await internal_provider_registry.tasks[provider.full_name].kiq(
                    request_id=rhash,
                    provider_id=str(provider.id),
                    params_json=params_json,
                    token=worker_token,
                )
            else:
                provider_room = f"providers:{provider.full_name}"
                await emit(
                    tsio,
                    {
                        Emission(
                            ProviderRequest.from_dict_params(
                                request_id=rhash,
                                provider_name=provider.full_name,
                                params=params,
                            ),
                            provider_room,
                        )
                    },
                )
        except BaseException:
            await result_backend.release_inflight(inflight_key)
            raise

    # Long-poll: wait for result via pub/sub
    requested_wait = parse_prefer_wait(prefer)
    if requested_wait is not None:
        timeout = min(requested_wait, settings.provider_long_poll_max_seconds)
    else:
        timeout = settings.provider_long_poll_default_seconds

    result = await result_backend.wait_for_key(cache_key, timeout)
    if result is not None:
        result_status = await result_backend.get(status_key)
        headers = {}
        if requested_wait is not None:
            headers["Preference-Applied"] = f"wait={int(timeout)}"
        if result_status == b"error":
            try:
                problem = json.loads(result)
                status_code = int(problem.get("status", 500))
            except (ValueError, TypeError, KeyError):
                status_code = 500
            return Response(
                content=result,
                media_type=ProblemDetail.MEDIA_TYPE,
                status_code=status_code,
                headers=headers,
            )
        return Response(
            content=result, media_type=provider.content_type, headers=headers
        )

    # Timeout — RFC 9457 error with retry guidance
    raise ProviderTimeout.exception(
        detail=f"Provider '{provider_name}' did not respond within {timeout}s",
        headers={"Retry-After": "2"},
    )


@router.delete("/providers/{provider_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_provider(
    provider_id: UUID,
    session: SessionDep,
    user: CurrentUserDep,
    tsio: TsioDep,
):
    """Unregister a provider. Must be owned by authenticated user or superuser."""
    result = await session.exec(
        select(ProviderRecord).where(ProviderRecord.id == provider_id)
    )
    provider = result.one_or_none()
    if not provider:
        raise ProviderNotFound.exception(detail=f"Provider '{provider_id}' not found")

    if provider.room_id == "@internal":
        raise Forbidden.exception(detail="@internal providers cannot be deleted")

    if provider.user_id != user.id and not user.is_superuser:
        raise Forbidden.exception(detail="Provider belongs to different user")

    room_id = provider.room_id
    await session.delete(provider)
    await session.commit()
    await emit(tsio, {Emission(ProvidersInvalidate(), f"room:{room_id}")})


@router.post(
    "/providers/{provider_id}/results",
    status_code=status.HTTP_204_NO_CONTENT,
)
async def upload_provider_result(
    provider_id: UUID,
    request: Request,
    session_maker: SessionMakerDep,
    user: CurrentUserFactoryDep,
    result_backend: ResultBackendDep,
    settings: SettingsDep,
    tsio: TsioDep,
    x_request_hash: Annotated[str, Header()],
    x_result_status: Annotated[str | None, Header()] = None,
):
    """Provider worker uploads a read result."""
    async with session_maker() as session:
        result = await session.exec(
            select(ProviderRecord).where(ProviderRecord.id == provider_id)
        )
        provider = result.one_or_none()
        if not provider:
            raise ProviderNotFound.exception(
                detail=f"Provider '{provider_id}' not found"
            )

        if provider.user_id != user.id and not user.is_superuser:
            raise Forbidden.exception(
                detail="Not authorized to upload results for this provider"
            )

    cache_key = f"provider-result:{provider.full_name}:{x_request_hash}"
    status_key = f"{cache_key}:status"
    inflight_key = f"provider-inflight:{provider.full_name}:{x_request_hash}"

    # Store raw body as-is (JSON or binary depending on provider content_type)
    data = await request.body()
    # Write status FIRST so any reader that wakes up on cache_key sees
    # the status already populated.
    if x_result_status == "error":
        await result_backend.store(
            status_key,
            b"error",
            settings.provider_result_ttl_seconds,
        )
    await result_backend.store(
        cache_key,
        data,
        settings.provider_result_ttl_seconds,
    )

    # Release inflight lock
    await result_backend.release_inflight(inflight_key)

    # Wake long-polling waiters (Redis pub/sub)
    await result_backend.notify_key(cache_key)

    # Notify frontend (Socket.IO — UI refresh)
    await emit(
        tsio,
        {
            Emission(
                ProviderResultReady(
                    provider_name=provider.full_name,
                    request_hash=x_request_hash,
                ),
                f"room:{provider.room_id}",
            )
        },
    )
