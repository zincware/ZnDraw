# src/zndraw_joblib/__init__.py
"""ZnDraw Job Management Library."""

from zndraw._version import __version__ as __version__

from zndraw_joblib.client import Category, ClaimedTask, Extension, JobManager
from zndraw_joblib.dependencies import (
    FrameRoomCleanup,
    FrameRoomCleanupDep,
    JobLibSettingsDep,
    ResultBackend,
    ResultBackendDep,
    WritableRoomDep,
    get_frame_room_cleanup,
    get_joblib_settings,
    get_result_backend,
    get_tsio,
    request_hash,
    validate_room_id,
    verify_writable_room,
)
from zndraw_joblib.events import (
    Emission,
    FrozenEvent,
    JobsInvalidate,
    JoinJobRoom,
    JoinProviderRoom,
    LeaveJobRoom,
    LeaveProviderRoom,
    ProviderRequest,
    ProviderResultReady,
    ProvidersInvalidate,
    TaskAvailable,
    TaskStatusEvent,
    build_task_status_emission,
    emit,
)
from zndraw_joblib.exceptions import (
    Forbidden,
    InternalJobNotConfigured,
    InvalidCategory,
    InvalidRoomId,
    InvalidTaskTransition,
    JobNotFound,
    ProblemException,
    ProviderNotFound,
    SchemaConflict,
    TaskNotFound,
    WorkerNotFound,
    problem_exception_handler,
)
from zndraw_joblib.exceptions import (
    ProviderTimeoutError as ProviderTimeoutError,
)
from zndraw_joblib.models import (
    Job,
    ProviderRecord,
    Task,
    TaskStatus,
    Worker,
    WorkerJobLink,
)
from zndraw_joblib.provider import Provider
from zndraw_joblib.registry import (
    InternalExecutor,
    InternalRegistry,
    register_internal_jobs,
    register_internal_tasks,
)
from zndraw_joblib.router import router
from zndraw_joblib.schemas import PaginatedResponse
from zndraw_joblib.settings import JobLibSettings
from zndraw_joblib.sweeper import (
    cleanup_stale_workers,
    cleanup_stuck_internal_tasks,
    cleanup_worker,
    run_sweeper,
)

__all__ = [
    # Router
    "router",
    # Models
    "Job",
    "Worker",
    "Task",
    "WorkerJobLink",
    "TaskStatus",
    "ProviderRecord",
    # Provider
    "Provider",
    # Dependencies
    "get_joblib_settings",
    "JobLibSettingsDep",
    "get_tsio",
    "get_result_backend",
    "ResultBackend",
    "ResultBackendDep",
    "get_frame_room_cleanup",
    "FrameRoomCleanup",
    "FrameRoomCleanupDep",
    "request_hash",
    "verify_writable_room",
    "WritableRoomDep",
    "validate_room_id",
    # Exceptions
    "ProblemException",
    "problem_exception_handler",
    "JobNotFound",
    "SchemaConflict",
    "InvalidCategory",
    "WorkerNotFound",
    "TaskNotFound",
    "InvalidTaskTransition",
    "InvalidRoomId",
    "Forbidden",
    "InternalJobNotConfigured",
    "ProviderNotFound",
    "ProviderTimeoutError",
    # Schemas
    "PaginatedResponse",
    # Settings
    "JobLibSettings",
    # Client
    "JobManager",
    "ClaimedTask",
    "Extension",
    "Category",
    # Internal registry
    "register_internal_jobs",
    "register_internal_tasks",
    "InternalExecutor",
    "InternalRegistry",
    # Sweeper
    "run_sweeper",
    "cleanup_stale_workers",
    "cleanup_stuck_internal_tasks",
    "cleanup_worker",
    # Events
    "FrozenEvent",
    "JobsInvalidate",
    "ProvidersInvalidate",
    "ProviderRequest",
    "ProviderResultReady",
    "TaskAvailable",
    "TaskStatusEvent",
    "JoinJobRoom",
    "LeaveJobRoom",
    "JoinProviderRoom",
    "LeaveProviderRoom",
    "Emission",
    "build_task_status_emission",
    "emit",
]
