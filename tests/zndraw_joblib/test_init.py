# tests/test_init.py
"""Tests for public API exports."""


def test_public_api_exports():
    """All public API should be importable from the main package."""
    from zndraw_joblib import (
        Category,
        # Models
        Job,
        JobManager,
        # Router
        router,
    )

    assert router is not None
    assert Job is not None
    assert JobManager is not None
    assert Category.MODIFIER.value == "modifiers"


def test_all_exports_in_dunder_all():
    """__all__ should contain all public exports."""
    import zndraw_joblib

    expected = [
        "router",
        "Job",
        "Worker",
        "Task",
        "WorkerJobLink",
        "TaskStatus",
        "get_joblib_settings",
        "JobLibSettingsDep",
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
        "JobLibSettings",
        "JobManager",
        "ClaimedTask",
        "Extension",
        "Category",
        "InternalJobNotConfigured",
        "register_internal_jobs",
        "register_internal_tasks",
        "InternalExecutor",
        "InternalRegistry",
        "run_sweeper",
        "cleanup_stale_workers",
        "cleanup_stuck_internal_tasks",
        "get_tsio",
        "JobsInvalidate",
        "TaskAvailable",
        "TaskStatusEvent",
        "JoinJobRoom",
        "LeaveJobRoom",
        "Emission",
    ]

    for name in expected:
        assert name in zndraw_joblib.__all__, f"{name} not in __all__"
        assert hasattr(zndraw_joblib, name), f"{name} not accessible"


def test_internal_registry_exports():
    from zndraw_joblib import (
        InternalExecutor,
        InternalJobNotConfigured,
        InternalRegistry,
        cleanup_stuck_internal_tasks,
        register_internal_jobs,
        register_internal_tasks,
    )

    assert callable(register_internal_jobs)
    assert callable(register_internal_tasks)
    assert callable(cleanup_stuck_internal_tasks)
    assert InternalExecutor is not None
    assert InternalRegistry is not None
    assert InternalJobNotConfigured is not None
