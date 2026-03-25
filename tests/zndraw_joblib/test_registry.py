"""Tests for internal job registry."""

from typing import Any, ClassVar
from unittest.mock import AsyncMock, MagicMock

from sqlalchemy import select

from zndraw_joblib.client import Category, Extension
from zndraw_joblib.models import Job
from zndraw_joblib.registry import (
    InternalRegistry,
    register_internal_jobs,
    register_internal_tasks,
)


class ConcreteExtension(Extension):
    """Non-abstract extension base for tests that don't need run()."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self) -> None:
        pass


class Rotate(ConcreteExtension):
    category: ClassVar[Category] = Category.MODIFIER
    angle: float = 0.0


class Scale(ConcreteExtension):
    category: ClassVar[Category] = Category.MODIFIER
    factor: float = 1.0


def make_mock_broker():
    """Create a mock broker with register_task that returns a mock task handle."""
    broker = MagicMock()
    task_handles = {}

    def mock_register_task(fn, task_name, **kwargs):
        handle = MagicMock()
        handle.kiq = AsyncMock()
        task_handles[task_name] = handle
        return handle

    broker.register_task = mock_register_task
    broker._task_handles = task_handles
    return broker


async def mock_executor(
    extension_cls: type[Extension],
    payload: dict[str, Any],
    room_id: str,
    task_id: str,
) -> None:
    pass


def test_register_internal_tasks_returns_registry():
    broker = make_mock_broker()
    registry = register_internal_tasks(broker, [Rotate, Scale], executor=mock_executor)

    assert isinstance(registry, InternalRegistry)
    assert "@internal:modifiers:Rotate" in registry.tasks
    assert "@internal:modifiers:Scale" in registry.tasks
    assert len(registry.tasks) == 2


def test_register_internal_tasks_stores_extension_classes():
    broker = make_mock_broker()
    registry = register_internal_tasks(broker, [Rotate], executor=mock_executor)

    assert registry.extensions["@internal:modifiers:Rotate"] is Rotate


def test_register_internal_tasks_registers_on_broker():
    broker = make_mock_broker()
    register_internal_tasks(broker, [Rotate], executor=mock_executor)

    assert "@internal:modifiers:Rotate" in broker._task_handles


def test_register_internal_tasks_empty_list():
    broker = make_mock_broker()
    registry = register_internal_tasks(broker, [], executor=mock_executor)

    assert len(registry.tasks) == 0
    assert len(registry.extensions) == 0


def test_internal_registry_executor_stored():
    broker = make_mock_broker()
    registry = register_internal_tasks(broker, [Rotate], executor=mock_executor)

    assert registry.executor is mock_executor


async def test_register_internal_jobs_creates_db_rows(async_session_factory):
    """register_internal_jobs creates Job rows in the DB."""
    broker = make_mock_broker()
    app = MagicMock()
    app.state = MagicMock()

    await register_internal_jobs(
        app,
        broker,
        [Rotate, Scale],
        executor=mock_executor,
        session_factory=async_session_factory,
    )

    async with async_session_factory() as session:
        result = await session.execute(select(Job).where(Job.room_id == "@internal"))
        jobs = result.scalars().all()

    assert len(jobs) == 2
    names = {j.full_name for j in jobs}
    assert "@internal:modifiers:Rotate" in names
    assert "@internal:modifiers:Scale" in names
    assert all(not j.deleted for j in jobs)


async def test_register_internal_jobs_sets_app_state(async_session_factory):
    """register_internal_jobs stores registry on app.state."""
    broker = make_mock_broker()
    app = MagicMock()
    app.state = MagicMock()

    await register_internal_jobs(
        app,
        broker,
        [Rotate],
        executor=mock_executor,
        session_factory=async_session_factory,
    )

    assert hasattr(app.state, "internal_registry")
    registry = app.state.internal_registry
    assert isinstance(registry, InternalRegistry)
    assert "@internal:modifiers:Rotate" in registry.tasks


async def test_register_internal_jobs_reactivates_deleted(async_session_factory):
    """register_internal_jobs reactivates soft-deleted jobs."""
    # First: create a deleted job manually
    async with async_session_factory() as session:
        job = Job(
            room_id="@internal",
            category="modifiers",
            name="Rotate",
            schema_={},
            deleted=True,
        )
        session.add(job)
        await session.commit()

    broker = make_mock_broker()
    app = MagicMock()
    app.state = MagicMock()

    await register_internal_jobs(
        app,
        broker,
        [Rotate],
        executor=mock_executor,
        session_factory=async_session_factory,
    )

    async with async_session_factory() as session:
        result = await session.execute(
            select(Job).where(Job.room_id == "@internal", Job.name == "Rotate")
        )
        job = result.scalar_one()
        assert not job.deleted
