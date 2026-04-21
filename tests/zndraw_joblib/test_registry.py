"""Tests for internal job registry."""

import asyncio
import uuid
from typing import Any, ClassVar
from unittest.mock import AsyncMock, MagicMock

import pytest
from sqlmodel import select

from zndraw_auth import User
from zndraw_joblib.client import Category, Extension
from zndraw_joblib.models import Job, ProviderRecord
from zndraw_joblib.provider import Provider
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
    token: str,
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
        result = await session.exec(select(Job).where(Job.room_id == "@internal"))
        jobs = result.all()

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
        result = await session.exec(
            select(Job).where(Job.room_id == "@internal", Job.name == "Rotate")
        )
        job = result.one()
        assert not job.deleted


# --- Provider registry ---------------------------------------------------


class _DummyProvider(Provider):
    category: ClassVar[str] = "filesystem"
    path: str = "/"

    def read(self, handler):
        return handler.ls(self.path, detail=True)


async def test_register_internal_providers_registers_taskiq_task():
    """Each provider class becomes a task at @internal:<cat>:<ClassName>."""
    from taskiq import InMemoryBroker

    from zndraw_joblib.registry import register_internal_providers

    broker = InMemoryBroker()

    calls = []

    async def executor(cls, params_json, provider_id, request_id, token):
        calls.append((cls, params_json, provider_id, request_id, token))

    reg = register_internal_providers(broker, [_DummyProvider], executor)

    assert "@internal:filesystem:_DummyProvider" in reg.tasks
    assert reg.providers["@internal:filesystem:_DummyProvider"] is _DummyProvider


async def test_register_internal_providers_task_invokes_executor():
    """Calling the registered task fan-outs to the executor with forwarded args."""
    from taskiq import InMemoryBroker

    from zndraw_joblib.registry import register_internal_providers

    broker = InMemoryBroker()
    calls = []

    async def executor(cls, params_json, provider_id, request_id, token):
        calls.append((cls.__name__, params_json, provider_id, request_id, token))

    reg = register_internal_providers(broker, [_DummyProvider], executor)
    task = reg.tasks["@internal:filesystem:_DummyProvider"]

    # Directly invoke the underlying coroutine (InMemoryBroker doesn't run it)
    await task.original_func(
        request_id="abc",
        provider_id="11111111-1111-1111-1111-111111111111",
        params_json='{"path": "/data"}',
        token="tok",
    )

    assert calls == [
        (
            "_DummyProvider",
            '{"path": "/data"}',
            "11111111-1111-1111-1111-111111111111",
            "abc",
            "tok",
        )
    ]


async def test_ensure_internal_providers_creates_rows(async_session_factory):
    """Creates a ProviderRecord per provider at room_id='@internal'."""
    from zndraw_joblib.registry import ensure_internal_providers

    # Seed internal user
    user_id = uuid.uuid4()
    async with async_session_factory() as session:
        user = User(
            id=user_id,
            email="internal@test",
            hashed_password="x",
            is_active=True,
            is_superuser=True,
            is_verified=True,
        )
        session.add(user)
        await session.commit()

    await ensure_internal_providers(
        [_DummyProvider], async_session_factory, user_id=user_id
    )

    async with async_session_factory() as session:
        result = await session.exec(
            select(ProviderRecord).where(ProviderRecord.room_id == "@internal")
        )
        rows = result.all()
    assert len(rows) == 1
    assert rows[0].category == "filesystem"
    assert rows[0].name == "_DummyProvider"
    assert rows[0].user_id == user_id
    assert rows[0].worker_id is None


async def test_ensure_internal_providers_idempotent(async_session_factory):
    """Running twice leaves exactly one row (upsert on room+category+name)."""
    from zndraw_joblib.registry import ensure_internal_providers

    user_id = uuid.uuid4()
    async with async_session_factory() as session:
        session.add(
            User(
                id=user_id,
                email="internal@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
        )
        await session.commit()

    await ensure_internal_providers(
        [_DummyProvider], async_session_factory, user_id=user_id
    )
    await ensure_internal_providers(
        [_DummyProvider], async_session_factory, user_id=user_id
    )

    async with async_session_factory() as session:
        result = await session.exec(
            select(ProviderRecord).where(ProviderRecord.room_id == "@internal")
        )
        rows = result.all()
    assert len(rows) == 1


async def test_ensure_internal_providers_concurrent_startup_safe(tmp_path):
    """Two concurrent ensure_internal_providers calls produce exactly one row.

    Models the multi-replica startup race on ProviderRecord's unique constraint
    (room_id, category, name). The IntegrityError on the loser must be caught
    and the row updated, not crash startup.

    Uses NullPool (a new connection per checkout) so that the two concurrent
    sessions get independent connections, as they would in a real multi-replica
    deploy.  StaticPool shares a single connection across all sessions, which
    makes the rollback on the loser undo the winner's commit — not a realistic
    scenario and not the behaviour we're guarding against.
    """
    from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
    from sqlalchemy.pool import NullPool
    from sqlmodel import SQLModel, select
    from sqlmodel.ext.asyncio.session import AsyncSession

    from zndraw.providers import BUNDLED_PROVIDERS
    from zndraw_joblib.models import ProviderRecord
    from zndraw_joblib.registry import ensure_internal_providers

    if not BUNDLED_PROVIDERS:
        pytest.skip("no providers bundled")
    prov_cls = next(iter(BUNDLED_PROVIDERS))

    db_path = tmp_path / "test_concurrent.db"
    engine = create_async_engine(
        f"sqlite+aiosqlite:///{db_path}",
        connect_args={"check_same_thread": False},
        poolclass=NullPool,
    )
    async with engine.begin() as conn:
        await conn.run_sync(SQLModel.metadata.create_all)
    maker = async_sessionmaker(engine, class_=AsyncSession, expire_on_commit=False)
    user_id = uuid.uuid4()
    try:
        await asyncio.gather(
            ensure_internal_providers([prov_cls], maker, user_id=user_id),
            ensure_internal_providers([prov_cls], maker, user_id=user_id),
        )
        async with maker() as s:
            rows = (
                await s.exec(
                    select(ProviderRecord).where(ProviderRecord.room_id == "@internal")
                )
            ).all()
        assert len(rows) == 1, f"expected exactly 1 row, got {len(rows)}: {rows}"
    finally:
        await engine.dispose()
