"""Internal job registry for taskiq-based server-side execution."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Protocol

if TYPE_CHECKING:
    from collections.abc import Callable
    from contextlib import AbstractAsyncContextManager
    from uuid import UUID

    from fastapi import FastAPI
    from sqlmodel.ext.asyncio.session import AsyncSession
    from taskiq import AsyncBroker

    from zndraw_joblib.client import Extension
    from zndraw_joblib.provider import Provider

logger = logging.getLogger(__name__)


class InternalExecutor(Protocol):
    """Protocol for the host-provided executor callback.

    The server base URL is captured at creation time.
    A per-task JWT token is passed at call time.
    """

    async def __call__(
        self,
        extension_cls: type[Extension],
        payload: dict[str, Any],
        room_id: str,
        task_id: str,
        token: str,
    ) -> None: ...


@dataclass
class InternalRegistry:
    """Holds taskiq task handles and extension class mappings."""

    tasks: dict[str, Any] = field(default_factory=dict)
    extensions: dict[str, type[Extension]] = field(default_factory=dict)
    executor: InternalExecutor | None = None


def register_internal_tasks(
    broker: AsyncBroker,
    extensions: list[type[Extension]],
    executor: InternalExecutor,
) -> InternalRegistry:
    """Register extension classes as taskiq tasks on the broker.

    For external taskiq worker processes that have no FastAPI app or DB.
    Also used internally by register_internal_jobs.
    """
    registry = InternalRegistry(executor=executor)

    for ext_cls in extensions:
        category = ext_cls.category.value
        name = ext_cls.__name__
        full_name = f"@internal:{category}:{name}"

        def _make_task_fn(
            cls: type[Extension] = ext_cls, ex: InternalExecutor = executor
        ):
            async def _execute(
                task_id: str, room_id: str, payload: dict[str, Any], token: str
            ) -> None:
                await ex(cls, payload, room_id, task_id, token)

            return _execute

        task_handle = broker.register_task(
            _make_task_fn(),
            task_name=full_name,
        )

        registry.tasks[full_name] = task_handle
        registry.extensions[full_name] = ext_cls
        logger.debug("Registered internal task: %s", full_name)

    logger.info("Registered %d internal task(s)", len(extensions))
    return registry


class InternalProviderExecutor(Protocol):
    """Protocol for the host-provided provider executor.

    The server base URL and any handler configuration are captured at
    creation time. Per-request data (params, ids, token) is passed at
    call time.
    """

    async def __call__(
        self,
        provider_cls: type[Provider],
        params_json: str,
        provider_id: str,
        request_id: str,
        token: str,
    ) -> None: ...


@dataclass
class InternalProviderRegistry:
    """Holds taskiq task handles and provider class mappings."""

    tasks: dict[str, Any] = field(default_factory=dict)
    providers: dict[str, type[Provider]] = field(default_factory=dict)
    executor: InternalProviderExecutor | None = None


def register_internal_providers(
    broker: AsyncBroker,
    providers: list[type[Provider]],
    executor: InternalProviderExecutor,
) -> InternalProviderRegistry:
    """Register Provider classes as taskiq tasks on the broker.

    Each class becomes a task named ``@internal:<category>:<ClassName>``.
    The task handler forwards ``(cls, params_json, provider_id, request_id,
    token)`` to *executor*.
    """
    registry = InternalProviderRegistry(executor=executor)

    for prov_cls in providers:
        category = prov_cls.category
        name = prov_cls.__name__
        full_name = f"@internal:{category}:{name}"

        def _make_task_fn(
            cls: type[Provider] = prov_cls,
            ex: InternalProviderExecutor = executor,
        ):
            async def _execute(
                request_id: str,
                provider_id: str,
                params_json: str,
                token: str,
            ) -> None:
                await ex(cls, params_json, provider_id, request_id, token)

            return _execute

        task_handle = broker.register_task(
            _make_task_fn(),
            task_name=full_name,
        )

        registry.tasks[full_name] = task_handle
        registry.providers[full_name] = prov_cls
        logger.debug("Registered internal provider task: %s", full_name)

    logger.info("Registered %d internal provider task(s)", len(providers))
    return registry


async def ensure_internal_providers(
    providers: list[type[Provider]],
    session_factory: Callable[[], AbstractAsyncContextManager[AsyncSession]],
    *,
    user_id: UUID,
) -> None:
    """Create or update @internal ProviderRecord rows.

    @internal providers are server-owned (dispatched by the in-process
    taskiq worker) and have no Worker row — ``worker_id`` is always
    ``None``. Idempotent and concurrent-startup-safe: an IntegrityError
    on the ``unique_provider`` constraint is caught and the existing
    row is updated instead.
    """
    from sqlalchemy.exc import IntegrityError
    from sqlmodel import select

    from zndraw_joblib.models import ProviderRecord

    for prov_cls in providers:
        category = prov_cls.category
        name = prov_cls.__name__
        schema = prov_cls.model_json_schema()
        content_type = prov_cls.content_type

        async with session_factory() as session:
            result = await session.exec(
                select(ProviderRecord).where(
                    ProviderRecord.room_id == "@internal",
                    ProviderRecord.category == category,
                    ProviderRecord.name == name,
                )
            )
            existing = result.one_or_none()

            if existing is not None:
                existing.schema_ = schema
                existing.content_type = content_type
                existing.user_id = user_id
                existing.worker_id = None
                await session.commit()
                continue

            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category=category,
                    name=name,
                    schema_=schema,
                    content_type=content_type,
                    user_id=user_id,
                    worker_id=None,
                )
            )
            try:
                await session.commit()
            except IntegrityError:
                await session.rollback()
                # Lost the race; requery + update.
                result = await session.exec(
                    select(ProviderRecord).where(
                        ProviderRecord.room_id == "@internal",
                        ProviderRecord.category == category,
                        ProviderRecord.name == name,
                    )
                )
                existing = result.one()
                existing.schema_ = schema
                existing.content_type = content_type
                existing.user_id = user_id
                existing.worker_id = None
                await session.commit()

    logger.info("Ensured %d @internal provider row(s) in DB", len(providers))


async def ensure_internal_jobs(
    extensions: list[type[Extension]],
    session_factory: Callable[[], AbstractAsyncContextManager[AsyncSession]],
) -> None:
    """Create or update @internal Job rows in the database.

    Idempotent — safe to call on every startup. For production with multiple
    replicas, call once from ``init_database()`` / ``zndraw-db`` to avoid
    race conditions on the ``unique_job`` constraint.
    """
    from sqlmodel import select

    from zndraw_joblib.models import Job

    async with session_factory() as session:
        for ext_cls in extensions:
            category = ext_cls.category.value
            name = ext_cls.__name__
            schema = ext_cls.model_json_schema()

            result = await session.exec(
                select(Job).where(
                    Job.room_id == "@internal",
                    Job.category == category,
                    Job.name == name,
                )
            )
            existing = result.one_or_none()

            if existing and existing.deleted:
                existing.deleted = False
                existing.schema_ = schema
            elif existing:
                existing.schema_ = schema
            else:
                session.add(
                    Job(
                        room_id="@internal",
                        category=category,
                        name=name,
                        schema_=schema,
                    )
                )

        await session.commit()

    logger.info("Ensured %d @internal job row(s) in DB", len(extensions))


async def register_internal_jobs(
    app: FastAPI,
    broker: AsyncBroker,
    extensions: list[type[Extension]],
    executor: InternalExecutor,
    session_factory: Callable[[], AbstractAsyncContextManager[AsyncSession]],
) -> None:
    """Register internal extensions for server-side execution.

    Does three things:
    1. Registers each Extension as a taskiq task on the broker
    2. Creates/reactivates @internal:category:name Job rows in the DB
    3. Stores the InternalRegistry on app.state.internal_registry
    """
    registry = register_internal_tasks(broker, extensions, executor)
    await ensure_internal_jobs(extensions, session_factory)
    app.state.internal_registry = registry
