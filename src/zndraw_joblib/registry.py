"""Internal job registry for taskiq-based server-side execution."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Protocol

if TYPE_CHECKING:
    from collections.abc import Callable
    from contextlib import AbstractAsyncContextManager

    from fastapi import FastAPI
    from sqlalchemy.ext.asyncio import AsyncSession

from taskiq import AsyncBroker

from zndraw_joblib.client import Extension

logger = logging.getLogger(__name__)


class InternalExecutor(Protocol):
    """Protocol for the host-provided executor callback.

    Deployment config (base URL, credentials) is captured at creation time.
    Only per-task data is passed at call time.
    """

    async def __call__(
        self,
        extension_cls: type[Extension],
        payload: dict[str, Any],
        room_id: str,
        task_id: str,
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
                task_id: str, room_id: str, payload: dict[str, Any]
            ) -> None:
                await ex(cls, payload, room_id, task_id)

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


async def ensure_internal_jobs(
    extensions: list[type[Extension]],
    session_factory: "Callable[[], AbstractAsyncContextManager[AsyncSession]]",
) -> None:
    """Create or update @internal Job rows in the database.

    Idempotent — safe to call on every startup. For production with multiple
    replicas, call once from ``init_database()`` / ``zndraw-db`` to avoid
    race conditions on the ``unique_job`` constraint.
    """
    from sqlalchemy import select

    from zndraw_joblib.models import Job

    async with session_factory() as session:
        for ext_cls in extensions:
            category = ext_cls.category.value
            name = ext_cls.__name__
            schema = ext_cls.model_json_schema()

            result = await session.execute(
                select(Job).where(
                    Job.room_id == "@internal",
                    Job.category == category,
                    Job.name == name,
                )
            )
            existing = result.scalar_one_or_none()

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
    app: "FastAPI",
    broker: AsyncBroker,
    extensions: list[type[Extension]],
    executor: InternalExecutor,
    session_factory: "Callable[[], AbstractAsyncContextManager[AsyncSession]]",
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
