"""Regression: the internal @internal worker must survive the sweeper.

The internal worker is seeded at startup (``ensure_internal_worker_row``) and
never heartbeats. The sweeper's ``cleanup_stale_workers`` runs every 30s by
default and considers any worker with ``last_heartbeat < now - 60s`` stale;
when it finds one, ``cleanup_worker`` DELETEs every provider owned by that
worker — including the seeded ``@internal:filesystem:FilesystemRead`` row.
Without a guard, the default FS provider disappears from every fresh server
after ~60 seconds and GETs start returning 404.
"""

from __future__ import annotations

import uuid
from datetime import UTC, datetime, timedelta

import pytest
from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel, select
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw_auth import User
from zndraw_joblib.models import ProviderRecord, Worker
from zndraw_joblib.sweeper import cleanup_stale_workers


@pytest.mark.asyncio
async def test_sweeper_does_not_delete_internal_worker_or_its_providers():
    """An internal worker with stale heartbeat MUST survive the sweeper."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    async with engine.begin() as conn:
        await conn.run_sync(SQLModel.metadata.create_all)
    maker = async_sessionmaker(engine, class_=AsyncSession, expire_on_commit=False)
    try:
        async with maker() as session:
            user = User(
                id=uuid.uuid4(),
                email="internal@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.commit()

            # Stale internal worker (heartbeat 10 minutes ago).
            worker = Worker(
                id=uuid.uuid4(),
                user_id=user.id,
                last_heartbeat=datetime.now(UTC) - timedelta(minutes=10),
            )
            session.add(worker)
            await session.commit()
            await session.refresh(worker)

            provider = ProviderRecord(
                room_id="@internal",
                category="filesystem",
                name="FilesystemRead",
                schema_={},
                content_type="application/json",
                user_id=user.id,
                worker_id=worker.id,
            )
            session.add(provider)
            await session.commit()

        async with maker() as session:
            count, _, _ = await cleanup_stale_workers(
                session, stale_after=timedelta(seconds=60)
            )
            assert count == 0, (
                "Internal worker was swept; its @internal providers would be deleted"
            )

        # @internal provider row must still exist.
        async with maker() as session:
            rows = (
                await session.exec(
                    select(ProviderRecord).where(ProviderRecord.room_id == "@internal")
                )
            ).all()
            assert len(rows) == 1, (
                f"Expected @internal provider to survive; got {len(rows)} rows"
            )
    finally:
        await engine.dispose()
