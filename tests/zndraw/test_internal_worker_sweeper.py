"""Regression: @internal providers are decoupled from the Worker table.

@internal providers are server-owned (dispatched by the in-process
taskiq worker, not a remote client). They have ``worker_id=None`` by
construction, mirroring how @internal jobs have no ``WorkerJobLink``
entries.

This test pins down two invariants the sweeper must respect:

1. A stale REMOTE worker (with a valid UUID and an owned provider) IS
   swept, and the remote provider is cascade-deleted. No regression to
   the heartbeat-based cleanup that governs real external clients.
2. A server-owned @internal provider (``worker_id=None``) IS NOT
   swept under any circumstance — it has no Worker to be stale.
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
async def test_sweeper_deletes_remote_worker_but_preserves_internal_provider():
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
            remote_user = User(
                id=uuid.uuid4(),
                email="remote@test",
                hashed_password="x",
                is_active=True,
                is_superuser=False,
                is_verified=True,
            )
            internal_user = User(
                id=uuid.uuid4(),
                email="internal@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add_all([remote_user, internal_user])
            await session.commit()

            # Stale remote worker with a real worker_id and an owned provider.
            stale_remote_worker = Worker(
                id=uuid.uuid4(),
                user_id=remote_user.id,
                last_heartbeat=datetime.now(UTC) - timedelta(minutes=10),
            )
            session.add(stale_remote_worker)
            await session.commit()
            await session.refresh(stale_remote_worker)

            remote_provider = ProviderRecord(
                room_id="room-alpha",
                category="filesystem",
                name="FilesystemRead",
                schema_={},
                content_type="application/json",
                user_id=remote_user.id,
                worker_id=stale_remote_worker.id,
            )
            # @internal provider — server-owned, worker_id is None.
            internal_provider = ProviderRecord(
                room_id="@internal",
                category="filesystem",
                name="FilesystemRead",
                schema_={},
                content_type="application/json",
                user_id=internal_user.id,
                worker_id=None,
            )
            session.add_all([remote_provider, internal_provider])
            await session.commit()

        async with maker() as session:
            count, _emissions, _frame_rooms = await cleanup_stale_workers(
                session, stale_after=timedelta(seconds=0)
            )
            assert count == 1, (
                f"Expected the stale remote worker to be swept; count={count}"
            )

        async with maker() as session:
            rows = (await session.exec(select(ProviderRecord))).all()
            by_room = {r.room_id: r for r in rows}
            assert "room-alpha" not in by_room, (
                "Remote provider should have been cascade-deleted with its worker"
            )
            assert "@internal" in by_room, (
                "@internal provider (worker_id=None) must survive the sweep"
            )
            assert by_room["@internal"].worker_id is None
    finally:
        await engine.dispose()
