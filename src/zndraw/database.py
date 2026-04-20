"""Async database setup and application lifespan management.

Uses app.state pattern to store resources:
- redis: Redis async client
- frame_storage: FrameStorage (room-scoped AsyncBlobIO registry)
- fake_server: TcpFakeServer instance (when REDIS_URL not configured)
- tsio: zndraw-socketio typed wrapper
"""

import asyncio
import contextlib
import importlib
import logging
import socket
import threading
import uuid
from collections.abc import AsyncIterator

import redis.asyncio as redis_client
import socketio as socketio_lib
from fastapi import FastAPI
from fastapi_users.password import PasswordHelper
from sqlalchemy.ext.asyncio import AsyncEngine, async_sessionmaker
from sqlmodel import SQLModel, select
from sqlmodel.ext.asyncio.session import (
    AsyncSession,
    AsyncSession as SQLModelAsyncSession,
)
from taskiq.api import run_receiver_task
from taskiq_redis import ListQueueBroker

import zndraw.models  # noqa: F401 - registers Room, Message, etc.
import zndraw_joblib.models  # noqa: F401 - registers Job, Worker, Task
from zndraw.config import Settings
from zndraw.executor import InternalExtensionExecutor
from zndraw.extensions.analysis import analysis
from zndraw.extensions.modifiers import modifiers
from zndraw.extensions.selections import selections
from zndraw.providers import BUNDLED_PROVIDERS
from zndraw.redis import RedisKey
from zndraw.socket_events import FramesInvalidate
from zndraw.socketio import tsio
from zndraw_auth import User
from zndraw_auth.db import create_engine_for_url, ensure_default_admin
from zndraw_auth.settings import AuthSettings
from zndraw_joblib import (
    JoinJobRoom,
    JoinProviderRoom,
    LeaveJobRoom,
    LeaveProviderRoom,
    register_internal_jobs,
    run_sweeper,
)
from zndraw_joblib.registry import ensure_internal_jobs, register_internal_tasks
from zndraw_joblib.settings import JobLibSettings


def _get_free_port() -> int:
    """Find a free port on localhost using socket binding."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _collect_extensions() -> list[type]:
    """Collect all built-in extension classes."""
    return [
        *modifiers.values(),
        *selections.values(),
        *analysis.values(),
    ]


log = logging.getLogger(__name__)


def _is_sqlite(database_url: str) -> bool:
    """Check if database is SQLite (sync or async)."""
    return database_url.startswith(("sqlite://", "sqlite+aiosqlite://"))


def _apply_sqlite_locking(app: FastAPI) -> None:
    """Apply locking wrapper to session_maker for SQLite databases.

    Wraps the existing app.state.session_maker with an asyncio.Lock so that
    concurrent writes are serialized.
    """
    base_maker = app.state.session_maker
    db_lock = asyncio.Lock()

    @contextlib.asynccontextmanager
    async def locked_session_maker():
        async with db_lock, base_maker() as session:
            yield session

    app.state.session_maker = locked_session_maker


async def ensure_internal_worker(
    session: AsyncSession,
    email: str,
) -> None:
    """Create or update the internal worker superuser.

    Idempotent — safe to call on every startup. The password is a random
    UUID generated each time — it is never used for login.

    Parameters
    ----------
    session
        Async database session.
    email
        Internal worker email from ``Settings.internal_worker_email``.
    """
    password_helper = PasswordHelper()

    result = await session.exec(
        select(User).where(User.email == email)  # type: ignore[arg-type]
    )
    existing = result.one_or_none()

    hashed = password_helper.hash(str(uuid.uuid4()))

    if existing is None:
        worker = User(
            email=email,
            hashed_password=hashed,
            is_active=True,
            is_superuser=True,
            is_verified=True,
        )
        session.add(worker)
        await session.commit()
        log.info("Created internal worker user: %s", email)
    else:
        existing.hashed_password = hashed
        existing.is_active = True
        existing.is_superuser = True
        existing.is_verified = True
        await session.commit()
        log.debug("Updated internal worker user: %s", email)


async def init_database(
    engine: AsyncEngine | None = None,
    settings: Settings | None = None,
) -> None:
    """Initialize ALL tables for zndraw-fastapi and dependencies.

    Imports all models to register in metadata, then creates tables.
    Idempotent - safe to call multiple times (CREATE TABLE IF NOT EXISTS).

    Parameters
    ----------
    engine : AsyncEngine | None
        Optional engine to use (app startup). If not provided, creates one
        from settings and disposes it after (CLI use).
    settings : Settings | None
        Optional pre-built settings (e.g. from CLI overrides). If not
        provided, a fresh ``Settings()`` is created.
    """
    # Models already imported at module level
    if settings is None:
        settings = Settings()
    auth_settings = AuthSettings()

    # Create engine if not provided (CLI mode) or use existing (app startup)
    own_engine = engine is None
    if own_engine:
        engine = create_engine_for_url(settings.database_url)

    async with engine.begin() as conn:
        await conn.run_sync(SQLModel.metadata.create_all)

    # Create default users
    session_maker = async_sessionmaker(
        engine, class_=AsyncSession, expire_on_commit=False
    )
    async with session_maker() as session:
        await ensure_default_admin(session, auth_settings)
    async with session_maker() as session:
        await ensure_internal_worker(session, settings.internal_worker_email)

    # Seed @internal Job rows for built-in extensions
    await ensure_internal_jobs(_collect_extensions(), session_maker)

    # Seed @internal provider rows. @internal providers are server-owned —
    # they have no Worker row, so we only need the internal worker User.id
    # (used at mint time for JWTs in get_worker_token).
    if settings.filebrowser_enabled:
        from zndraw_joblib.registry import ensure_internal_providers

        async with session_maker() as session:
            result = await session.exec(
                select(User).where(User.email == settings.internal_worker_email)
            )
            internal_user = result.one()
            internal_user_id = internal_user.id

        await ensure_internal_providers(
            list(BUNDLED_PROVIDERS),
            session_maker,
            user_id=internal_user_id,
        )
    else:
        # Feature disabled — clean up any stale @internal provider rows from a
        # previous run that had it enabled.
        from sqlalchemy import delete

        from zndraw_joblib.models import ProviderRecord

        async with session_maker() as session:
            await session.exec(
                delete(ProviderRecord).where(ProviderRecord.room_id == "@internal")
            )
            await session.commit()

    # Only dispose if we created the engine (CLI mode)
    if own_engine:
        await engine.dispose()


@contextlib.asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncIterator[None]:
    """FastAPI lifespan context manager for all application resources."""
    settings = Settings(**app.state.settings_overrides)  # type: ignore[arg-type]
    auth_settings = AuthSettings()
    joblib_settings = JobLibSettings(
        allowed_categories=["modifiers", "selections", "analysis", "filesystem"],
    )

    app.state.settings = settings
    app.state.auth_settings = auth_settings
    app.state.joblib_settings = joblib_settings

    # Database - create engine and session_maker
    engine = create_engine_for_url(settings.database_url)
    app.state.engine = engine
    app.state.session_maker = async_sessionmaker(
        engine, class_=SQLModelAsyncSession, expire_on_commit=False
    )

    # Pre-initialize all app.state attrs consumed by request-time dependencies.
    # Downstream deps access these directly (no getattr fallback); lifespan
    # sets them to None now and overwrites with real values below.
    app.state.internal_worker_user = None
    app.state.internal_provider_registry = None
    app.state.internal_registry = None

    try:
        # Initialize tables on startup (dev mode)
        if settings.init_db_on_startup:
            await init_database(engine=engine, settings=settings)

        # Cache the internal worker User for get_worker_token (avoids the
        # SQLite-lock deadlock from re-querying mid-request when the route
        # already holds a yield-based SessionDep).
        async with app.state.session_maker() as session:
            result = await session.exec(
                select(User).where(User.email == settings.internal_worker_email)  # type: ignore[arg-type]
            )
            app.state.internal_worker_user = result.one_or_none()

        # SQLite locking (if needed)
        if _is_sqlite(settings.database_url):
            _apply_sqlite_locking(app)

        # Redis - auto-start TcpFakeServer if REDIS_URL not configured
        app.state.fake_server = None
        app.state.fake_server_thread = None
        redis_protocol: int | None = None
        if settings.redis_url is None:
            from fakeredis import TcpFakeServer

            port = _get_free_port()
            server_address = ("127.0.0.1", port)
            fake_server = TcpFakeServer(server_address, server_type="redis")
            fake_server_thread = threading.Thread(
                target=fake_server.serve_forever, daemon=True
            )
            fake_server_thread.start()
            app.state.fake_server = fake_server
            app.state.fake_server_thread = fake_server_thread
            redis_url = f"redis://127.0.0.1:{port}"
            redis_protocol = 3  # TcpFakeServer speaks RESP3
        else:
            redis_url = settings.redis_url

        app.state.redis_url = redis_url
        redis_kwargs: dict = {}
        if redis_protocol is not None:
            redis_kwargs["protocol"] = redis_protocol
        app.state.redis = redis_client.from_url(
            redis_url, decode_responses=True, **redis_kwargs
        )

        # Frame storage: room-scoped AsyncBlobIO registry
        from zndraw.storage import FrameStorage

        app.state.frame_storage = FrameStorage(
            uri=settings.storage, redis=app.state.redis
        )

        # Socket.IO with AsyncRedisManager
        client_manager = socketio_lib.AsyncRedisManager(redis_url)
        client_manager.set_server(tsio)
        tsio.manager = client_manager
        tsio.app = app  # Enable DI in socket handlers (resolved at event time)
        app.state.tsio = tsio

        # zndraw-joblib: Socket.IO handlers for job room management
        @tsio.on(JoinJobRoom)
        async def handle_join_job_room(sid: str, data: JoinJobRoom) -> None:
            await tsio.enter_room(sid, f"jobs:{data.job_name}")
            session = await tsio.get_session(sid)
            session["worker_id"] = data.worker_id
            await tsio.save_session(sid, session)

        @tsio.on(LeaveJobRoom)
        async def handle_leave_job_room(sid: str, data: LeaveJobRoom) -> None:
            await tsio.leave_room(sid, f"jobs:{data.job_name}")

        @tsio.on(JoinProviderRoom)
        async def handle_join_provider_room(sid: str, data: JoinProviderRoom) -> None:
            await tsio.enter_room(sid, f"providers:{data.provider_name}")
            session = await tsio.get_session(sid)
            session.setdefault("worker_id", data.worker_id)
            await tsio.save_session(sid, session)

        @tsio.on(LeaveProviderRoom)
        async def handle_leave_provider_room(sid: str, data: LeaveProviderRoom) -> None:
            await tsio.leave_room(sid, f"providers:{data.provider_name}")

        # zndraw-joblib: register internal extensions for TaskIQ dispatch
        broker = ListQueueBroker(redis_url, queue_name=settings.task_queue_name)
        # Internal executor connects back to the same server — resolve
        # 0.0.0.0 (bind-all) to 127.0.0.1 (loopback) for the client URL.
        executor_host = "127.0.0.1" if settings.host == "0.0.0.0" else settings.host
        from zndraw.providers.bootstrap import build_internal_providers_resolver

        executor = InternalExtensionExecutor(
            base_url=f"http://{executor_host}:{settings.port}",
            providers_resolver=build_internal_providers_resolver(settings),
        )

        if settings.init_db_on_startup:
            await register_internal_jobs(
                app,
                broker,
                extensions=_collect_extensions(),
                executor=executor,
                session_factory=app.state.session_maker,
            )
        else:
            # Production: db-init already seeded Job rows — only register
            # broker tasks (no DB writes, avoids unique_job race).
            registry = register_internal_tasks(broker, _collect_extensions(), executor)
            app.state.internal_registry = registry

        # Register @internal provider tasks
        from zndraw.providers.bootstrap import register_filebrowser_providers

        provider_registry = register_filebrowser_providers(
            broker,
            base_url=f"http://{executor_host}:{settings.port}",
            settings=settings,
        )
        app.state.internal_provider_registry = provider_registry

        await broker.startup()

        # zndraw-joblib: wire ResultBackend for provider caching
        # Frame data is too large for Redis — route to the same storage
        # backend used for frame storage.  Everything else stays in Redis.
        from zndraw.result_backends import (
            CompositeResultBackend,
            RedisResultBackend,
            StorageResultBackend,
        )
        from zndraw_joblib.dependencies import (
            get_frame_room_cleanup,
            get_result_backend,
        )

        redis_raw = redis_client.from_url(redis_url, **redis_kwargs)
        redis_backend = RedisResultBackend(
            redis_raw, key_prefix=settings.result_backend_key_prefix
        )
        frame_cache = StorageResultBackend(
            app.state.frame_storage,
            key_prefix=settings.result_backend_key_prefix,
        )
        result_backend = CompositeResultBackend(redis=redis_backend, frames=frame_cache)
        app.state.result_backend = result_backend
        app.dependency_overrides[get_result_backend] = lambda: result_backend

        # zndraw-joblib: frame room cleanup callback
        # Deletes Redis provider frame count keys and emits FramesInvalidate
        # when frame providers are removed (DELETE worker or sweeper).
        async def frame_room_cleanup(room_ids: set[str]) -> None:
            for rid in room_ids:
                await app.state.redis.delete(  # type: ignore[misc]
                    RedisKey.provider_frame_count(rid)
                )
                await tsio.emit(
                    FramesInvalidate(
                        room_id=rid,
                        action="clear",
                        count=0,
                        reason="provider_disconnected",
                    ),
                    room=f"room:{rid}",
                )

        app.dependency_overrides[get_frame_room_cleanup] = lambda: frame_room_cleanup

        # Spawn in-process TaskIQ worker (disabled in Docker — dedicated containers)
        worker_task = None
        if settings.worker_enabled:
            worker_task = asyncio.create_task(run_receiver_task(broker))

        # zndraw-joblib: background sweeper for stale workers
        async def get_session():
            async with app.state.session_maker() as session:
                yield session

        sweeper_task = asyncio.create_task(
            run_sweeper(
                get_session=get_session,
                settings=joblib_settings,
                tsio=tsio,
                on_frame_rooms=frame_room_cleanup,
            )
        )

        # Warm heavy optional imports in a background thread so first
        # request doesn't pay the import cost.  Failures are silently
        # ignored — the route-level import will raise if truly missing.
        def _warm_imports() -> None:
            for mod in ("molify", "rdkit", "rdkit.Chem", "rdkit.Chem.Draw"):
                with contextlib.suppress(ImportError):
                    importlib.import_module(mod)

        warmup_task = asyncio.create_task(asyncio.to_thread(_warm_imports))

        yield

        # Shutdown in reverse order
        if worker_task is not None:
            worker_task.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await worker_task

        await broker.shutdown()

        sweeper_task.cancel()
        with contextlib.suppress(asyncio.CancelledError):
            await sweeper_task

        warmup_task.cancel()
        with contextlib.suppress(asyncio.CancelledError):
            await warmup_task

        await app.state.frame_storage.close()
        await redis_raw.aclose()
        await app.state.redis.aclose()
        if app.state.fake_server is not None:
            app.state.fake_server.shutdown()
            app.state.fake_server.server_close()
            if app.state.fake_server_thread is not None:
                app.state.fake_server_thread.join(timeout=1.0)
    finally:
        await engine.dispose()
