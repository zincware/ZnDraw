"""Async database setup and application lifespan management.

Uses app.state pattern to store resources:
- redis: Redis async client
- frame_storage: Frame storage backend (InMemory, LMDB, etc.)
- fake_server: TcpFakeServer instance (when REDIS_URL not configured)
- tsio: zndraw-socketio typed wrapper
"""

import asyncio
import contextlib
import importlib
import logging
import socket
import threading
from collections.abc import AsyncIterator

import redis.asyncio as redis_client
import socketio as socketio_lib
import zndraw_joblib.models  # noqa: F401 - registers Job, Worker, Task
from fastapi import FastAPI
from fastapi_users.password import PasswordHelper
from pydantic import SecretStr
from sqlalchemy.ext.asyncio import AsyncEngine, AsyncSession, async_sessionmaker
from sqlmodel import SQLModel, select
from sqlmodel.ext.asyncio.session import AsyncSession as SQLModelAsyncSession
from taskiq.api import run_receiver_task
from taskiq_redis import ListQueueBroker
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
from zndraw_joblib.settings import JobLibSettings

import zndraw.models  # noqa: F401 - registers Room, Message, etc.
from zndraw.config import (
    LMDBStorage as LMDBStorageConfig,
    MemoryStorage,
    MongoDBStorage,
    Settings,
    StorageConfig,
)
from zndraw.executor import InternalExtensionExecutor
from zndraw.extensions.analysis import analysis
from zndraw.extensions.modifiers import modifiers
from zndraw.extensions.selections import selections
from zndraw.socketio import tsio
from zndraw.storage import InMemoryStorage, LMDBStorage as LMDBStorageBackend
from zndraw.storage.base import StorageBackend


def _get_free_port() -> int:
    """Find a free port on localhost using socket binding."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _create_storage_backend(config: StorageConfig) -> StorageBackend:
    """Create a storage backend based on configuration."""
    if isinstance(config, MemoryStorage):
        return InMemoryStorage()
    if isinstance(config, LMDBStorageConfig):
        return LMDBStorageBackend(path=config.path, map_size=config.map_size)
    if isinstance(config, MongoDBStorage):
        raise NotImplementedError("MongoDB storage is not yet implemented")
    raise ValueError(f"Unknown storage config type: {type(config)}")


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


WORKER_EMAIL = "worker@internal.user"


async def ensure_internal_worker(
    session: AsyncSession,
    password: SecretStr,
) -> None:
    """Create or update the internal worker superuser.

    Idempotent — safe to call on every startup.

    Parameters
    ----------
    session
        Async database session.
    password
        Worker password from ``Settings.worker_password``.
    """
    password_helper = PasswordHelper()

    result = await session.execute(
        select(User).where(User.email == WORKER_EMAIL)  # type: ignore[arg-type]
    )
    existing = result.scalar_one_or_none()

    hashed = password_helper.hash(password.get_secret_value())

    if existing is None:
        worker = User(
            email=WORKER_EMAIL,
            hashed_password=hashed,
            is_active=True,
            is_superuser=True,
            is_verified=True,
        )
        session.add(worker)
        await session.commit()
        log.info("Created internal worker user: %s", WORKER_EMAIL)
    else:
        existing.hashed_password = hashed
        existing.is_superuser = True
        await session.commit()
        log.debug("Updated internal worker user: %s", WORKER_EMAIL)


async def init_database(engine: AsyncEngine | None = None) -> None:
    """Initialize ALL tables for zndraw-fastapi and dependencies.

    Imports all models to register in metadata, then creates tables.
    Idempotent - safe to call multiple times (CREATE TABLE IF NOT EXISTS).

    Parameters
    ----------
    engine : AsyncEngine | None
        Optional engine to use (app startup). If not provided, creates one
        from settings and disposes it after (CLI use).
    """
    # Models already imported at module level
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
        await ensure_internal_worker(session, settings.worker_password)

    # Only dispose if we created the engine (CLI mode)
    if own_engine:
        await engine.dispose()


@contextlib.asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncIterator[None]:
    """FastAPI lifespan context manager for all application resources."""
    # Create all settings and store in app.state
    settings = Settings()
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

    try:
        # Initialize tables on startup (dev mode)
        if settings.init_db_on_startup:
            await init_database(engine=engine)

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

        # Frame storage (wrapped in StorageRouter for virtual mount support)
        # The underlying backend is shared with the provider cache so frame
        # provider results automatically use the same storage engine.
        from zndraw.storage.router import StorageRouter

        default_storage = _create_storage_backend(settings.storage)
        app.state.frame_storage = StorageRouter(
            default=default_storage,
            redis=app.state.redis,
        )

        # Socket.IO with AsyncRedisManager
        client_manager = socketio_lib.AsyncRedisManager(redis_url)
        client_manager.set_server(tsio)
        tsio.manager = client_manager
        tsio.manager_initialized = True
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
        broker = ListQueueBroker(redis_url)
        # Internal executor connects back to the same server — resolve
        # 0.0.0.0 (bind-all) to 127.0.0.1 (loopback) for the client URL.
        executor_host = "127.0.0.1" if settings.host == "0.0.0.0" else settings.host
        executor = InternalExtensionExecutor(
            base_url=f"http://{executor_host}:{settings.port}",
            worker_email=WORKER_EMAIL,
            worker_password=settings.worker_password,
        )

        await register_internal_jobs(
            app,
            broker,
            extensions=_collect_extensions(),
            executor=executor,
            session_factory=app.state.session_maker,
        )
        await broker.startup()

        # zndraw-joblib: wire ResultBackend for provider caching
        # Frame data is too large for Redis — route to the same storage
        # backend used for frame storage.  Everything else stays in Redis.
        from zndraw_joblib.dependencies import get_result_backend

        from zndraw.result_backends import (
            CompositeResultBackend,
            RedisResultBackend,
            StorageResultBackend,
        )

        redis_raw = redis_client.from_url(redis_url, **redis_kwargs)
        redis_backend = RedisResultBackend(redis_raw)
        frame_cache = StorageResultBackend(default_storage)
        result_backend = CompositeResultBackend(redis=redis_backend, frames=frame_cache)
        app.state.result_backend = result_backend
        app.dependency_overrides[get_result_backend] = lambda: result_backend

        # Spawn in-process TaskIQ worker
        worker_task = asyncio.create_task(run_receiver_task(broker))

        # zndraw-joblib: background sweeper for stale workers
        async def get_session():
            async with app.state.session_maker() as session:
                yield session

        sweeper_task = asyncio.create_task(
            run_sweeper(get_session=get_session, settings=joblib_settings, tsio=tsio)
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
