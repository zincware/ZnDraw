"""Database models and session management."""

import logging
import uuid
from collections.abc import AsyncIterator
from datetime import datetime
from typing import Annotated

from fastapi import Depends, Request
from fastapi_users.db import SQLAlchemyBaseUserTableUUID, SQLAlchemyUserDatabase
from fastapi_users.password import PasswordHelper
from sqlalchemy import select
from sqlalchemy.ext.asyncio import (
    AsyncEngine,
    AsyncSession,
    async_sessionmaker,
    create_async_engine,
)
from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.pool import NullPool, StaticPool
from sqlmodel import Field, SQLModel

from zndraw_auth.settings import AuthSettings

log = logging.getLogger(__name__)


class Base(DeclarativeBase):
    """SQLAlchemy declarative base, sharing metadata with SQLModel."""

    metadata = SQLModel.metadata


class User(SQLAlchemyBaseUserTableUUID, Base):
    """User model for authentication.

    Inherits from fastapi-users base which provides:
    - id: UUID (primary key)
    - email: str (unique, indexed)
    - hashed_password: str
    - is_active: bool (default True)
    - is_superuser: bool (default False)
    - is_verified: bool (default False)
    """

    pass


class CLILoginChallenge(SQLModel, table=True):
    """Challenge for device-code style CLI login flow.

    Lifecycle: pending -> approved -> redeemed
    On redeem: token and secret are nulled, row kept for audit.
    """

    id: int | None = Field(default=None, primary_key=True)
    code: str = Field(index=True, unique=True)
    secret: str | None = None
    status: str = "pending"
    token: str | None = None
    user_id: uuid.UUID | None = Field(default=None, foreign_key="user.id")
    created_at: datetime
    expires_at: datetime


def create_engine_for_url(database_url: str) -> AsyncEngine:
    """Create engine with appropriate connection pooling.

    Strategy:
    - In-memory SQLite: StaticPool (single shared connection)
    - File SQLite: NullPool (connection per checkout, avoids locks)
    - PostgreSQL: QueuePool (default connection pool)
    """
    if database_url == "sqlite+aiosqlite://":
        return create_async_engine(
            database_url,
            connect_args={"check_same_thread": False},
            poolclass=StaticPool,
        )
    if database_url.startswith("sqlite"):
        return create_async_engine(database_url, poolclass=NullPool)
    return create_async_engine(database_url)


def get_engine(request: Request) -> AsyncEngine:
    """Retrieve engine from app.state.

    Override point #1 (advanced use cases).
    """
    return request.app.state.engine


def get_session_maker(request: Request) -> async_sessionmaker[AsyncSession]:
    """Retrieve session maker from app.state.

    Override point #2 (PRIMARY for tests).
    Tests override by setting app.state.session_maker directly.
    Can also be overridden via app.dependency_overrides[get_session_maker].

    Returns factory, not session, because:
    - Long-polling needs multiple sessions per request
    - Socket.IO needs session per event
    - TaskIQ needs session per task
    """
    return request.app.state.session_maker


async def get_session(
    session_maker: Annotated[
        async_sessionmaker[AsyncSession], Depends(get_session_maker)
    ],
) -> AsyncIterator[AsyncSession]:
    """Yield request-scoped session.

    Override point #3 (rare - specific session mocking).

    Session lifecycle:
    - Created at request time
    - Must be explicitly committed by caller
    - Rolled back on exception
    - Closed after request
    """
    async with session_maker() as session:
        yield session


# Type alias for convenience
SessionDep = Annotated[AsyncSession, Depends(get_session)]


async def ensure_default_admin(
    session: AsyncSession,
    settings: AuthSettings,
) -> None:
    """Create or promote default admin user.

    Called by host app during initialization with a session.
    Idempotent - safe to call multiple times.

    If default_admin_email is None, runs in dev mode (all users are superusers).

    Parameters
    ----------
    session : AsyncSession
        Database session to use for operations.
    settings : AuthSettings
        Authentication settings with admin credentials.
    """
    if settings.default_admin_email is None:
        log.info("No default admin configured - running in dev mode")
        return

    if settings.default_admin_password is None:
        log.warning(
            "DEFAULT_ADMIN_EMAIL is set but DEFAULT_ADMIN_PASSWORD is not - "
            "skipping admin creation"
        )
        return

    password_helper = PasswordHelper()

    # Check if user exists
    result = await session.execute(
        select(User).where(User.email == settings.default_admin_email)  # type: ignore[arg-type]
    )
    existing = result.scalar_one_or_none()

    if existing is None:
        # Create admin user
        hashed = password_helper.hash(
            settings.default_admin_password.get_secret_value()
        )
        admin = User(
            email=settings.default_admin_email,
            hashed_password=hashed,
            is_active=True,
            is_superuser=True,
            is_verified=True,
        )
        session.add(admin)
        await session.commit()
        log.info(f"Created default admin user: {settings.default_admin_email}")
    elif not existing.is_superuser:
        # Promote to superuser
        existing.is_superuser = True
        await session.commit()
        log.info(f"Promoted user to superuser: {settings.default_admin_email}")
    else:
        log.debug(f"Default admin already exists: {settings.default_admin_email}")


async def get_user_db(
    session: SessionDep,
) -> AsyncIterator[SQLAlchemyUserDatabase[User, uuid.UUID]]:
    """FastAPI dependency that yields the user database adapter."""
    yield SQLAlchemyUserDatabase[User, uuid.UUID](session, User)
