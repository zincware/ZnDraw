"""Test the is_guest column on User."""

import pytest
from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw_auth import User


@pytest.mark.asyncio
async def test_user_is_guest_defaults_false() -> None:
    """A freshly-created User has is_guest=False unless explicitly set."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    async with engine.begin() as conn:
        await conn.run_sync(SQLModel.metadata.create_all)

    session_maker = async_sessionmaker(
        engine, class_=AsyncSession, expire_on_commit=False
    )
    async with session_maker() as session:
        user = User(email="u@test", hashed_password="x")
        session.add(user)
        await session.commit()
        await session.refresh(user)
        assert user.is_guest is False

    await engine.dispose()


@pytest.mark.asyncio
async def test_user_is_guest_can_be_true() -> None:
    """A guest user can be created with is_guest=True."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    async with engine.begin() as conn:
        await conn.run_sync(SQLModel.metadata.create_all)

    session_maker = async_sessionmaker(
        engine, class_=AsyncSession, expire_on_commit=False
    )
    async with session_maker() as session:
        user = User(email="g@test", hashed_password="x", is_guest=True)
        session.add(user)
        await session.commit()
        await session.refresh(user)
        assert user.is_guest is True

    await engine.dispose()
