"""Resolve a polymorphic owner_id to either a user or a group."""

from __future__ import annotations

from typing import TYPE_CHECKING
from uuid import uuid4

import pytest_asyncio

from zndraw.dependencies import OwnerKind, resolve_owner
from zndraw.models import Group
from zndraw_auth import User

if TYPE_CHECKING:
    from sqlmodel.ext.asyncio.session import AsyncSession


@pytest_asyncio.fixture
async def alice(session: AsyncSession) -> User:
    user = User(email="alice@example.com", hashed_password="x")
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user


async def test_resolve_user(session: AsyncSession, alice: User) -> None:
    result = await resolve_owner(session, alice.id)
    assert result is not None
    kind, label = result
    assert kind == OwnerKind.USER
    assert label == "alice@example.com"


async def test_resolve_group(session: AsyncSession, alice: User) -> None:
    group = Group(name="researchers", created_by_id=alice.id)
    session.add(group)
    await session.commit()
    await session.refresh(group)

    result = await resolve_owner(session, group.id)
    assert result is not None
    kind, label = result
    assert kind == OwnerKind.GROUP
    assert label == "researchers"


async def test_resolve_unknown(session: AsyncSession) -> None:
    assert await resolve_owner(session, uuid4()) is None
