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
    user = User(
        email="alice@example.com",
        hashed_password="x",
        display_name="alice-the-explorer",
    )
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user


async def test_resolve_user(session: AsyncSession, alice: User) -> None:
    result = await resolve_owner(session, alice.id)
    assert result is not None
    kind, label = result
    assert kind == OwnerKind.USER
    assert label == "alice-the-explorer"


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


async def test_resolve_owner_returns_display_name_for_user(
    session: AsyncSession,
) -> None:
    user = User(
        email="ada@example.com",
        hashed_password="x",
        display_name="ada-lovelace-coder",
    )
    session.add(user)
    await session.commit()
    await session.refresh(user)

    result = await resolve_owner(session, user.id)
    assert result is not None
    kind, label = result
    assert kind == OwnerKind.USER
    assert label == "ada-lovelace-coder"


async def test_resolve_owner_returns_group_name_unchanged(
    session: AsyncSession,
) -> None:
    creator = User(
        email="creator@example.com",
        hashed_password="x",
        display_name="creator-of-groups",
    )
    session.add(creator)
    await session.commit()
    await session.refresh(creator)

    group = Group(name="my-group", created_by_id=creator.id)
    session.add(group)
    await session.commit()
    await session.refresh(group)

    result = await resolve_owner(session, group.id)
    assert result is not None
    kind, label = result
    assert kind == OwnerKind.GROUP
    assert label == "my-group"


async def test_get_owner_uuid_from_segment_resolves_user(
    session: AsyncSession,
) -> None:
    from zndraw.dependencies import get_owner_uuid_from_segment

    user = User(
        email="seg@example.com",
        hashed_password="x",
        display_name="seg-test-user",
    )
    session.add(user)
    await session.commit()
    await session.refresh(user)

    assert (
        await get_owner_uuid_from_segment(session, "seg-test-user")
    ) == user.id


async def test_get_owner_uuid_from_segment_resolves_group(
    session: AsyncSession,
) -> None:
    from zndraw.dependencies import get_owner_uuid_from_segment

    creator = User(
        email="gcreator@example.com",
        hashed_password="x",
        display_name="g-creator-display",
    )
    session.add(creator)
    await session.commit()
    await session.refresh(creator)

    group = Group(name="visible-group", created_by_id=creator.id)
    session.add(group)
    await session.commit()
    await session.refresh(group)

    assert (
        await get_owner_uuid_from_segment(session, "visible-group")
    ) == group.id


async def test_get_owner_uuid_from_segment_unknown_raises_user_not_found(
    session: AsyncSession,
) -> None:
    import pytest

    from zndraw.dependencies import get_owner_uuid_from_segment
    from zndraw.exceptions import ProblemError

    with pytest.raises(ProblemError) as excinfo:
        await get_owner_uuid_from_segment(session, "no-such-owner-here")
    assert excinfo.value.problem.status == 404
