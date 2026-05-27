"""CHECK-constraint tests for Room ownership invariants."""

from uuid import uuid4

import pytest
from sqlalchemy.exc import IntegrityError
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import Visibility
from zndraw.models import Group, Room


@pytest.mark.asyncio
async def test_room_requires_exactly_one_owner(session: AsyncSession) -> None:
    """Room with both owner_user_id and owner_group_id set must fail."""
    owner_uid = uuid4()
    g = Group(name="g1", created_by_id=owner_uid)
    session.add(g)
    await session.commit()

    room = Room(
        id=str(uuid4()),
        owner_user_id=owner_uid,
        owner_group_id=g.id,
        visibility=Visibility.PUBLIC,
    )
    session.add(room)
    with pytest.raises(IntegrityError):
        await session.commit()


@pytest.mark.asyncio
async def test_room_requires_at_least_one_owner_unless_public(
    session: AsyncSession,
) -> None:
    """Room with neither owner must fail (XOR violation)."""
    room = Room(
        id=str(uuid4()),
        owner_user_id=None,
        owner_group_id=None,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    with pytest.raises(IntegrityError):
        await session.commit()


@pytest.mark.asyncio
async def test_visibility_private_requires_user_owner(session: AsyncSession) -> None:
    """PRIVATE visibility with only a group owner must fail."""
    uid = uuid4()
    g = Group(name="g2", created_by_id=uid)
    session.add(g)
    await session.commit()

    room = Room(
        id=str(uuid4()),
        owner_group_id=g.id,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    with pytest.raises(IntegrityError):
        await session.commit()


@pytest.mark.asyncio
async def test_visibility_group_requires_group_owner(session: AsyncSession) -> None:
    """GROUP visibility with only a user owner must fail."""
    uid = uuid4()
    room = Room(
        id=str(uuid4()),
        owner_user_id=uid,
        visibility=Visibility.GROUP,
    )
    session.add(room)
    with pytest.raises(IntegrityError):
        await session.commit()
