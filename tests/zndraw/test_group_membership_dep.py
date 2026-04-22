"""Integration test for group membership helpers."""
from uuid import uuid4

import pytest
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import GroupRole
from zndraw.dependencies import fetch_group_role, fetch_my_group_ids
from zndraw.models import Group, GroupMembership
from zndraw_auth import User


@pytest.mark.asyncio
async def test_fetch_my_group_ids(session: AsyncSession) -> None:
    u = User(email="a@x", hashed_password="x")
    session.add(u)
    await session.commit()
    g1 = Group(name="g1", created_by_id=u.id)
    g2 = Group(name="g2", created_by_id=u.id)
    session.add_all([g1, g2])
    await session.commit()
    session.add_all([
        GroupMembership(group_id=g1.id, user_id=u.id, role=GroupRole.MEMBER),
        GroupMembership(group_id=g2.id, user_id=u.id, role=GroupRole.VIEWER),
    ])
    await session.commit()

    ids = await fetch_my_group_ids(session, u.id)
    assert set(ids) == {g1.id, g2.id}


@pytest.mark.asyncio
async def test_fetch_group_role(session: AsyncSession) -> None:
    u = User(email="b@x", hashed_password="x")
    session.add(u)
    await session.commit()
    g = Group(name="g3", created_by_id=u.id)
    session.add(g)
    await session.commit()
    session.add(GroupMembership(group_id=g.id, user_id=u.id, role=GroupRole.ADMIN))
    await session.commit()

    role = await fetch_group_role(session, u.id, g.id)
    assert role is GroupRole.ADMIN

    other = uuid4()
    assert await fetch_group_role(session, other, g.id) is None
