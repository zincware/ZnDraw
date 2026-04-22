"""Unit tests for the AccessReadDep/AccessEditDep/AccessManageDep composites.

Uses the session fixture + direct calls to the async composites, bypassing
FastAPI's DI machinery (deferred to Tasks 10/11 integration tests).
"""

import pytest
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import GroupRole, ShareAccess, ShareContext, Visibility
from zndraw.dependencies import (
    get_editable_room,
    get_manageable_room,
    get_readable_room,
)
from zndraw.exceptions import ProblemError
from zndraw.models import Group, GroupMembership, Room
from zndraw_auth import User


async def _make_user(
    session: AsyncSession, email: str, *, superuser: bool = False
) -> User:
    u = User(email=email, hashed_password="x", is_superuser=superuser)
    session.add(u)
    await session.commit()
    await session.refresh(u)
    return u


@pytest.mark.asyncio
async def test_readable_private_room_raises_404_for_stranger(
    session: AsyncSession,
) -> None:
    owner = await _make_user(session, "o@x")
    stranger = await _make_user(session, "s@x")
    room = Room(
        id="r1",
        owner_user_id=owner.id,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    await session.commit()

    with pytest.raises(ProblemError) as exc:
        await get_readable_room(
            session=session,
            current_user=stranger,
            share=None,
            room_id="r1",
        )
    assert exc.value.problem.status == 404


@pytest.mark.asyncio
async def test_readable_private_room_visible_to_owner(
    session: AsyncSession,
) -> None:
    owner = await _make_user(session, "o2@x")
    room = Room(
        id="r2",
        owner_user_id=owner.id,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    await session.commit()

    ctx = await get_readable_room(
        session=session,
        current_user=owner,
        share=None,
        room_id="r2",
    )
    assert ctx.room.id == "r2"
    assert ctx.share is None
    assert ctx.group_role is None


@pytest.mark.asyncio
async def test_readable_group_room_sees_role(session: AsyncSession) -> None:
    admin = await _make_user(session, "ga@x")
    member = await _make_user(session, "gm@x")
    g = Group(name="t1", created_by_id=admin.id)
    session.add(g)
    await session.commit()
    session.add(
        GroupMembership(group_id=g.id, user_id=member.id, role=GroupRole.VIEWER)
    )
    room = Room(id="r3", owner_group_id=g.id, visibility=Visibility.GROUP)
    session.add(room)
    await session.commit()

    ctx = await get_readable_room(
        session=session,
        current_user=member,
        share=None,
        room_id="r3",
    )
    assert ctx.group_role is GroupRole.VIEWER


@pytest.mark.asyncio
async def test_editable_blocked_for_viewer(session: AsyncSession) -> None:
    admin = await _make_user(session, "ea@x")
    viewer = await _make_user(session, "ev@x")
    g = Group(name="t2", created_by_id=admin.id)
    session.add(g)
    await session.commit()
    session.add(
        GroupMembership(group_id=g.id, user_id=viewer.id, role=GroupRole.VIEWER)
    )
    room = Room(id="r4", owner_group_id=g.id, visibility=Visibility.GROUP)
    session.add(room)
    await session.commit()

    # Viewer can read
    ctx = await get_readable_room(
        session=session, current_user=viewer, share=None, room_id="r4"
    )
    # but cannot edit — composite should raise Forbidden
    with pytest.raises(ProblemError) as exc:
        await get_editable_room(ctx=ctx, current_user=viewer)
    assert exc.value.problem.status == 403


@pytest.mark.asyncio
async def test_manageable_blocked_for_member(session: AsyncSession) -> None:
    admin = await _make_user(session, "ma@x")
    member = await _make_user(session, "mm@x")
    g = Group(name="t3", created_by_id=admin.id)
    session.add(g)
    await session.commit()
    session.add(
        GroupMembership(group_id=g.id, user_id=member.id, role=GroupRole.MEMBER)
    )
    room = Room(id="r5", owner_group_id=g.id, visibility=Visibility.GROUP)
    session.add(room)
    await session.commit()

    ctx = await get_readable_room(
        session=session, current_user=member, share=None, room_id="r5"
    )
    # can edit, but not manage
    await get_editable_room(ctx=ctx, current_user=member)
    with pytest.raises(ProblemError) as exc:
        await get_manageable_room(ctx=ctx, current_user=member)
    assert exc.value.problem.status == 403


@pytest.mark.asyncio
async def test_share_view_grants_read_not_edit(session: AsyncSession) -> None:
    owner = await _make_user(session, "so@x")
    stranger = await _make_user(session, "ss@x")
    room = Room(id="r6", owner_user_id=owner.id, visibility=Visibility.PRIVATE)
    session.add(room)
    await session.commit()

    share = ShareContext(room_id="r6", access=ShareAccess.VIEW)
    ctx = await get_readable_room(
        session=session, current_user=stranger, share=share, room_id="r6"
    )
    assert ctx.share is share

    with pytest.raises(ProblemError) as exc:
        await get_editable_room(ctx=ctx, current_user=stranger)
    assert exc.value.problem.status == 403
