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

_user_counter = 0


def _next_display_name(prefix: str = "test-user") -> str:
    global _user_counter
    _user_counter += 1
    return f"{prefix}-{_user_counter:04d}"


async def _make_user(
    session: AsyncSession,
    email: str,
    *,
    superuser: bool = False,
    display_name: str | None = None,
) -> User:
    u = User(
        email=email,
        hashed_password="x",
        is_superuser=superuser,
        display_name=display_name or _next_display_name("acc-user"),
    )
    session.add(u)
    await session.commit()
    await session.refresh(u)
    return u


@pytest.mark.asyncio
async def test_readable_private_room_raises_404_for_stranger(
    session: AsyncSession,
) -> None:
    owner = await _make_user(session, "o@x", display_name="acc-owner-priv1")
    stranger = await _make_user(session, "s@x", display_name="acc-stranger-1")
    room = Room(
        room_name="r1",
        owner_user_id=owner.id,
        created_by_id=owner.id,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    await session.commit()

    with pytest.raises(ProblemError) as exc:
        await get_readable_room(
            session=session,
            current_user=stranger,
            share=None,
            owner=owner.display_name,
            room_name="r1",
        )
    assert exc.value.problem.status == 404


@pytest.mark.asyncio
async def test_readable_private_room_visible_to_owner(
    session: AsyncSession,
) -> None:
    owner = await _make_user(session, "o2@x", display_name="acc-owner-priv2")
    room = Room(
        room_name="r2",
        owner_user_id=owner.id,
        created_by_id=owner.id,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    await session.commit()

    ctx = await get_readable_room(
        session=session,
        current_user=owner,
        share=None,
        owner=owner.display_name,
        room_name="r2",
    )
    assert ctx.room.room_name == "r2"
    assert ctx.room.owner_user_id == owner.id
    assert ctx.share is None
    assert ctx.group_role is None


@pytest.mark.asyncio
async def test_readable_group_room_sees_role(session: AsyncSession) -> None:
    admin = await _make_user(session, "ga@x", display_name="acc-grp-admin1")
    member = await _make_user(session, "gm@x", display_name="acc-grp-member1")
    g = Group(name="acc-group-t1", created_by_id=admin.id)
    session.add(g)
    await session.commit()
    session.add(
        GroupMembership(group_id=g.id, user_id=member.id, role=GroupRole.VIEWER)
    )
    room = Room(
        room_name="r3",
        owner_group_id=g.id,
        created_by_id=admin.id,
        visibility=Visibility.GROUP,
    )
    session.add(room)
    await session.commit()

    ctx = await get_readable_room(
        session=session,
        current_user=member,
        share=None,
        owner=g.name,
        room_name="r3",
    )
    assert ctx.group_role is GroupRole.VIEWER


@pytest.mark.asyncio
async def test_editable_blocked_for_viewer(session: AsyncSession) -> None:
    admin = await _make_user(session, "ea@x", display_name="acc-edit-admin")
    viewer = await _make_user(session, "ev@x", display_name="acc-edit-viewer")
    g = Group(name="acc-group-t2", created_by_id=admin.id)
    session.add(g)
    await session.commit()
    session.add(
        GroupMembership(group_id=g.id, user_id=viewer.id, role=GroupRole.VIEWER)
    )
    room = Room(
        room_name="r4",
        owner_group_id=g.id,
        created_by_id=admin.id,
        visibility=Visibility.GROUP,
    )
    session.add(room)
    await session.commit()

    # Viewer can read
    ctx = await get_readable_room(
        session=session,
        current_user=viewer,
        share=None,
        owner=g.name,
        room_name="r4",
    )
    # but cannot edit — composite should raise Forbidden
    with pytest.raises(ProblemError) as exc:
        await get_editable_room(ctx=ctx, current_user=viewer)
    assert exc.value.problem.status == 403


@pytest.mark.asyncio
async def test_manageable_blocked_for_member(session: AsyncSession) -> None:
    admin = await _make_user(session, "ma@x", display_name="acc-mng-admin")
    member = await _make_user(session, "mm@x", display_name="acc-mng-member")
    g = Group(name="acc-group-t3", created_by_id=admin.id)
    session.add(g)
    await session.commit()
    session.add(
        GroupMembership(group_id=g.id, user_id=member.id, role=GroupRole.MEMBER)
    )
    room = Room(
        room_name="r5",
        owner_group_id=g.id,
        created_by_id=admin.id,
        visibility=Visibility.GROUP,
    )
    session.add(room)
    await session.commit()

    ctx = await get_readable_room(
        session=session,
        current_user=member,
        share=None,
        owner=g.name,
        room_name="r5",
    )
    # can edit, but not manage
    await get_editable_room(ctx=ctx, current_user=member)
    with pytest.raises(ProblemError) as exc:
        await get_manageable_room(ctx=ctx, current_user=member)
    assert exc.value.problem.status == 403


@pytest.mark.asyncio
async def test_share_view_grants_read_not_edit(session: AsyncSession) -> None:
    owner = await _make_user(session, "so@x", display_name="acc-share-owner")
    stranger = await _make_user(session, "ss@x", display_name="acc-share-stranger")
    room = Room(
        room_name="r6",
        owner_user_id=owner.id,
        created_by_id=owner.id,
        visibility=Visibility.PRIVATE,
    )
    session.add(room)
    await session.commit()

    share = ShareContext(room_id=room.id, access=ShareAccess.VIEW)
    ctx = await get_readable_room(
        session=session,
        current_user=stranger,
        share=share,
        owner=owner.display_name,
        room_name="r6",
    )
    assert ctx.share is share

    with pytest.raises(ProblemError) as exc:
        await get_editable_room(ctx=ctx, current_user=stranger)
    assert exc.value.problem.status == 403
