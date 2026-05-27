"""Group CRUD and membership endpoints."""

from uuid import UUID

from fastapi import APIRouter, status
from sqlalchemy.exc import IntegrityError
from sqlmodel import select

from zndraw.access import GroupRole
from zndraw.dependencies import (
    CurrentUserDep,
    SessionDep,
    fetch_group_role,
)
from zndraw.exceptions import (
    GroupHasRooms,
    GroupNameTaken,
    GroupNotFound,
    LastGroupAdmin,
    NotGroupAdmin,
    NotGroupMember,
    UserNotFound,
    problem_responses,
)
from zndraw.models import Group, GroupMembership, Room
from zndraw.schemas import (
    CollectionResponse,
    GroupCreate,
    GroupMemberCreateRequest,
    GroupMemberPatchRequest,
    GroupMemberResponse,
    GroupPatchRequest,
    GroupResponse,
)
from zndraw_auth import User

router = APIRouter(prefix="/v1/groups", tags=["groups"])


async def _require_admin(session: SessionDep, user_id: UUID, group_id: UUID) -> None:
    """Raise NotGroupMember or NotGroupAdmin if user is not a group admin.

    Parameters
    ----------
    session
        Async database session.
    user_id
        The user to check.
    group_id
        The group to check admin role in.
    """
    role = await fetch_group_role(session, user_id, group_id)
    if role is None:
        raise NotGroupMember.exception("You are not a member of this group")
    if role != GroupRole.ADMIN:
        raise NotGroupAdmin.exception("Admin role required")


async def _count_admins(session: SessionDep, group_id: UUID) -> int:
    """Return the number of admins in the given group.

    Parameters
    ----------
    session
        Async database session.
    group_id
        The group to count admins in.
    """
    result = await session.exec(
        select(GroupMembership).where(
            GroupMembership.group_id == group_id,
            GroupMembership.role == GroupRole.ADMIN,
        )
    )
    return len(list(result.all()))


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(GroupNameTaken),
)
async def create_group(
    session: SessionDep,
    current_user: CurrentUserDep,
    payload: GroupCreate,
) -> GroupResponse:
    """Create a new group and make the creator an admin.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user creating the group.
    payload
        Group creation request body.
    """
    group = Group(
        name=payload.name,
        description=payload.description,
        created_by_id=current_user.id,
    )
    session.add(group)
    try:
        await session.flush()
    except IntegrityError as exc:
        await session.rollback()
        raise GroupNameTaken.exception(
            f"Group name '{payload.name}' is already in use"
        ) from exc

    session.add(
        GroupMembership(
            group_id=group.id, user_id=current_user.id, role=GroupRole.ADMIN
        )
    )
    await session.commit()
    await session.refresh(group)
    return GroupResponse.model_validate(
        {**group.model_dump(), "my_role": GroupRole.ADMIN}
    )


@router.get("")
async def list_my_groups(
    session: SessionDep, current_user: CurrentUserDep
) -> CollectionResponse[GroupResponse]:
    """List all groups the current user is a member of.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user.
    """
    result = await session.exec(
        select(Group, GroupMembership.role)
        .join(GroupMembership, GroupMembership.group_id == Group.id)
        .where(GroupMembership.user_id == current_user.id)
    )
    items = [
        GroupResponse.model_validate({**group.model_dump(), "my_role": role})
        for group, role in result.all()
    ]
    return CollectionResponse(items=items)


@router.get(
    "/{group_id}",
    responses=problem_responses(GroupNotFound, NotGroupMember),
)
async def get_group(
    session: SessionDep, current_user: CurrentUserDep, group_id: UUID
) -> GroupResponse:
    """Get a group by ID.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user.
    group_id
        UUID of the group to retrieve.
    """
    group = await session.get(Group, group_id)
    if group is None:
        raise GroupNotFound.exception("Group not found")
    role = await fetch_group_role(session, current_user.id, group_id)
    if role is None and not current_user.is_superuser:
        raise NotGroupMember.exception("You are not a member of this group")
    return GroupResponse.model_validate({**group.model_dump(), "my_role": role})


@router.patch(
    "/{group_id}",
    responses=problem_responses(GroupNotFound, NotGroupAdmin, GroupNameTaken),
)
async def update_group(
    session: SessionDep,
    current_user: CurrentUserDep,
    group_id: UUID,
    payload: GroupPatchRequest,
) -> GroupResponse:
    """Update group metadata (name or description).

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user (must be group admin or superuser).
    group_id
        UUID of the group to update.
    payload
        Partial update request body.
    """
    group = await session.get(Group, group_id)
    if group is None:
        raise GroupNotFound.exception("Group not found")
    if not current_user.is_superuser:
        await _require_admin(session, current_user.id, group_id)
    if payload.name is not None:
        group.name = payload.name
    if payload.description is not None:
        group.description = payload.description
    try:
        await session.commit()
    except IntegrityError as exc:
        await session.rollback()
        raise GroupNameTaken.exception("Group name already in use") from exc
    role = await fetch_group_role(session, current_user.id, group_id)
    await session.refresh(group)
    return GroupResponse.model_validate({**group.model_dump(), "my_role": role})


@router.delete(
    "/{group_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    responses=problem_responses(GroupNotFound, NotGroupAdmin, GroupHasRooms),
)
async def delete_group(
    session: SessionDep, current_user: CurrentUserDep, group_id: UUID
) -> None:
    """Delete a group if it owns no rooms.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user (must be group admin or superuser).
    group_id
        UUID of the group to delete.
    """
    group = await session.get(Group, group_id)
    if group is None:
        raise GroupNotFound.exception("Group not found")
    if not current_user.is_superuser:
        await _require_admin(session, current_user.id, group_id)
    rooms = await session.exec(select(Room.id).where(Room.owner_group_id == group_id))
    if rooms.first() is not None:
        raise GroupHasRooms.exception(
            "Group still owns one or more rooms; reassign or delete them first"
        )
    memberships = await session.exec(
        select(GroupMembership).where(GroupMembership.group_id == group_id)
    )
    for m in memberships.all():
        await session.delete(m)
    await session.delete(group)
    await session.commit()


@router.get(
    "/{group_id}/members",
    responses=problem_responses(GroupNotFound, NotGroupMember),
)
async def list_members(
    session: SessionDep, current_user: CurrentUserDep, group_id: UUID
) -> CollectionResponse[GroupMemberResponse]:
    """List all members of a group.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user (must be a group member or superuser).
    group_id
        UUID of the group to list members for.
    """
    group = await session.get(Group, group_id)
    if group is None:
        raise GroupNotFound.exception("Group not found")
    if not current_user.is_superuser:
        role = await fetch_group_role(session, current_user.id, group_id)
        if role is None:
            raise NotGroupMember.exception("You are not a member of this group")

    result = await session.exec(
        select(GroupMembership, User)
        .join(User, User.id == GroupMembership.user_id)
        .where(GroupMembership.group_id == group_id)
    )
    return CollectionResponse(
        items=[
            GroupMemberResponse(
                user_id=m.user_id,
                display_name=u.display_name,
                role=m.role,
                joined_at=m.joined_at,
            )
            for m, u in result.all()
        ]
    )


@router.post(
    "/{group_id}/members",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(GroupNotFound, NotGroupAdmin, UserNotFound),
)
async def add_member(
    session: SessionDep,
    current_user: CurrentUserDep,
    group_id: UUID,
    payload: GroupMemberCreateRequest,
) -> GroupMemberResponse:
    """Add a user to a group.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user (must be group admin or superuser).
    group_id
        UUID of the group to add a member to.
    payload
        Member creation request specifying user_id and optional role.
    """
    group = await session.get(Group, group_id)
    if group is None:
        raise GroupNotFound.exception("Group not found")
    if not current_user.is_superuser:
        await _require_admin(session, current_user.id, group_id)
    user = await session.get(User, payload.user_id)
    if user is None:
        raise UserNotFound.exception("User not found")
    existing = await fetch_group_role(session, payload.user_id, group_id)
    if existing is not None:
        # Idempotent: return existing membership
        result = await session.exec(
            select(GroupMembership).where(
                GroupMembership.group_id == group_id,
                GroupMembership.user_id == payload.user_id,
            )
        )
        m = result.one()
        return GroupMemberResponse(
            user_id=m.user_id,
            display_name=user.display_name,
            role=m.role,
            joined_at=m.joined_at,
        )
    membership = GroupMembership(
        group_id=group_id, user_id=payload.user_id, role=payload.role
    )
    session.add(membership)
    await session.commit()
    await session.refresh(membership)
    return GroupMemberResponse(
        user_id=membership.user_id,
        display_name=user.display_name,
        role=membership.role,
        joined_at=membership.joined_at,
    )


@router.patch(
    "/{group_id}/members/{user_id}",
    responses=problem_responses(
        GroupNotFound, NotGroupAdmin, LastGroupAdmin, UserNotFound
    ),
)
async def update_member_role(
    session: SessionDep,
    current_user: CurrentUserDep,
    group_id: UUID,
    user_id: UUID,
    payload: GroupMemberPatchRequest,
) -> GroupMemberResponse:
    """Update a group member's role.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user (must be group admin or superuser).
    group_id
        UUID of the group.
    user_id
        UUID of the member whose role to update.
    payload
        Role update request body.
    """
    if not current_user.is_superuser:
        await _require_admin(session, current_user.id, group_id)
    result = await session.exec(
        select(GroupMembership).where(
            GroupMembership.group_id == group_id,
            GroupMembership.user_id == user_id,
        )
    )
    membership = result.one_or_none()
    if membership is None:
        raise UserNotFound.exception("User is not a member of this group")

    if (
        membership.role == GroupRole.ADMIN
        and payload.role != GroupRole.ADMIN
        and await _count_admins(session, group_id) <= 1
    ):
        raise LastGroupAdmin.exception(
            "Cannot demote the last admin; promote another member first"
        )

    membership.role = payload.role
    await session.commit()

    user = await session.get(User, user_id)
    return GroupMemberResponse(
        user_id=membership.user_id,
        display_name=user.display_name if user else None,
        role=membership.role,
        joined_at=membership.joined_at,
    )


@router.delete(
    "/{group_id}/members/{user_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    responses=problem_responses(
        GroupNotFound, NotGroupAdmin, LastGroupAdmin, UserNotFound
    ),
)
async def remove_member(
    session: SessionDep,
    current_user: CurrentUserDep,
    group_id: UUID,
    user_id: UUID,
) -> None:
    """Remove a member from a group.

    Parameters
    ----------
    session
        Async database session.
    current_user
        Authenticated user (must be group admin, superuser, or self-removing).
    group_id
        UUID of the group.
    user_id
        UUID of the member to remove.
    """
    if current_user.id != user_id and not current_user.is_superuser:
        await _require_admin(session, current_user.id, group_id)
    result = await session.exec(
        select(GroupMembership).where(
            GroupMembership.group_id == group_id,
            GroupMembership.user_id == user_id,
        )
    )
    membership = result.one_or_none()
    if membership is None:
        raise UserNotFound.exception("User is not a member of this group")

    if (
        membership.role == GroupRole.ADMIN
        and await _count_admins(session, group_id) <= 1
    ):
        raise LastGroupAdmin.exception(
            "Cannot remove the last admin; promote another member first"
        )

    await session.delete(membership)
    await session.commit()
