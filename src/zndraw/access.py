"""Access-control primitives for rooms, groups, and share links."""

from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol
from uuid import UUID


class Visibility(StrEnum):
    PRIVATE = "private"
    GROUP = "group"
    PUBLIC = "public"


class GroupRole(StrEnum):
    VIEWER = "viewer"
    MEMBER = "member"
    ADMIN = "admin"


class ShareAccess(StrEnum):
    VIEW = "view"
    EDIT = "edit"


@dataclass(frozen=True, slots=True)
class ShareContext:
    """A validated share-token context attached to a request.

    The resolver guarantees instances correspond to a token that is
    non-revoked, non-expired, and scoped to ``room_id``. Permission
    checks treat ``None`` as "no share".
    """

    room_id: str
    access: ShareAccess


class _RoomLike(Protocol):
    visibility: Visibility
    owner_user_id: UUID | None
    owner_group_id: UUID | None


class _UserLike(Protocol):
    id: UUID
    is_superuser: bool


def can_read(
    user: _UserLike,
    room: _RoomLike,
    share: ShareContext | None,
    *,
    group_role: GroupRole | None,
) -> bool:
    """Return True if ``user`` may read ``room``.

    Parameters
    ----------
    user
        The requesting user (must be authenticated — anonymous access is
        removed in this refactor).
    room
        The target room.
    share
        Resolver-validated share context, or ``None``. When not ``None``,
        it is guaranteed to match ``room.id``.
    group_role
        Caller's role in ``room.owner_group_id``, or ``None`` when the
        room is user-owned or the caller is not a member of the group.
    """
    if user.is_superuser:
        return True
    if room.visibility == Visibility.PUBLIC:
        return True
    if room.owner_user_id is not None and room.owner_user_id == user.id:
        return True
    if room.owner_group_id is not None and group_role is not None:
        return True
    return share is not None


def can_edit(
    user: _UserLike,
    room: _RoomLike,
    share: ShareContext | None,
    *,
    group_role: GroupRole | None,
) -> bool:
    """Return True if ``user`` may edit ``room`` content.

    Parameters
    ----------
    user
        The requesting user.
    room
        The target room.
    share
        Resolver-validated share context, or ``None``.
    group_role
        Caller's role in ``room.owner_group_id``, or ``None``.
    """
    if user.is_superuser:
        return True
    if room.owner_user_id is not None and room.owner_user_id == user.id:
        return True
    if room.owner_group_id is not None and group_role in (
        GroupRole.MEMBER,
        GroupRole.ADMIN,
    ):
        return True
    if room.owner_user_id is not None and room.visibility == Visibility.PUBLIC:
        return True  # chaotic-edit default for user-owned public rooms
    return share is not None and share.access == ShareAccess.EDIT


def can_manage(
    user: _UserLike,
    room: _RoomLike,
    *,
    group_role: GroupRole | None,
) -> bool:
    """Return True if ``user`` may delete, transfer, or change visibility.

    Share tokens never grant manage rights.

    Parameters
    ----------
    user
        The requesting user.
    room
        The target room.
    group_role
        Caller's role in ``room.owner_group_id``, or ``None``.
    """
    if user.is_superuser:
        return True
    if room.owner_user_id is not None and room.owner_user_id == user.id:
        return True
    return room.owner_group_id is not None and group_role == GroupRole.ADMIN
