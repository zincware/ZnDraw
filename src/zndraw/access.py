"""Access-control primitives for rooms, groups, and share links."""

from dataclasses import dataclass
from enum import StrEnum


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
