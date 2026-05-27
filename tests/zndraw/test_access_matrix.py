"""Pure-predicate tests for can_read / can_edit / can_manage."""

from uuid import UUID, uuid4

from zndraw.access import (
    GroupRole,
    ShareAccess,
    ShareContext,
    Visibility,
    can_edit,
    can_manage,
    can_read,
)


class _Room:
    def __init__(
        self,
        *,
        visibility: Visibility,
        owner_user_id: UUID | None = None,
        owner_group_id: UUID | None = None,
    ) -> None:
        self.visibility = visibility
        self.owner_user_id = owner_user_id
        self.owner_group_id = owner_group_id


class _User:
    def __init__(self, uid: UUID, *, is_superuser: bool = False) -> None:
        self.id = uid
        self.is_superuser = is_superuser


def test_public_room_readable_by_anyone() -> None:
    room = _Room(visibility=Visibility.PUBLIC, owner_user_id=uuid4())
    stranger = _User(uuid4())
    assert can_read(stranger, room, None, group_role=None) is True


def test_private_room_hidden_from_non_owner() -> None:
    owner_id = uuid4()
    room = _Room(visibility=Visibility.PRIVATE, owner_user_id=owner_id)
    stranger = _User(uuid4())
    assert can_read(stranger, room, None, group_role=None) is False


def test_private_room_visible_to_owner() -> None:
    owner_id = uuid4()
    room = _Room(visibility=Visibility.PRIVATE, owner_user_id=owner_id)
    owner = _User(owner_id)
    assert can_read(owner, room, None, group_role=None) is True


def test_group_room_visible_to_member() -> None:
    gid = uuid4()
    room = _Room(visibility=Visibility.GROUP, owner_group_id=gid)
    member = _User(uuid4())
    assert can_read(member, room, None, group_role=GroupRole.VIEWER) is True


def test_group_room_hidden_from_non_member() -> None:
    gid = uuid4()
    room = _Room(visibility=Visibility.GROUP, owner_group_id=gid)
    stranger = _User(uuid4())
    assert can_read(stranger, room, None, group_role=None) is False


def test_superuser_reads_everything() -> None:
    room = _Room(visibility=Visibility.PRIVATE, owner_user_id=uuid4())
    admin = _User(uuid4(), is_superuser=True)
    assert can_read(admin, room, None, group_role=None) is True


def test_viewer_role_cannot_edit_group_room() -> None:
    gid = uuid4()
    room = _Room(visibility=Visibility.GROUP, owner_group_id=gid)
    viewer = _User(uuid4())
    assert can_edit(viewer, room, None, group_role=GroupRole.VIEWER) is False


def test_member_role_can_edit_group_room() -> None:
    gid = uuid4()
    room = _Room(visibility=Visibility.GROUP, owner_group_id=gid)
    member = _User(uuid4())
    assert can_edit(member, room, None, group_role=GroupRole.MEMBER) is True


def test_admin_role_can_manage_group_room() -> None:
    gid = uuid4()
    room = _Room(visibility=Visibility.GROUP, owner_group_id=gid)
    admin = _User(uuid4())
    assert can_manage(admin, room, group_role=GroupRole.ADMIN) is True


def test_member_role_cannot_manage_group_room() -> None:
    gid = uuid4()
    room = _Room(visibility=Visibility.GROUP, owner_group_id=gid)
    member = _User(uuid4())
    assert can_manage(member, room, group_role=GroupRole.MEMBER) is False


def test_public_user_owned_room_is_chaotic_edit() -> None:
    owner_id = uuid4()
    room = _Room(visibility=Visibility.PUBLIC, owner_user_id=owner_id)
    stranger = _User(uuid4())
    assert can_edit(stranger, room, None, group_role=None) is True


def test_public_group_owned_room_is_read_only_for_non_members() -> None:
    gid = uuid4()
    room = _Room(visibility=Visibility.PUBLIC, owner_group_id=gid)
    stranger = _User(uuid4())
    assert can_read(stranger, room, None, group_role=None) is True
    assert can_edit(stranger, room, None, group_role=None) is False


def test_share_view_grants_read_only() -> None:
    room = _Room(visibility=Visibility.PRIVATE, owner_user_id=uuid4())
    stranger = _User(uuid4())
    share = ShareContext(room_id="r", access=ShareAccess.VIEW)
    assert can_read(stranger, room, share, group_role=None) is True
    assert can_edit(stranger, room, share, group_role=None) is False


def test_share_edit_grants_edit() -> None:
    room = _Room(visibility=Visibility.PRIVATE, owner_user_id=uuid4())
    stranger = _User(uuid4())
    share = ShareContext(room_id="r", access=ShareAccess.EDIT)
    assert can_edit(stranger, room, share, group_role=None) is True


def test_share_never_grants_manage() -> None:
    room = _Room(visibility=Visibility.PRIVATE, owner_user_id=uuid4())
    stranger = _User(uuid4())
    assert can_manage(stranger, room, group_role=None) is False
