"""TDD tests for new ZnDraw client API surface.

Tests for classmethods (list_rooms, login),
instance properties (locked, chat, screenshots, extensions, tasks).
"""

import uuid
import warnings

from zndraw import ZnDraw
from zndraw.extensions import Category, Extension
from zndraw.schemas import MessageResponse, ScreenshotListItem

# =============================================================================
# ZnDraw.list_rooms classmethod
# =============================================================================


def test_list_rooms_returns_list(server: str):
    """list_rooms returns a list (possibly empty)."""
    rooms = ZnDraw.list_rooms(url=server)
    assert isinstance(rooms, list)


def test_list_rooms_contains_created_room(server: str):
    """A room created via ZnDraw appears in list_rooms."""
    vis = ZnDraw(url=server)
    room_id = vis.room
    vis.disconnect()

    rooms = ZnDraw.list_rooms(url=server)
    room_ids = [r["id"] for r in rooms]
    assert room_id in room_ids


def test_list_rooms_search_filters(server: str):
    """list_rooms with search only returns matching rooms."""
    room_name = f"searchable-{uuid.uuid4().hex[:8]}"
    vis = ZnDraw(url=server, room=room_name)
    vis.disconnect()

    rooms = ZnDraw.list_rooms(url=server, search="searchable")
    room_ids = [r["id"] for r in rooms]
    assert room_name in room_ids


# =============================================================================
# ZnDraw.login classmethod
# =============================================================================


def test_login_returns_token(server_auth: str):
    """login with valid credentials returns a token string."""
    token = ZnDraw.login(
        url=server_auth,
        username="admin@local.test",
        password="adminpassword",
    )
    assert isinstance(token, str)
    assert len(token) > 0


# =============================================================================
# ZnDraw.locked property
# =============================================================================


def test_locked_default_false(server: str):
    """New rooms are unlocked by default."""
    vis = ZnDraw(url=server)
    assert vis.locked is False
    vis.disconnect()


def test_locked_roundtrip(server: str):
    """Setting locked=True locks the room, False unlocks."""
    vis = ZnDraw(url=server)
    vis.locked = True
    assert vis.locked is True

    vis.locked = False
    assert vis.locked is False
    vis.disconnect()


# =============================================================================
# ZnDraw.chat Sequence property
# =============================================================================


def test_chat_empty_room(server: str):
    """Chat on a new room is empty."""
    vis = ZnDraw(url=server)
    assert len(vis.chat) == 0
    vis.disconnect()


def test_chat_send_and_read(server: str):
    """Sending a message via chat.send makes it appear as MessageResponse."""
    vis = ZnDraw(url=server)
    vis.chat.send("hello world")
    assert len(vis.chat) >= 1
    msg = vis.chat[-1]
    assert isinstance(msg, MessageResponse)
    assert msg.content == "hello world"
    vis.disconnect()


def test_chat_log_delegates_to_chat_send(server: str):
    """vis.log() is deprecated and delegates to vis.chat.send()."""
    vis = ZnDraw(url=server)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        vis.log("via log")
    assert len(vis.chat) >= 1
    assert vis.chat[-1].content == "via log"
    vis.disconnect()


def test_chat_log_emits_deprecation_warning(server: str):
    """vis.log() emits a DeprecationWarning."""
    vis = ZnDraw(url=server)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        vis.log("test")
        assert any(issubclass(x.category, DeprecationWarning) for x in w)
    vis.disconnect()


def test_chat_slicing(server: str):
    """Chat supports slicing for pagination."""
    vis = ZnDraw(url=server)
    for i in range(5):
        vis.chat.send(f"msg-{i}")
    last_two = vis.chat[-2:]
    assert len(last_two) == 2
    assert last_two[-1].content == "msg-4"
    vis.disconnect()


def test_chat_iteration(server: str):
    """Chat is iterable."""
    vis = ZnDraw(url=server)
    vis.chat.send("first")
    vis.chat.send("second")
    messages = list(vis.chat)
    assert len(messages) == 2
    assert all(isinstance(m, MessageResponse) for m in messages)
    vis.disconnect()


# =============================================================================
# ZnDraw.screenshots Mapping property
# =============================================================================


def test_screenshots_empty_room(server: str):
    """New room has no screenshots."""
    vis = ZnDraw(url=server)
    assert len(vis.screenshots) == 0
    assert list(vis.screenshots) == []
    vis.disconnect()


# =============================================================================
# ZnDraw.extensions Mapping property
# =============================================================================


def test_extensions_is_mapping(server: str):
    """vis.extensions is iterable and has len."""
    vis = ZnDraw(url=server)
    ext_names = list(vis.extensions)
    assert len(ext_names) > 0
    assert len(vis.extensions) == len(ext_names)
    vis.disconnect()


def test_extensions_contains_internal(server: str):
    """Internal extensions like Delete are present."""
    vis = ZnDraw(url=server)
    matches = [k for k in vis.extensions if "Delete" in k]
    assert len(matches) == 1
    vis.disconnect()


def test_extensions_getitem_returns_extension_subclass(server: str):
    """vis.extensions[name] returns an Extension subclass (not instance)."""
    vis = ZnDraw(url=server)
    key = next(k for k in vis.extensions if "Delete" in k)
    cls = vis.extensions[key]
    assert isinstance(cls, type)
    assert issubclass(cls, Extension)
    vis.disconnect()


def test_extensions_class_has_category(server: str):
    """Returned extension class has the correct category."""
    vis = ZnDraw(url=server)
    key = next(k for k in vis.extensions if "Delete" in k)
    cls = vis.extensions[key]
    assert cls.category == Category.MODIFIER
    vis.disconnect()


def test_extensions_class_has_correct_name(server: str):
    """Returned extension class __name__ matches the server name."""
    vis = ZnDraw(url=server)
    key = next(k for k in vis.extensions if "Delete" in k)
    cls = vis.extensions[key]
    assert cls.__name__ == "Delete"
    vis.disconnect()


def test_extensions_class_is_instantiable(server: str):
    """Returned extension class can be instantiated."""
    vis = ZnDraw(url=server)
    key = next(k for k in vis.extensions if "Delete" in k)
    cls = vis.extensions[key]
    instance = cls()
    assert isinstance(instance, Extension)
    vis.disconnect()


def test_extensions_keyerror_on_missing(server: str):
    """Accessing a non-existent extension raises KeyError."""
    vis = ZnDraw(url=server)
    import pytest

    with pytest.raises(KeyError):
        vis.extensions["nonexistent:fake:NoSuchExtension"]
    vis.disconnect()


# =============================================================================
# ZnDraw.tasks Mapping property
# =============================================================================


def test_tasks_empty_room(server: str):
    """New room has no tasks."""
    vis = ZnDraw(url=server)
    assert len(vis.tasks) == 0
    assert list(vis.tasks) == []
    vis.disconnect()


def test_tasks_is_mapping(server: str):
    """vis.tasks is iterable and has len."""
    vis = ZnDraw(url=server)
    task_ids = list(vis.tasks)
    assert isinstance(task_ids, list)
    assert len(vis.tasks) == len(task_ids)
    vis.disconnect()


def test_tasks_callable_filter(server: str):
    """vis.tasks(status='pending') returns a filtered view."""
    vis = ZnDraw(url=server)
    filtered = vis.tasks(status="pending")
    assert len(filtered) >= 0
    vis.disconnect()
