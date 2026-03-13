"""Smoke tests for the exception hierarchy."""

import pytest

from zndraw.client.exceptions import NotConnectedError
from zndraw.exceptions import RoomLockedError, ZnDrawError


@pytest.mark.parametrize(
    ("cls", "parent"),
    [
        (ZnDrawError, Exception),
        (RoomLockedError, ZnDrawError),
        (NotConnectedError, ZnDrawError),
    ],
    ids=[
        "ZnDrawError<-Exception",
        "RoomLockedError<-ZnDrawError",
        "NotConnectedError<-ZnDrawError",
    ],
)
def test_inheritance(cls, parent):
    """Each exception is a subclass of its expected parent."""
    assert issubclass(cls, parent)


@pytest.mark.parametrize(
    "cls",
    [ZnDrawError, RoomLockedError, NotConnectedError],
    ids=["ZnDrawError", "RoomLockedError", "NotConnectedError"],
)
def test_instantiation_with_message(cls):
    """Each exception can be instantiated with a message string."""
    exc = cls("test message")
    assert exc.args[0] == "test message"


@pytest.mark.parametrize(
    "cls",
    [RoomLockedError, NotConnectedError],
    ids=["RoomLockedError", "NotConnectedError"],
)
def test_catch_as_zndraw_error(cls):
    """Subclass exceptions can be caught as ZnDrawError."""
    with pytest.raises(ZnDrawError):
        raise cls("test")
