class ZnDrawException(Exception):
    """Base ZnDraw exception."""


class RoomNotFound(Exception):
    """Raised when a room is not found on the server."""


class RoomLockedError(ZnDrawException):
    """Raised when tried to modify a locked room."""
