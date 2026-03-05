"""Client-specific exception classes."""

from zndraw.exceptions import ZnDrawError


class NotConnectedError(ZnDrawError):
    """Raised when operation requires connection but client is disconnected."""
