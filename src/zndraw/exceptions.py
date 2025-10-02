"""ZnDraw exception classes."""


class ZnDrawException(Exception):
    """Base exception for all ZnDraw errors."""

    pass


class LockError(ZnDrawException):
    """Raised when an operation fails due to a lock being held by another client."""

    pass
