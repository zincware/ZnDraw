"""Built-in Provider subclasses shipped with the server.

These are registered at ``@internal`` scope on startup and dispatched
via the taskiq worker (see ``zndraw.providers.executor``).
"""

from zndraw.providers.filesystem import FilesystemRead

BUNDLED_PROVIDERS = [FilesystemRead]

__all__ = ["BUNDLED_PROVIDERS", "FilesystemRead"]
