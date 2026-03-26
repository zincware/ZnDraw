"""Provider base model for client-side read handlers."""

from __future__ import annotations

from typing import Any, ClassVar

from pydantic import BaseModel


class Provider(BaseModel):
    """Base model for provider read handlers.

    Host apps subclass this to define typed read requests for a specific
    provider kind (e.g., filesystem, frame_source). The ``category`` ClassVar
    groups providers so the router can dispatch to the correct handler,
    and ``read()`` performs the actual work against a caller-supplied
    handler object.

    Example::

        class FilesystemRead(Provider):
            category: ClassVar[str] = "filesystem"
            path: str

            def read(self, handler: Any) -> Any:
                return handler.cat(self.path)
    """

    category: ClassVar[str]
    content_type: ClassVar[str] = "application/json"

    def read(self, handler: Any) -> Any:
        """Process a read request against *handler*.

        Subclasses must override this method. *handler* is the
        provider-specific backend object (e.g., an fsspec filesystem).
        The return value must be JSON-serializable (for ``application/json``
        providers) or ``bytes`` (for binary content types).
        """
        raise NotImplementedError
