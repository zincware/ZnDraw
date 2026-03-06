"""Custom FastAPI response classes for binary data formats.

Provides MessagePackResponse for efficient binary serialization of frame data
and other large datasets.
"""

from typing import Any

import msgpack
from fastapi import Response


class MessagePackResponse(Response):
    """FastAPI response class for MessagePack binary data.

    Automatically serializes Python objects to MessagePack format.
    If content is already bytes, it is returned as-is (useful for pre-packed data).

    Usage:
        @app.get("/data", response_class=MessagePackResponse)
        async def get_data():
            return {"key": "value", "numbers": [1, 2, 3]}

        # Or with pre-packed bytes:
        @app.get("/raw", response_class=MessagePackResponse)
        async def get_raw():
            return MessagePackResponse(content=raw_bytes)
    """

    media_type = "application/x-msgpack"

    def render(self, content: Any) -> bytes:
        """Render content to MessagePack bytes.

        Args:
            content: Python object to serialize, or raw bytes to return as-is.

        Returns:
            MessagePack-encoded bytes.
        """
        # If content is already bytes (pre-packed data), return as-is
        if isinstance(content, bytes):
            return content
        # Otherwise, serialize to MessagePack
        return msgpack.packb(content)
