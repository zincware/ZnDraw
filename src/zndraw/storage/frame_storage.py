"""Room-scoped AsyncBlobIO registry with provider frame count support."""

import base64
from typing import Any

import msgpack
from asebytes import AsyncBlobIO
from asebytes._async_adapters import AsyncObjectToBlobReadWriteAdapter
from asebytes._async_backends import AsyncReadWriteBackend
from redis.asyncio import Redis as AsyncRedis

from zndraw.redis import RedisKey

# Type alias for raw frame data: dict with bytes keys and bytes values
# Keys are like b"arrays.positions", b"cell", etc.
# Values are msgpack-numpy encoded bytes
RawFrame = dict[bytes, bytes]


def to_raw_frame(frame: dict[str, Any] | RawFrame) -> RawFrame:
    """Convert input frame dict to raw bytes format.

    Handles two input formats:
    1. Already raw: dict[bytes, bytes] -- pass through
    2. Base64 encoded: dict[str, str] with keys like "b64:..." -- decode

    Parameters
    ----------
    frame
        Input frame data (dict[str, Any] or dict[bytes, bytes]).
    """
    if not frame:
        return {}

    first_key = next(iter(frame.keys()))

    if isinstance(first_key, bytes):
        return frame  # type: ignore[return-value]

    result: RawFrame = {}
    for k, v in frame.items():
        key_bytes: bytes
        if isinstance(k, str) and k.startswith("b64:"):
            key_bytes = base64.b64decode(k[4:])
        elif isinstance(k, str):
            key_bytes = k.encode()
        elif isinstance(k, bytes):
            key_bytes = k
        else:
            key_bytes = str(k).encode()

        val_bytes: bytes
        if isinstance(v, str):
            val_bytes = base64.b64decode(v)
        elif isinstance(v, bytes):
            val_bytes = v
        else:
            packed = msgpack.packb(v)
            val_bytes = packed if packed is not None else b""

        result[key_bytes] = val_bytes

    return result


def _create_blob_backend(
    uri: str, group: str | None = None
) -> AsyncReadWriteBackend[bytes, bytes]:
    """Create an async blob backend from a URI string.

    All backends from asebytes registries are object-level (str keys).
    This function wraps them to blob-level (bytes keys) via
    AsyncObjectToBlobReadWriteAdapter, converting sync backends to
    async first when needed.

    Parameters
    ----------
    uri
        Storage backend URI string.
    group
        Optional group/namespace for data isolation (e.g., room_id).
    """
    from asebytes._registry import get_async_backend_cls

    cls = get_async_backend_cls(uri, readonly=False)
    if hasattr(cls, "from_uri"):
        backend = cls.from_uri(uri, group=group)
    else:
        backend = cls(uri, group=group)

    # Ensure async
    if not isinstance(backend, AsyncReadWriteBackend):
        from asebytes._async_backends import sync_to_async

        backend = sync_to_async(backend)  # type: ignore[arg-type]

    # All registry backends are object-level (str keys) → wrap to blob
    return AsyncObjectToBlobReadWriteAdapter(backend)  # type: ignore[arg-type]


class FrameStorage:
    """Room-scoped AsyncBlobIO registry with provider frame count support.

    Provides lazy-created ``AsyncBlobIO`` instances per room via subscript
    access. Routes use the ``AsyncBlobIO`` pandas-like API directly rather
    than going through wrapper methods.

    Parameters
    ----------
    uri
        Storage backend URI string (e.g. ``memory://``, ``path/to/data.lmdb``).
    redis
        Async Redis client for provider frame count metadata.
    """

    def __init__(self, uri: str, redis: AsyncRedis) -> None:  # type: ignore[type-arg]
        self._uri = uri
        self._redis = redis
        self._rooms: dict[str, AsyncBlobIO] = {}

    def __getitem__(self, room_id: str) -> AsyncBlobIO:
        """Get or lazily create an AsyncBlobIO for a room.

        Parameters
        ----------
        room_id
            Room identifier.
        """
        if room_id not in self._rooms:
            backend = _create_blob_backend(self._uri, group=room_id)
            self._rooms[room_id] = AsyncBlobIO(backend)
        return self._rooms[room_id]

    async def get_length(self, room_id: str) -> int:
        """Get frame count for a room.

        Checks the AsyncBlobIO length first, then falls back to the
        Redis provider frame count key (set by mounted providers).

        Parameters
        ----------
        room_id
            Room identifier.
        """
        io = self[room_id]
        length = await io.len()
        if length > 0:
            return length
        cached = await self._redis.get(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id)
        )
        return int(cached) if cached else 0

    async def has_mount(self, room_id: str) -> bool:
        """Check if a room has a provider mount (read-only).

        Parameters
        ----------
        room_id
            Room identifier.
        """
        return await self._redis.exists(RedisKey.provider_frame_count(room_id)) > 0  # type: ignore[misc]

    async def set_frame_count(self, room_id: str, count: int) -> None:
        """Store provider frame count in Redis.

        Parameters
        ----------
        room_id
            Room identifier.
        count
            Number of frames the provider exposes.
        """
        await self._redis.set(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id), count
        )

    async def clear_frame_count(self, room_id: str) -> None:
        """Remove provider frame count from Redis.

        Parameters
        ----------
        room_id
            Room identifier.
        """
        await self._redis.delete(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id)
        )

    async def clear(self, room_id: str) -> None:
        """Clear all frame data and provider frame count for a room.

        Parameters
        ----------
        room_id
            Room identifier.
        """
        io = self[room_id]
        await io.clear()
        await self.clear_frame_count(room_id)

    async def close(self) -> None:
        """Release in-memory handles without deleting stored data."""
        self._rooms.clear()
