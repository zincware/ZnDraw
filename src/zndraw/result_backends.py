"""Result backends for the zndraw-joblib provider system.

- ``RedisResultBackend`` — general-purpose, for small JSON payloads.
- ``StorageResultBackend`` — adapts any ``StorageBackend`` to the
  ``ResultBackend`` protocol.  Reuses the same storage infrastructure
  as the main frame storage (memory, LMDB, MongoDB, …).
- ``CompositeResultBackend`` — routes by cache-key pattern: frame
  provider results go to the storage backend, everything else to Redis.
  Inflight locks always use Redis (atomic SET NX semantics).
"""

from __future__ import annotations

import asyncio
from collections.abc import Awaitable, Callable

from redis.asyncio import Redis
from redis.asyncio.client import PubSub

from zndraw.storage.base import StorageBackend

NOTIFY_PREFIX = "notify:"


async def _pubsub_wait(
    pubsub_factory: Callable[[], PubSub],
    get_fn: Callable[[str], Awaitable[bytes | None]],
    key: str,
    timeout: float,
) -> bytes | None:
    """Wait for a cache key via Redis pub/sub.

    Subscribe first, then check cache (race-safe).
    Loop handles subscribe confirmation messages that
    ``ignore_subscribe_messages=True`` filters as None.
    """
    channel = f"{NOTIFY_PREFIX}{key}"
    loop = asyncio.get_running_loop()
    async with pubsub_factory() as pubsub:
        await pubsub.subscribe(channel)

        # Check cache AFTER subscribing — handles the race where
        # the result landed between our first check and the subscribe.
        cached = await get_fn(key)
        if cached is not None:
            return cached

        deadline = loop.time() + timeout
        while (remaining := deadline - loop.time()) > 0:
            msg = await pubsub.get_message(
                ignore_subscribe_messages=True, timeout=remaining
            )
            if msg is not None:
                return await get_fn(key)

        return None


class RedisResultBackend:
    """Store provider results in Redis with TTL."""

    def __init__(self, redis: Redis) -> None:
        self._redis = redis

    async def store(self, key: str, data: bytes, ttl: int) -> None:
        await self._redis.set(key, data, ex=ttl)

    async def get(self, key: str) -> bytes | None:
        return await self._redis.get(key)

    async def delete(self, key: str) -> None:
        await self._redis.delete(key)

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        return bool(await self._redis.set(key, b"1", nx=True, ex=ttl))

    async def release_inflight(self, key: str) -> None:
        await self._redis.delete(key)

    def pubsub(self) -> PubSub:
        """Return a pub/sub instance for this Redis connection."""
        return self._redis.pubsub()

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:
        """Wait for a cache key via Redis pub/sub (race-safe)."""
        return await _pubsub_wait(self.pubsub, self.get, key, timeout)

    async def notify_key(self, key: str) -> None:
        """Publish notification that a cache key has been populated."""
        await self._redis.publish(f"{NOTIFY_PREFIX}{key}", b"1")


class StorageResultBackend:
    """Adapt any ``StorageBackend`` to the ``ResultBackend`` protocol.

    Uses cache keys as room_ids and stores raw bytes as a single-entry
    frame (``{b"_": data}``).  This means the provider cache
    automatically uses whatever storage the server is configured with
    (in-memory, LMDB, MongoDB, …) — no storage-specific code needed.

    TTL is not enforced at this layer; entries are cleaned up when
    providers disconnect.

    Parameters
    ----------
    storage
        The storage backend to delegate to (same instance as the main
        frame storage is fine — cache keys are namespaced).
    """

    def __init__(self, storage: StorageBackend) -> None:
        self._storage = storage

    async def store(self, key: str, data: bytes, ttl: int) -> None:
        await self._storage.clear(key)
        await self._storage.extend(key, [{b"_": data}])

    async def get(self, key: str) -> bytes | None:
        frame = await self._storage.get(key, 0)
        if frame is None:
            return None
        return frame.get(b"_")

    async def delete(self, key: str) -> None:
        await self._storage.clear(key)

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        raise NotImplementedError("Use Redis for inflight locks")

    async def release_inflight(self, key: str) -> None:
        raise NotImplementedError("Use Redis for inflight locks")

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:
        raise NotImplementedError("Use CompositeResultBackend for long-polling")

    async def notify_key(self, key: str) -> None:
        raise NotImplementedError("Use CompositeResultBackend for long-polling")


class CompositeResultBackend:
    """Route provider results to storage (frames) or Redis (everything else).

    Cache keys follow ``provider-result:{room}:{category}:{name}:{hash}``.
    Keys containing ``:frames:`` are routed to the storage backend; all
    others go to Redis.  Inflight locks always use Redis for atomic SET NX.
    Pub/sub notifications always go through Redis.

    Parameters
    ----------
    redis
        Backend for JSON/small payloads and all inflight locks.
    frames
        Backend for large binary (frame) payloads.
    """

    def __init__(
        self,
        redis: RedisResultBackend,
        frames: StorageResultBackend,
    ) -> None:
        self._redis = redis
        self._frames = frames

    def _backend_for(self, key: str) -> RedisResultBackend | StorageResultBackend:
        """Select the backend based on cache-key pattern."""
        return self._frames if ":frames:" in key else self._redis

    async def store(self, key: str, data: bytes, ttl: int) -> None:
        await self._backend_for(key).store(key, data, ttl)

    async def get(self, key: str) -> bytes | None:
        return await self._backend_for(key).get(key)

    async def delete(self, key: str) -> None:
        await self._backend_for(key).delete(key)

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        return await self._redis.acquire_inflight(key, ttl)

    async def release_inflight(self, key: str) -> None:
        await self._redis.release_inflight(key)

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:
        """Wait via Redis pub/sub, read from the correct backend."""
        return await _pubsub_wait(self._redis.pubsub, self.get, key, timeout)

    async def notify_key(self, key: str) -> None:
        await self._redis.notify_key(key)
