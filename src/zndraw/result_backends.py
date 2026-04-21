"""Result backends for the zndraw-joblib provider system.

- ``RedisResultBackend`` — general-purpose, for small JSON payloads.
- ``StorageResultBackend`` — adapts ``FrameStorage`` to the
  ``ResultBackend`` protocol.  Reuses the same storage infrastructure
  as the main frame storage (memory, LMDB, MongoDB, …).
- ``CompositeResultBackend`` — routes by cache-key pattern: frame
  provider results go to the storage backend, everything else to Redis.
  Inflight locks always use Redis (atomic SET NX semantics).
"""

from __future__ import annotations

import asyncio
from typing import TYPE_CHECKING

import msgpack

if TYPE_CHECKING:
    from collections.abc import Awaitable, Callable

    from redis.asyncio import Redis
    from redis.asyncio.client import PubSub

    from zndraw.storage import FrameStorage

NOTIFY_PREFIX = "notify:"


async def _pubsub_wait_prefixed(
    pubsub_factory: Callable[[], PubSub],
    get_fn: Callable[[str], Awaitable[bytes | None]],
    bare_key: str,
    channel_key: str,
    timeout: float,  # noqa: ASYNC109
) -> bytes | None:
    """Wait for a cache key via pub/sub, subscribing on the prefixed channel
    but reading the bare (un-prefixed) key through ``get_fn``.

    Subscribe first, then check cache (race-safe).
    Loop handles subscribe confirmation messages that
    ``ignore_subscribe_messages=True`` filters as None.

    Parameters
    ----------
    pubsub_factory
        Callable returning a Redis PubSub instance (used as context manager).
    get_fn
        Callable that reads the cache by bare key.
    bare_key
        The key passed to ``get_fn`` (pre-transformed by the backend's ``_k``
        method, i.e. the key the backend's public API understands).
    channel_key
        The fully-prefixed key used to construct the pub/sub channel.
    timeout
        Maximum seconds to wait.
    """
    channel = f"{NOTIFY_PREFIX}{channel_key}"
    loop = asyncio.get_running_loop()
    async with pubsub_factory() as pubsub:
        await pubsub.subscribe(channel)

        # Check cache AFTER subscribing — handles the race where
        # the result landed between our first check and the subscribe.
        cached = await get_fn(bare_key)
        if cached is not None:
            return cached

        deadline = loop.time() + timeout
        while (remaining := deadline - loop.time()) > 0:
            msg = await pubsub.get_message(
                ignore_subscribe_messages=True, timeout=remaining
            )
            if msg is not None:
                return await get_fn(bare_key)

        return None


class RedisResultBackend:
    """Store provider results in Redis with TTL.

    Parameters
    ----------
    redis
        Async Redis client instance.
    key_prefix
        Optional namespace prefix.  When non-empty every Redis key and the
        pub/sub notification channel are stored as ``{key_prefix}:{key}``.
        Use this to isolate multiple server instances sharing one Redis.
    """

    def __init__(self, redis: Redis, key_prefix: str = "") -> None:
        self._redis = redis
        self._key_prefix = key_prefix

    def _namespaced_key(self, key: str) -> str:
        """Return the namespaced key, prepending the prefix when set."""
        if self._key_prefix:
            return f"{self._key_prefix}:{key}"
        return key

    async def store(self, key: str, data: bytes, ttl: int) -> None:
        await self._redis.set(self._namespaced_key(key), data, ex=ttl)

    async def get(self, key: str) -> bytes | None:
        return await self._redis.get(self._namespaced_key(key))

    async def delete(self, key: str) -> None:
        await self._redis.delete(self._namespaced_key(key))

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        return bool(
            await self._redis.set(self._namespaced_key(key), b"1", nx=True, ex=ttl)
        )

    async def release_inflight(self, key: str) -> None:
        await self._redis.delete(self._namespaced_key(key))

    def pubsub(self) -> PubSub:
        """Return a pub/sub instance for this Redis connection."""
        return self._redis.pubsub()

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:  # noqa: ASYNC109
        """Wait for a cache key via Redis pub/sub (race-safe)."""
        return await _pubsub_wait_prefixed(
            self.pubsub, self.get, key, self._namespaced_key(key), timeout
        )

    async def notify_key(self, key: str) -> None:
        """Publish notification that a cache key has been populated."""
        await self._redis.publish(f"{NOTIFY_PREFIX}{self._namespaced_key(key)}", b"1")


class StorageResultBackend:
    """Adapt ``FrameStorage`` to the ``ResultBackend`` protocol.

    Uses cache keys as room_ids and stores raw bytes as a single-entry
    frame (``{b"_": data}``).  This means the provider cache
    automatically uses whatever storage the server is configured with
    (memory, LMDB, MongoDB, …) — no storage-specific code needed.

    TTL is not enforced at this layer; entries are cleaned up when
    providers disconnect.

    Parameters
    ----------
    storage
        The FrameStorage registry to delegate to (same instance as the
        main frame storage is fine — cache keys are namespaced via
        ``key_prefix``).
    key_prefix
        Optional namespace prefix.  When non-empty every storage key is
        stored as ``{key_prefix}:{key}``. Required for multi-server
        deployments sharing a single FrameStorage instance.
    """

    def __init__(self, storage: FrameStorage, key_prefix: str = "") -> None:
        self._storage = storage
        self._key_prefix = key_prefix

    def _namespaced_key(self, key: str) -> str:
        """Return the namespaced key, prepending the prefix when set."""
        if self._key_prefix:
            return f"{self._key_prefix}:{key}"
        return key

    async def store(self, key: str, data: bytes, ttl: int) -> None:  # noqa: ARG002
        io = self._storage[self._namespaced_key(key)]
        await io.clear()
        # Wrap in msgpack so the blob↔object adapter round-trip works
        packed = msgpack.packb(data)
        assert packed is not None
        await io.extend([{b"_": packed}])

    async def get(self, key: str) -> bytes | None:
        io = self._storage[self._namespaced_key(key)]
        if await io.len() == 0:
            return None
        packed = await io[b"_"][0]
        if packed is None:
            return None
        return msgpack.unpackb(packed)

    async def delete(self, key: str) -> None:
        await self._storage[self._namespaced_key(key)].clear()

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        raise NotImplementedError("Use Redis for inflight locks")

    async def release_inflight(self, key: str) -> None:
        raise NotImplementedError("Use Redis for inflight locks")

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:  # noqa: ASYNC109
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

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:  # noqa: ASYNC109
        """Wait via Redis pub/sub, read from the correct backend."""
        return await _pubsub_wait_prefixed(
            self._redis.pubsub,
            self.get,
            key,
            self._redis._namespaced_key(key),  # noqa: SLF001 — same-module helper
            timeout,
        )

    async def notify_key(self, key: str) -> None:
        await self._redis.notify_key(key)
