"""StorageRouter: routes storage calls by room_id with provider support."""

from typing import Any

from redis.asyncio import Redis as AsyncRedis

from zndraw.exceptions import RoomReadOnly
from zndraw.redis import RedisKey
from zndraw.storage.base import RawFrame, StorageBackend


class StorageRouter(StorageBackend):
    """Wraps a default StorageBackend with provider-backed room support.

    Provider-backed rooms (via ``set_frame_count``) store their frame count
    in Redis.  Write methods reject modifications to these rooms.
    """

    def __init__(self, default: StorageBackend, redis: AsyncRedis) -> None:  # type: ignore[type-arg]
        self._default = default
        self._redis = redis

    async def has_mount(self, room_id: str) -> bool:
        """Check if a room has provider-backed frames (read-only)."""
        return await self._redis.exists(RedisKey.provider_frame_count(room_id)) > 0  # type: ignore[misc]

    # -- StorageBackend interface: delegate to default backend -----------------

    async def get(self, room_id: str, index: int) -> RawFrame | None:
        return await self._default.get(room_id, index)

    async def get_range(
        self, room_id: str, start: int, stop: int | None
    ) -> list[RawFrame | None]:
        total = await self.get_length(room_id)
        actual_stop = min(stop, total) if stop is not None else total
        result = await self._default.get_range(room_id, start, actual_stop)
        expected = max(0, actual_stop - start)
        if len(result) < expected and await self.has_mount(room_id):
            result.extend([None] * (expected - len(result)))
        return result

    async def get_many(self, room_id: str, indices: list[int]) -> list[RawFrame | None]:
        return await self._default.get_many(room_id, indices)

    async def get_length(self, room_id: str) -> int:
        """Get frame count: default storage first, then Redis fallback."""
        length = await self._default.get_length(room_id)
        if length > 0:
            return length
        # Fallback: check Redis for provider-set frame count
        cached = await self._redis.get(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id)
        )
        return int(cached) if cached else 0

    # -- Write methods: reject for provider-backed rooms ----------------------

    async def extend(
        self, room_id: str, frames: list[dict[str, Any]] | list[RawFrame]
    ) -> int:
        if await self.has_mount(room_id):
            raise RoomReadOnly.exception("Room is provider-backed (read-only)")
        return await self._default.extend(room_id, frames)

    async def set_item(
        self, room_id: str, index: int, frame: dict[str, Any] | RawFrame
    ) -> None:
        if await self.has_mount(room_id):
            raise RoomReadOnly.exception("Room is provider-backed (read-only)")
        await self._default.set_item(room_id, index, frame)

    async def merge_item(self, room_id: str, index: int, partial: RawFrame) -> None:
        if await self.has_mount(room_id):
            raise RoomReadOnly.exception("Room is provider-backed (read-only)")
        await self._default.merge_item(room_id, index, partial)

    async def delete_range(self, room_id: str, start: int, stop: int) -> None:
        if await self.has_mount(room_id):
            raise RoomReadOnly.exception("Room is provider-backed (read-only)")
        await self._default.delete_range(room_id, start, stop)

    # -- Provider frame count management --------------------------------------

    async def set_frame_count(self, room_id: str, count: int) -> None:
        """Store provider frame count in Redis."""
        await self._redis.set(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id), count
        )

    async def clear_frame_count(self, room_id: str) -> None:
        """Remove provider frame count from Redis."""
        await self._redis.delete(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id)
        )

    async def clear(self, room_id: str) -> None:
        await self._default.clear(room_id)
        await self.clear_frame_count(room_id)

    async def reserve(self, room_id: str, count: int) -> None:
        await self._default.reserve(room_id, count)

    async def remove_items(self, room_id: str, indices: list[int]) -> None:
        await self._default.remove_items(room_id, indices)

    async def close(self) -> None:
        await self._default.close()
