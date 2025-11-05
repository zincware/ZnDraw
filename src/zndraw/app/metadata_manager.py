"""Room metadata manager for Redis operations."""

import logging

from .redis_keys import RoomKeys

log = logging.getLogger(__name__)


class RoomMetadataManager:
    """Manages room metadata in Redis.

    Provides CRUD operations for room metadata with proper validation
    and atomic operations. All metadata values are stored as strings
    in a Redis hash at `room:{room_id}:metadata`.

    Parameters
    ----------
    redis_client : Redis
        Redis client instance
    room_id : str
        Room identifier

    Examples
    --------
    >>> manager = RoomMetadataManager(redis_client, "my_room")
    >>> manager.set("file_path", "data/file.xyz")
    >>> manager.get("file_path")
    'data/file.xyz'
    >>> manager.get_all()
    {'file_path': 'data/file.xyz'}
    """

    def __init__(self, redis_client, room_id: str):
        self.redis = redis_client
        self.room_id = room_id
        room_keys = RoomKeys(room_id)
        self.key = room_keys.metadata()

    def get_all(self) -> dict[str, str]:
        """Get all metadata for the room.

        Returns
        -------
        dict[str, str]
            Dictionary of all metadata key-value pairs.
            Returns empty dict if no metadata exists.
        """
        metadata = self.redis.hgetall(self.key)
        return metadata if metadata else {}

    def get(self, field: str) -> str | None:
        """Get specific metadata field.

        Parameters
        ----------
        field : str
            Metadata field name

        Returns
        -------
        str | None
            Field value or None if field doesn't exist
        """
        return self.redis.hget(self.key, field)

    def set(self, field: str, value: str) -> None:
        """Set a metadata field.

        Parameters
        ----------
        field : str
            Metadata field name
        value : str
            Field value (must be string)
        """
        if not isinstance(value, str):
            raise ValueError(f"Value must be string, got {type(value).__name__}")
        self.redis.hset(self.key, field, value)

    def update(self, data: dict[str, str]) -> None:
        """Update multiple metadata fields atomically.

        Parameters
        ----------
        data : dict[str, str]
            Dictionary of field-value pairs to update

        Raises
        ------
        ValueError
            If any value is not a string
        """
        if not data:
            return

        # Validate all values are strings
        for key, value in data.items():
            if not isinstance(value, str):
                raise ValueError(
                    f"All values must be strings. Field '{key}' has type {type(value).__name__}"
                )

        # Update atomically
        self.redis.hset(self.key, mapping=data)

    def delete(self, field: str) -> bool:
        """Delete a metadata field.

        Parameters
        ----------
        field : str
            Metadata field name to delete

        Returns
        -------
        bool
            True if field was deleted, False if field didn't exist
        """
        result = self.redis.hdel(self.key, field)
        return result > 0

    def clear(self) -> None:
        """Clear all metadata for the room."""
        self.redis.delete(self.key)

    def exists(self) -> bool:
        """Check if metadata exists for the room.

        Returns
        -------
        bool
            True if any metadata exists, False otherwise
        """
        return self.redis.exists(self.key) > 0
