"""Client session metadata management (JWT-aware).

Note: JWT tokens are created/validated by auth.py.
This service only manages client metadata in Redis.
"""

import logging

from redis import Redis

log = logging.getLogger(__name__)


class ClientService:
    """Handles client session metadata in Redis.

    JWT authentication is handled separately by auth.py.
    This service focuses on storing client metadata like currentRoom.

    Parameters
    ----------
    redis_client : Redis
        Redis client instance for data operations
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def update_user_room(self, user_name: str, room_id: str) -> None:
        """Update user's current room in Redis.

        Parameters
        ----------
        user_name : str
            User name (from JWT sub claim)
        room_id : str
            Room the user is joining
        """
        user_key = f"user:{user_name}"
        self.r.hset(user_key, "currentRoom", room_id)
        log.info(f"User {user_name} updated room to {room_id}")

    def add_user_to_room(self, room_id: str, user_name: str) -> None:
        """Add user to room's user set.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            User name
        """
        self.r.sadd(f"room:{room_id}:users", user_name)

    def remove_user_from_room(self, room_id: str, user_name: str) -> None:
        """Remove user from room's user set.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            User name
        """
        self.r.srem(f"room:{room_id}:users", user_name)

    def get_room_users(self, room_id: str) -> set[str]:
        """Get all userNames currently in a room.

        Parameters
        ----------
        room_id : str
            Room identifier

        Returns
        -------
        set[str]
            Set of userNames
        """
        members = self.r.smembers(f"room:{room_id}:users")
        return {m.decode() if isinstance(m, bytes) else m for m in members}

    def update_user_and_room_membership(self, user_name: str, room_id: str) -> None:
        """Update user room and add to room membership atomically.

        Uses Redis pipeline for atomic operations.

        Parameters
        ----------
        user_name : str
            User name
        room_id : str
            Room identifier
        """
        pipe = self.r.pipeline()
        pipe.hset(f"user:{user_name}", "currentRoom", room_id)
        pipe.sadd(f"room:{room_id}:users", user_name)
        pipe.execute()
        log.info(f"User {user_name} joined room {room_id}")
