"""Client session metadata management (JWT-aware).

Note: JWT tokens are created/validated by auth.py.
This service only manages client metadata in Redis.
"""

import logging

from redis import Redis

from zndraw.app.redis_keys import RoomKeys, UserKeys

log = logging.getLogger(__name__)


class ClientService:
    """Handles client session metadata in Redis.

    JWT authentication is handled separately by auth.py.
    This service manages room membership and visit tracking.

    Parameters
    ----------
    redis_client : Redis
        Redis client instance for data operations
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def remove_user_from_room(self, room_id: str, user_name: str) -> None:
        """Remove user from room's user set.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            User name
        """
        self.r.srem(RoomKeys(room_id).users(), user_name)

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
        members = self.r.smembers(RoomKeys(room_id).users())
        return {m.decode() if isinstance(m, bytes) else m for m in members}

    def update_user_and_room_membership(self, user_name: str, room_id: str) -> None:
        """Add user to room membership and track visit atomically.

        Uses Redis pipeline for atomic operations.

        Parameters
        ----------
        user_name : str
            User name
        room_id : str
            Room identifier
        """
        keys = UserKeys(user_name)
        pipe = self.r.pipeline()
        pipe.sadd(RoomKeys(room_id).users(), user_name)
        pipe.sadd(keys.visited_rooms(), room_id)
        pipe.execute()
        log.debug(f"User {user_name} joined room {room_id}")

    def get_visited_rooms(self, user_name: str) -> set[str]:
        """Get all room IDs that a user has visited.

        Parameters
        ----------
        user_name : str
            User name

        Returns
        -------
        set[str]
            Set of room IDs the user has visited
        """
        keys = UserKeys(user_name)
        members = self.r.smembers(keys.visited_rooms())
        return {m.decode() if isinstance(m, bytes) else m for m in members}
