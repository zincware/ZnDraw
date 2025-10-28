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

    def update_client_room(self, client_id: str, room_id: str) -> None:
        """Update client's current room in Redis.

        Parameters
        ----------
        client_id : str
            Client identifier (from JWT sub claim)
        room_id : str
            Room the client is joining
        """
        client_key = f"client:{client_id}"
        self.r.hset(client_key, "currentRoom", room_id)
        log.info(f"Client {client_id} updated room to {room_id}")

    def add_client_to_room(self, room_id: str, client_id: str) -> None:
        """Add client to room's client set.

        Parameters
        ----------
        room_id : str
            Room identifier
        client_id : str
            Client identifier
        """
        self.r.sadd(f"room:{room_id}:clients", client_id)

    def remove_client_from_room(self, room_id: str, client_id: str) -> None:
        """Remove client from room's client set.

        Parameters
        ----------
        room_id : str
            Room identifier
        client_id : str
            Client identifier
        """
        self.r.srem(f"room:{room_id}:clients", client_id)

    def get_room_clients(self, room_id: str) -> set[str]:
        """Get all client IDs currently in a room.

        Parameters
        ----------
        room_id : str
            Room identifier

        Returns
        -------
        set[str]
            Set of client IDs
        """
        members = self.r.smembers(f"room:{room_id}:clients")
        return {m.decode() if isinstance(m, bytes) else m for m in members}

    def update_client_and_room_membership(self, client_id: str, room_id: str) -> None:
        """Update client room and add to room membership atomically.

        Uses Redis pipeline for atomic operations.

        Parameters
        ----------
        client_id : str
            Client identifier
        room_id : str
            Room identifier
        """
        pipe = self.r.pipeline()
        pipe.hset(f"client:{client_id}", "currentRoom", room_id)
        pipe.sadd(f"room:{room_id}:clients", client_id)
        pipe.execute()
        log.info(f"Client {client_id} joined room {room_id}")
