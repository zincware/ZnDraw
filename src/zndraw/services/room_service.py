"""Room lifecycle and data management.

This service handles room creation, validation, and data management
following SOLID principles.

Note: Pylance may show type errors for Redis operations suggesting they
return Awaitables. These are false positives - we use synchronous Redis.
"""

import json
import logging
import re

from redis import Redis  # type: ignore

log = logging.getLogger(__name__)


class RoomService:
    """Handles room creation, validation, and data management.

    Follows SOLID principles:
    - Single Responsibility: Only manages room lifecycle
    - Open/Closed: Extensible through inheritance
    - Interface Segregation: Focused interface

    Parameters
    ----------
    redis_client : Redis
        Redis client instance for data operations
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def room_exists(self, room_id: str) -> bool:
        """Check if room exists by checking for current_frame key.

        Parameters
        ----------
        room_id : str
            Unique room identifier

        Returns
        -------
        bool
            True if room exists, False otherwise
        """
        return self.r.exists(f"room:{room_id}:current_frame") > 0

    def create_room(
        self,
        room_id: str,
        user_name: str,
        description: str | None = None,
        copy_from: str | None = None,
    ) -> dict:
        """Create new room, optionally copying from existing room.

        Parameters
        ----------
        room_id : str
            Unique room identifier. Allowed characters: alphanumeric (a-Z, 0-9),
            underscores (_), hyphens (-), and dots (.). Colons (:) are not allowed
            as they are used internally for Redis key namespacing.
        user_name : str
            Username creating the room (for default settings)
        description : str | None
            Optional room description
        copy_from : str | None
            Optional source room to copy from

        Returns
        -------
        dict
            {"created": bool, "frameCount": int}

        Raises
        ------
        ValueError
            If room_id contains invalid characters or copy_from room doesn't exist
        """
        # Validate room ID
        if ":" in room_id:
            raise ValueError("Room ID cannot contain ':' character")

        # Reject any non-string/int room IDs (allow alphanumeric, underscores, hyphens, and dots)
        if not re.match(r"^[a-zA-Z0-9_\-\.]+$", room_id):
            raise ValueError("Room ID contains invalid characters")

        if copy_from:
            return self._create_room_from_copy(room_id, copy_from)
        return self._create_empty_room(room_id, user_name, description)

    def _create_empty_room(
        self, room_id: str, user_name: str, description: str | None = None
    ) -> dict:
        """Create empty room with defaults.

        Uses Redis pipeline for atomic multi-operation setup.

        Note: User-specific settings are NOT initialized here (YAGNI).
        Settings will be initialized lazily when first accessed.

        Parameters
        ----------
        room_id : str
            Unique room identifier
        user_name : str
            Username creating the room
        description : str | None
            Optional room description

        Returns
        -------
        dict
            {"created": True, "frameCount": 0}
        """
        pipe = self.r.pipeline()

        # Set description if provided
        if description:
            pipe.set(f"room:{room_id}:description", description)

        # Initialize metadata
        pipe.set(f"room:{room_id}:current_frame", 0)
        pipe.set(f"room:{room_id}:locked", 0)
        pipe.set(f"room:{room_id}:hidden", 0)

        # Create default geometries
        self._initialize_default_geometries_pipeline(room_id, pipe)

        # Execute all operations atomically
        pipe.execute()

        log.info(f"Created empty room '{room_id}'")
        return {"created": True, "frameCount": 0}

    def _create_room_from_copy(self, room_id: str, source_room: str) -> dict:
        """Copy room data from existing room.

        Uses Redis pipeline for efficient bulk copying.

        Parameters
        ----------
        room_id : str
            New room identifier
        source_room : str
            Source room to copy from

        Returns
        -------
        dict
            {"created": True, "frameCount": int}

        Raises
        ------
        ValueError
            If source room doesn't exist
        """
        source_indices_key = f"room:{source_room}:trajectory:indices"
        if not self.r.exists(source_indices_key):
            raise ValueError(f"Source room '{source_room}' not found")

        pipe = self.r.pipeline()

        # Copy trajectory indices (shares frame data)
        source_indices = self.r.zrange(source_indices_key, 0, -1, withscores=True)
        if source_indices:
            pipe.zadd(
                f"room:{room_id}:trajectory:indices",
                {member: score for member, score in source_indices},
            )

        # Copy geometries
        geometries = self.r.hgetall(f"room:{source_room}:geometries")
        if geometries:
            pipe.hset(f"room:{room_id}:geometries", mapping=geometries)

        # Copy bookmarks
        bookmarks = self.r.hgetall(f"room:{source_room}:bookmarks")
        if bookmarks:
            pipe.hset(f"room:{room_id}:bookmarks", mapping=bookmarks)

        # Initialize metadata
        pipe.set(f"room:{room_id}:current_frame", 0)
        pipe.set(f"room:{room_id}:locked", 0)
        pipe.set(f"room:{room_id}:hidden", 0)

        # Execute all operations atomically
        pipe.execute()

        log.info(
            f"Created room '{room_id}' from '{source_room}' with {len(source_indices)} frames"
        )
        return {"created": True, "frameCount": len(source_indices)}

    def _initialize_default_geometries_pipeline(
        self, room_id: str, pipe: Redis.pipeline
    ):
        """Initialize default geometries for new room using pipeline.

        Parameters
        ----------
        room_id : str
            Room identifier
        pipe : Redis.pipeline
            Redis pipeline to add operations to
        """
        from zndraw.geometries import Bond, Cell, Curve, Floor, Sphere

        defaults = {
            "particles": (
                Sphere,
                {"position": "arrays.positions", "color": "arrays.colors", "radius": "arrays.radii"}
            ),
            "bonds": (Bond, {"position": "arrays.positions", "color": "arrays.colors"}),
            "curve": (Curve, {}),
            "cell": (Cell, {}),
            "floor": (Floor, {}),
        }

        for key, (geometry_class, kwargs) in defaults.items():
            geometry_data = geometry_class(**kwargs).model_dump()
            pipe.hset(
                f"room:{room_id}:geometries",
                key,
                json.dumps({"type": geometry_class.__name__, "data": geometry_data}),
            )

    def get_frame_count(self, room_id: str) -> int:
        """Get number of frames in room.

        Parameters
        ----------
        room_id : str
            Room identifier

        Returns
        -------
        int
            Number of frames in the room's trajectory
        """
        return self.r.zcard(f"room:{room_id}:trajectory:indices")

    def get_current_frame(self, room_id: str) -> int:
        """Get current frame number, handling invalid values.

        Parameters
        ----------
        room_id : str
            Room identifier

        Returns
        -------
        int
            Current frame number (0 if invalid or not set)
        """
        step = self.r.get(f"room:{room_id}:current_frame")
        try:
            if step is not None:
                step_int = int(step)
                if step_int < 0:
                    log.warning(
                        f"Negative frame in Redis for room {room_id}: {step_int}, resetting to 0"
                    )
                    return 0
                return step_int
            return 0
        except (ValueError, TypeError) as e:
            log.error(f"Invalid step value in Redis for room {room_id}: {step} - {e}")
            return 0
