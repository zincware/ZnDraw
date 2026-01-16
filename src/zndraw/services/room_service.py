"""Room lifecycle and data management.

This service handles room creation, validation, and data management
following SOLID principles.

Note: Pylance may show type errors for Redis operations suggesting they
return Awaitables. These are false positives - we use synchronous Redis.
"""

import json
import logging
import re

import ase
from asebytes import encode
from redis import Redis  # type: ignore

from zndraw.app.frame_index_manager import FrameIndexManager
from zndraw.app.redis_keys import GlobalIndexKeys, RoomKeys
from zndraw.room_templates import TEMPLATES
from zndraw.storage import StorageBackend

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
        keys = RoomKeys(room_id)
        return self.r.exists(keys.current_frame()) > 0

    def validate_room_available(self, room_id: str) -> None:
        """Validate that a room ID is available (doesn't already exist).

        Parameters
        ----------
        room_id : str
            Room identifier to check

        Raises
        ------
        ValueError
            If room already exists
        """
        if self.room_exists(room_id):
            raise ValueError(f"Room '{room_id}' already exists")

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
            If room_id contains invalid characters, is reserved, or copy_from
            room doesn't exist
        """
        # Validate room ID using is_valid_room_id()
        valid, error = self.is_valid_room_id(room_id)
        if not valid:
            raise ValueError(error)

        if copy_from:
            return self._create_room_from_copy(room_id, copy_from, description)
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
        from zndraw.app.redis_keys import GlobalIndexKeys

        keys = RoomKeys(room_id)
        pipe = self.r.pipeline()

        # Set description if provided
        if description:
            pipe.set(keys.description(), description)

        # Initialize metadata
        pipe.set(keys.current_frame(), 0)
        pipe.set(keys.locked(), 0)

        # Add room to global index for O(1) room listing
        pipe.sadd(GlobalIndexKeys.rooms_index(), room_id)

        # Create default geometries
        self._initialize_default_geometries_pipeline(room_id, pipe)

        # Execute all operations atomically
        pipe.execute()

        log.debug(f"Created empty room '{room_id}'")
        return {"created": True, "frameCount": 0}

    def _create_room_from_copy(
        self, room_id: str, source_room: str, description: str | None = None
    ) -> dict:
        """Copy room data from existing room.

        Uses Redis pipeline for efficient bulk copying.

        Parameters
        ----------
        room_id : str
            New room identifier
        source_room : str
            Source room to copy from
        description : str | None
            Optional description for the new room

        Returns
        -------
        dict
            {"created": True, "frameCount": int}

        Raises
        ------
        ValueError
            If source room doesn't exist
        """
        from zndraw.app.redis_keys import GlobalIndexKeys

        source_keys = RoomKeys(source_room)
        new_keys = RoomKeys(room_id)

        if not self.r.exists(source_keys.trajectory_indices()):
            raise ValueError(f"Source room '{source_room}' not found")

        pipe = self.r.pipeline()

        # Set description if provided
        if description:
            pipe.set(new_keys.description(), description)

        # Copy trajectory indices (shares frame data)
        source_indices = self.r.zrange(
            source_keys.trajectory_indices(), 0, -1, withscores=True
        )
        if source_indices:
            pipe.zadd(
                new_keys.trajectory_indices(),
                {member: score for member, score in source_indices},
            )

        # Copy geometries
        geometries = self.r.hgetall(source_keys.geometries())
        if geometries:
            pipe.hset(new_keys.geometries(), mapping=geometries)

        # Copy bookmarks
        bookmarks = self.r.hgetall(source_keys.bookmarks())
        if bookmarks:
            pipe.hset(new_keys.bookmarks(), mapping=bookmarks)

        # Initialize metadata
        pipe.set(new_keys.current_frame(), 0)
        pipe.set(new_keys.locked(), 0)

        # Add room to global index for O(1) room listing
        pipe.sadd(GlobalIndexKeys.rooms_index(), room_id)

        # Execute all operations atomically
        pipe.execute()

        log.debug(
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
        from zndraw.geometries import (
            Bond,
            Cell,
            Curve,
            Floor,
            InteractionSettings,
            Sphere,
        )
        from zndraw.materials import MeshBasicMaterial
        from zndraw.transformations import InArrayTransform

        keys = RoomKeys(room_id)

        defaults = {
            "particles": (
                Sphere,
                {
                    "position": "arrays.positions",
                    "color": "arrays.colors",
                    "radius": "arrays.radii",
                    "scale": [(0.7, 0.7, 0.7)],
                },
            ),
            "bonds": (
                Bond,
                {
                    "position": "arrays.positions",
                    "color": "arrays.colors",
                    "scale": 0.15,
                },
            ),
            "curve": (Curve, {}),
            "cell": (Cell, {}),
            "floor": (Floor, {}),
            "constraints-fixed-atoms": (
                Sphere,
                {
                    "position": InArrayTransform(
                        source="constraints",
                        path="0.kwargs.indices",
                        filter="arrays.positions",
                    ),
                    "radius": InArrayTransform(
                        source="constraints",
                        path="0.kwargs.indices",
                        filter="arrays.radii",
                    ),
                    "color": ["#FF0000"],
                    "material": MeshBasicMaterial(wireframe=True),
                    "scale": [(0.71, 0.71, 0.71)],  # Larger to be clearly visible
                    "active": True,  # Active by default to visualize constraints
                    "selecting": InteractionSettings(enabled=False),
                    "hovering": InteractionSettings(enabled=False),
                },
            ),
        }

        for key, (geometry_class, kwargs) in defaults.items():
            geometry_data = geometry_class(**kwargs).model_dump()
            pipe.hset(
                keys.geometries(),
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
        keys = RoomKeys(room_id)
        return self.r.zcard(keys.trajectory_indices())

    def get_frame_counts_batch(self, room_ids: list[str]) -> dict[str, int]:
        """Get frame counts for multiple rooms efficiently.

        Parameters
        ----------
        room_ids : list[str]
            List of room identifiers

        Returns
        -------
        dict[str, int]
            Mapping of room_id to frame count

        Notes
        -----
        Uses Redis pipeline to batch ZCARD queries for better performance
        when fetching counts for many rooms.
        """
        if not room_ids:
            return {}

        # Use pipeline to batch ZCARD queries
        pipe = self.r.pipeline()
        for room_id in room_ids:
            keys = RoomKeys(room_id)
            pipe.zcard(keys.trajectory_indices())

        results = pipe.execute()

        # Build mapping
        return {room_id: count for room_id, count in zip(room_ids, results)}

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
        keys = RoomKeys(room_id)
        step = self.r.get(keys.current_frame())
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

    def is_valid_room_id(self, room_id: str) -> tuple[bool, str | None]:
        """Check if a room ID is valid.

        Parameters
        ----------
        room_id : str
            Room identifier to validate

        Returns
        -------
        tuple[bool, str | None]
            (True, None) if valid, (False, error_message) if invalid
        """
        if not room_id:
            return False, "Room ID cannot be empty"

        if ":" in room_id:
            return False, "Room ID cannot contain ':' character"

        if not re.match(r"^[a-zA-Z0-9_\-\.]+$", room_id):
            return False, "Room ID contains invalid characters"

        if room_id in TEMPLATES:
            return False, f"Room ID '{room_id}' is reserved (template name)"

        return True, None

    def apply_template(
        self, room_id: str, template_name: str, storage: StorageBackend
    ) -> None:
        """Apply a template to populate a room with initial frames.

        Parameters
        ----------
        room_id : str
            Room to populate
        template_name : str
            Template name: "empty" (1 frame) or "none" (0 frames)
        storage : StorageBackend
            Storage backend for frame data

        Raises
        ------
        ValueError
            If template name is not found
        """
        if template_name not in TEMPLATES:
            raise ValueError(
                f"Template '{template_name}' not found. Available: {list(TEMPLATES.keys())}"
            )

        writer = _RoomWriter(room_id, self.r, storage)
        TEMPLATES[template_name](writer)

    def get_default_room(self) -> str | None:
        """Get the admin-configured default room.

        Returns
        -------
        str | None
            Default room ID if set, None otherwise
        """
        return self.r.get("default_room")

    def create_room_with_defaults(
        self,
        room_id: str,
        user_name: str,
        storage: StorageBackend,
        description: str | None = None,
        copy_from: str | None = None,
        template: str | None = None,
    ) -> dict:
        """Create a room with server-side fallback logic.

        Fallback order:
        1. If copy_from specified → copy from that room
        2. Else if template specified → use that template
        3. Else if default_room is set → copy from default room
        4. Else → use "empty" template (1 empty frame)

        Parameters
        ----------
        room_id : str
            Unique room identifier
        user_name : str
            Username creating the room
        storage : StorageBackend
            Storage backend for frame data
        description : str | None
            Optional room description
        copy_from : str | None
            Optional source room to copy from (takes precedence)
        template : str | None
            Optional template name: "empty" (1 frame) or "none" (0 frames)

        Returns
        -------
        dict
            {"created": bool, "frameCount": int, "source": str}

        Raises
        ------
        ValueError
            If room_id is invalid, copy_from doesn't exist, or template not found
        """
        # Validate room ID
        valid, error = self.is_valid_room_id(room_id)
        if not valid:
            raise ValueError(error)

        if copy_from:
            # Explicit copyFrom takes precedence
            result = self.create_room(room_id, user_name, description, copy_from)
            result["source"] = f"room:{copy_from}"
            return result

        if template:
            # Explicit template takes precedence
            if template not in TEMPLATES:
                raise ValueError(
                    f"Template '{template}' not found. Available: {list(TEMPLATES.keys())}"
                )
            result = self.create_room(room_id, user_name, description)
            self.apply_template(room_id, template, storage)
            result["frameCount"] = self.get_frame_count(room_id)
            result["source"] = f"template:{template}"
            return result

        # No explicit source - check default_room
        default_room = self.get_default_room()
        if default_room and self.room_exists(default_room):
            result = self.create_room(room_id, user_name, description, default_room)
            result["source"] = f"room:{default_room}"
            log.debug(f"Created room '{room_id}' from default room '{default_room}'")
            return result

        # Final fallback: empty template
        result = self.create_room(room_id, user_name, description)
        self.apply_template(room_id, "empty", storage)
        result["frameCount"] = self.get_frame_count(room_id)
        result["source"] = "template:empty"
        return result

    def delete_room(self, room_id: str) -> bool:
        """Delete a room and all its data.

        Removes the room from the global index and deletes all room-related keys.

        Parameters
        ----------
        room_id : str
            Room identifier to delete

        Returns
        -------
        bool
            True if room was deleted, False if room didn't exist
        """
        keys = RoomKeys(room_id)

        # Check if room exists
        if not self.room_exists(room_id):
            return False

        # Delete all known room keys directly (no SCAN needed)
        pipe = self.r.pipeline()
        pipe.srem(GlobalIndexKeys.rooms_index(), room_id)
        pipe.delete(*keys.all_static_keys())
        pipe.execute()

        log.debug(f"Deleted room '{room_id}' and removed from index")
        return True


class _RoomWriter:
    """Minimal interface for templates to write to a room's storage.

    This class provides an append/extend interface that templates use to
    add frames to a room without going through the full ZnDraw client.

    Parameters
    ----------
    room_id : str
        Room identifier
    redis_client : Redis
        Redis client instance
    storage : StorageBackend
        Storage backend for frame data
    """

    def __init__(self, room_id: str, redis_client: Redis, storage: StorageBackend):
        self.room_id = room_id
        self.r = redis_client
        self.storage = storage

    def append(self, atoms) -> None:
        """Append a single frame.

        Parameters
        ----------
        atoms : ase.Atoms
            Atoms object to append
        """
        self.extend([atoms])

    def extend(self, frames: list[ase.Atoms]) -> None:
        """Extend with multiple frames.

        Parameters
        ----------
        frames : list[ase.Atoms]
            List of Atoms objects to add
        """
        room_keys = RoomKeys(self.room_id)
        manager = FrameIndexManager(self.r, room_keys.trajectory_indices())

        for atoms in frames:
            encoded = encode(atoms)
            physical_index = len(self.storage)
            self.storage.extend([encoded])
            manager.append(f"{self.room_id}:{physical_index}")
