"""Redis key management for ZnDraw.

Centralizes Redis key construction to avoid duplication and errors.
All key construction should go through these classes.
"""

from dataclasses import dataclass
from typing import ClassVar


@dataclass
class ExtensionKeys:
    """Redis keys for an extension.

    Centralizes key construction to avoid duplication and errors.
    """

    schema: str
    workers: str
    pending_jobs: str

    @classmethod
    def for_extension(
        cls, room_id: str, category: str, extension: str
    ) -> "ExtensionKeys":
        """Create ExtensionKeys for a specific extension.

        Args:
            room_id: The room identifier
            category: The extension category (e.g., 'modifiers', 'selections')
            extension: The extension name

        Returns:
            ExtensionKeys instance with all relevant Redis keys
        """
        base = f"room:{room_id}:extensions:{category}"
        return cls(
            schema=base,
            workers=f"{base}:{extension}:workers",
            pending_jobs=f"{base}:{extension}:pending_jobs",
        )

    @staticmethod
    def schema_key(room_id: str, category: str) -> str:
        """Get the schema hash key for a room and category.

        This is the Redis hash that stores all extension schemas for a category.

        Args:
            room_id: The room identifier
            category: The extension category

        Returns:
            Redis key for the schema hash
        """
        return f"room:{room_id}:extensions:{category}"

    @staticmethod
    def user_extensions_key(room_id: str, category: str, sid: str) -> str:
        """Get the reverse mapping key for a worker's registered extensions.

        This key maps from a session ID to the set of extension names
        that the worker can handle.

        Args:
            room_id: The room identifier
            category: The extension category
            sid: The session ID of the worker

        Returns:
            Redis key for the user's extensions set
        """
        return f"room:{room_id}:extensions:{category}:{sid}"

    @classmethod
    def for_global_extension(cls, category: str, extension: str) -> "ExtensionKeys":
        """Create ExtensionKeys for a global extension (not room-specific).

        Args:
            category: The extension category (e.g., 'modifiers', 'selections')
            extension: The extension name

        Returns:
            ExtensionKeys instance with all relevant Redis keys for global extensions
        """
        base = f"global:extensions:{category}"
        return cls(
            schema=base,
            workers=f"{base}:{extension}:workers",
            pending_jobs=f"{base}:{extension}:pending_jobs",
        )

    @staticmethod
    def global_schema_key(category: str) -> str:
        """Get the schema hash key for global extensions in a category.

        This is the Redis hash that stores all global extension schemas for a category.

        Args:
            category: The extension category

        Returns:
            Redis key for the global schema hash
        """
        return f"global:extensions:{category}"

    @staticmethod
    def global_user_extensions_key(category: str, sid: str) -> str:
        """Get the reverse mapping key for a worker's registered global extensions.

        This key maps from a session ID to the set of global extension names
        that the worker can handle.

        Args:
            category: The extension category
            sid: The session ID of the worker

        Returns:
            Redis key for the user's global extensions set
        """
        return f"global:extensions:{category}:{sid}"

    @staticmethod
    def worker_capacity_key(worker_id: str) -> str:
        """Get the worker capacity key for tracking available job slots.

        Args:
            worker_id: The worker identifier (socket sid or celery task id)

        Returns:
            Redis key for the worker's available capacity (INTEGER)
        """
        return f"worker:{worker_id}:available_capacity"

    @staticmethod
    def worker_active_jobs_key(worker_id: str) -> str:
        """Get the worker active jobs key for tracking currently assigned jobs.

        Args:
            worker_id: The worker identifier (socket sid or celery task id)

        Returns:
            Redis key for the worker's active jobs (SET)
        """
        return f"worker:{worker_id}:active_jobs"


@dataclass
class FilesystemKeys:
    """Redis keys for a filesystem provider.

    Centralizes key construction for filesystem-related Redis keys.
    """

    metadata: str
    worker: str

    @classmethod
    def for_filesystem(cls, room_id: str, fs_name: str) -> "FilesystemKeys":
        """Create FilesystemKeys for a specific filesystem.

        Parameters
        ----------
        room_id : str
            The room identifier
        fs_name : str
            The filesystem name

        Returns
        -------
        FilesystemKeys
            Instance with all relevant Redis keys
        """
        base = f"room:{room_id}:filesystems:{fs_name}"
        return cls(
            metadata=base,
            worker=f"{base}:worker",
        )

    @classmethod
    def for_global_filesystem(cls, fs_name: str) -> "FilesystemKeys":
        """Create FilesystemKeys for a global filesystem (not room-specific).

        Parameters
        ----------
        fs_name : str
            The filesystem name

        Returns
        -------
        FilesystemKeys
            Instance with all relevant Redis keys for global filesystems
        """
        base = f"global:filesystems:{fs_name}"
        return cls(
            metadata=base,
            worker=f"{base}:worker",
        )

    @staticmethod
    def user_filesystems_key(room_id: str, sid: str) -> str:
        """Get the reverse mapping key for a worker's registered filesystems.

        This key maps from a session ID to the set of filesystem names
        that the worker provides.

        Parameters
        ----------
        room_id : str
            The room identifier
        sid : str
            The session ID of the worker

        Returns
        -------
        str
            Redis key for the user's filesystems set
        """
        return f"room:{room_id}:filesystems:{sid}"

    @staticmethod
    def global_user_filesystems_key(sid: str) -> str:
        """Get the reverse mapping key for a worker's registered global filesystems.

        This key maps from a session ID to the set of global filesystem names
        that the worker provides.

        Parameters
        ----------
        sid : str
            The session ID of the worker

        Returns
        -------
        str
            Redis key for the user's global filesystems set
        """
        return f"global:filesystems:{sid}"


@dataclass(frozen=True)
class RoomKeys:
    """Redis keys for room-related data.

    Centralizes all room-scoped key construction.
    """

    room_id: str

    def description(self) -> str:
        """Room description string."""
        return f"room:{self.room_id}:description"

    def locked(self) -> str:
        """Room locked status (0 or 1)."""
        return f"room:{self.room_id}:locked"

    def locked_by(self) -> str:
        """Username of admin who locked the room (or None for non-admin locks)."""
        return f"room:{self.room_id}:locked_by"

    def current_frame(self) -> str:
        """Current frame index."""
        return f"room:{self.room_id}:current_frame"

    def presenter_lock(self) -> str:
        """Presenter lock (sid with TTL)."""
        return f"room:{self.room_id}:presenter_lock"

    def lock(self, target: str) -> str:
        """Generic lock key for a target.

        Parameters
        ----------
        target : str
            Lock target identifier

        Returns
        -------
        str
            Redis key for the lock
        """
        return f"room:{self.room_id}:lock:{target}"

    def lock_metadata(self, target: str) -> str:
        """Lock metadata for a target.

        Parameters
        ----------
        target : str
            Lock target identifier

        Returns
        -------
        str
            Redis key for lock metadata
        """
        return f"room:{self.room_id}:lock:{target}:metadata"

    def bookmarks(self) -> str:
        """Bookmarks hash (frame_index -> label)."""
        return f"room:{self.room_id}:bookmarks"

    def users(self) -> str:
        """Set of usernames currently in this room."""
        return f"room:{self.room_id}:users"

    def geometries(self) -> str:
        """Geometries hash."""
        return f"room:{self.room_id}:geometries"

    def figures(self) -> str:
        """Figures hash."""
        return f"room:{self.room_id}:figures"

    def selections(self) -> str:
        """Selections hash (geometry -> indices)."""
        return f"room:{self.room_id}:selections"

    def selection_groups(self) -> str:
        """Selection groups hash."""
        return f"room:{self.room_id}:selection_groups"

    def active_selection_group(self) -> str:
        """Active selection group name."""
        return f"room:{self.room_id}:active_selection_group"

    def session_settings(self, session_id: str) -> str:
        """Session-scoped rendering settings.

        Replaces user-scoped settings to allow different settings per browser window.

        Parameters
        ----------
        session_id : str
            The session identifier

        Returns
        -------
        str
            Redis key for session settings hash
        """
        return f"room:{self.room_id}:session_settings:{session_id}"

    def session_active_camera(self, session_id: str) -> str:
        """Active camera key for a session.

        Stores which camera geometry key the session is viewing through.

        Parameters
        ----------
        session_id : str
            The session identifier

        Returns
        -------
        str
            Redis key for session's active camera
        """
        return f"room:{self.room_id}:session:{session_id}:active_camera"

    def frontend_sessions(self) -> str:
        """Set of frontend session IDs in this room.

        Tracks which sessions are frontend (browser) vs Python clients.

        Returns
        -------
        str
            Redis key for frontend sessions set
        """
        return f"room:{self.room_id}:frontend_sessions"

    def trajectory_indices(self) -> str:
        """Trajectory indices sorted set."""
        return f"room:{self.room_id}:trajectory:indices"

    def chat_counter(self) -> str:
        """Chat message counter for generating IDs."""
        return f"room:{self.room_id}:chat:counter"

    def chat_data(self) -> str:
        """Chat message data hash."""
        return f"room:{self.room_id}:chat:data"

    def chat_index(self) -> str:
        """Chat message index (sorted set)."""
        return f"room:{self.room_id}:chat:index"

    def chat_message(self, message_id: str) -> str:
        """Individual chat message hash.

        Parameters
        ----------
        message_id : str
            Message identifier

        Returns
        -------
        str
            Redis key for message hash
        """
        return f"room:{self.room_id}:chat:message:{message_id}"

    def progress(self) -> str:
        """Active progress tracking hash (progressId -> progress data)."""
        return f"room:{self.room_id}:progress"

    def frame_selection(self, group: str = "default") -> str:
        """Frame selection for a specific group.

        Parameters
        ----------
        group : str
            Selection group name (default: "default")

        Returns
        -------
        str
            Redis key for frame selection
        """
        return f"room:{self.room_id}:frame_selection:{group}"

    def jobs_active(self) -> str:
        """Active jobs set (queued or running)."""
        return f"room:{self.room_id}:jobs:active"

    def jobs_inactive(self) -> str:
        """Inactive jobs set (completed or failed)."""
        return f"room:{self.room_id}:jobs:inactive"

    def jobs_by_time(self) -> str:
        """Jobs sorted set by timestamp."""
        return f"room:{self.room_id}:jobs:by_time"

    def extension_jobs(self, category: str, extension: str) -> str:
        """Jobs set for a specific extension.

        Parameters
        ----------
        category : str
            Extension category
        extension : str
            Extension name

        Returns
        -------
        str
            Redis key for extension jobs set
        """
        return f"room:{self.room_id}:extension:{category}:{extension}:jobs"

    def user_extension_data(self, username: str, category: str) -> str:
        """User-specific extension data hash.

        Parameters
        ----------
        username : str
            Username
        category : str
            Extension category

        Returns
        -------
        str
            Redis key for user extension data hash
        """
        return f"room:{self.room_id}:user:{username}:{category}"

    def room_extension_data(self, category: str) -> str:
        """Room-wide extension data hash (shared across all users).

        Used for settings and other data that should be synchronized across users.

        Parameters
        ----------
        category : str
            Extension category

        Returns
        -------
        str
            Redis key for room extension data hash
        """
        return f"room:{self.room_id}:extension:{category}"

    def metadata(self) -> str:
        """Room metadata hash.

        Returns
        -------
        str
            Redis key for room metadata hash
        """
        return f"room:{self.room_id}:metadata"

    def filesystems_pattern(self) -> str:
        """Pattern for scanning all filesystems in this room.

        Returns
        -------
        str
            Redis key pattern for scanning filesystems
        """
        return f"room:{self.room_id}:filesystems:*"

    def all_keys_pattern(self) -> str:
        """Pattern for scanning all keys for this room.

        Returns
        -------
        str
            Redis key pattern for scanning all room keys
        """
        return f"room:{self.room_id}:*"

    def all_static_keys(self) -> list[str]:
        """Return all static room keys (keys without dynamic parameters).

        Use this for bulk deletion. Does not include parameterized keys like
        settings(username), lock(target), chat_message(id), session_settings(id), etc.

        Returns
        -------
        list[str]
            All static Redis keys for this room
        """
        return [
            self.description(),
            self.locked(),
            self.locked_by(),
            self.current_frame(),
            self.presenter_lock(),
            self.bookmarks(),
            self.geometries(),
            self.figures(),
            self.selections(),
            self.selection_groups(),
            self.active_selection_group(),
            self.trajectory_indices(),
            self.chat_counter(),
            self.chat_data(),
            self.chat_index(),
            self.progress(),
            self.metadata(),
            self.jobs_active(),
            self.jobs_inactive(),
            self.jobs_by_time(),
            self.frontend_sessions(),
        ]


@dataclass(frozen=True)
class UserKeys:
    """Redis keys for user-related data.

    Key structure uses `users:{type}:{username}` format for efficient scanning:
    - `users:data:{username}` - user data hash
    - `users:admin:{username}` - admin status
    - `users:visited:{username}` - visited rooms set
    """

    # Key prefixes (single source of truth) - ClassVar so not treated as fields
    DATA_PREFIX: ClassVar[str] = "users:data:"
    ADMIN_PREFIX: ClassVar[str] = "users:admin:"
    VISITED_PREFIX: ClassVar[str] = "users:visited:"

    username: str

    def hash_key(self) -> str:
        """User data hash containing all user fields."""
        return f"{self.DATA_PREFIX}{self.username}"

    def admin_key(self) -> str:
        """Admin status key for this user."""
        return f"{self.ADMIN_PREFIX}{self.username}"

    def visited_rooms(self) -> str:
        """Set of room IDs this user has visited."""
        return f"{self.VISITED_PREFIX}{self.username}"

    @staticmethod
    def data_pattern() -> str:
        """Pattern for scanning all user data keys."""
        return f"{UserKeys.DATA_PREFIX}*"

    @staticmethod
    def admin_pattern() -> str:
        """Pattern for scanning all admin keys."""
        return f"{UserKeys.ADMIN_PREFIX}*"

    @staticmethod
    def username_from_data_key(key: str) -> str:
        """Extract username from a data key."""
        return key.removeprefix(UserKeys.DATA_PREFIX)

    @staticmethod
    def username_from_admin_key(key: str) -> str:
        """Extract username from an admin key."""
        return key.removeprefix(UserKeys.ADMIN_PREFIX)


@dataclass(frozen=True)
class SessionKeys:
    """Redis keys for session mapping."""

    sid: str

    def username(self) -> str:
        """Session ID to username mapping."""
        return f"sid:{self.sid}"

    def role(self) -> str:
        """Session role."""
        return f"sid:{self.sid}:role"

    def session_id(self) -> str:
        """Socket ID to session ID mapping (for cleanup on disconnect)."""
        return f"sid:{self.sid}:session"

    @staticmethod
    def session_data(session_id: str) -> str:
        """Session data storage key.

        Stores session metadata (userId, roomId, timestamps, etc.)

        Parameters
        ----------
        session_id : str
            The session identifier

        Returns
        -------
        str
            Redis key for session data storage
        """
        return f"session:{session_id}"

    @staticmethod
    def session_to_sid(session_id: str) -> str:
        """Session ID to socket ID reverse mapping.

        Parameters
        ----------
        session_id : str
            The session identifier

        Returns
        -------
        str
            Redis key for session ID to socket ID mapping
        """
        return f"session:{session_id}:sid"

    @staticmethod
    def session_locks(session_id: str) -> str:
        """Set of lock keys held by this session.

        Parameters
        ----------
        session_id : str
            The session identifier

        Returns
        -------
        str
            Redis key for session's lock set
        """
        return f"session:{session_id}:locks"


@dataclass(frozen=True)
class JobKeys:
    """Redis keys for job tracking."""

    job_id: str

    def hash_key(self) -> str:
        """Job data hash containing all job fields."""
        return f"job:{self.job_id}"


@dataclass(frozen=True)
class WorkerKeys:
    """Redis keys for worker-related data."""

    worker_id: str

    def active_jobs(self) -> str:
        """Set of active job IDs assigned to this worker."""
        return f"worker:{self.worker_id}:jobs"


class GlobalIndexKeys:
    """Redis keys for global indices."""

    # Room index: SET containing all room IDs
    ROOMS_INDEX = "rooms:index"

    # User indices
    USERS_INDEX = "users:index"  # SET of all usernames
    ADMINS_INDEX = "admins:index"  # SET of all admin usernames

    @staticmethod
    def rooms_index() -> str:
        """Set of all room IDs.
        When a room is created, its ID is added to this set.
        When a room is deleted, its ID is removed from this set.

        Returns
        -------
        str
            Redis key for the rooms index set
        """
        return GlobalIndexKeys.ROOMS_INDEX

    @staticmethod
    def users_index() -> str:
        """Set of all usernames.

        Returns
        -------
        str
            Redis key for the users index set
        """
        return GlobalIndexKeys.USERS_INDEX

    @staticmethod
    def admins_index() -> str:
        """Set of all admin usernames.

        Returns
        -------
        str
            Redis key for the admins index set
        """
        return GlobalIndexKeys.ADMINS_INDEX
