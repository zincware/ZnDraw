"""Constants for the ZnDraw application."""


class SocketEvents:
    """Socket.IO event names."""

    TASK_RUN = "task:run"
    TASK_FINISHED = "task:finished"
    QUEUE_UPDATE = "queue:update"
    INVALIDATE = "invalidate"
    INVALIDATE_SCHEMA = "invalidate:schema"
    REGISTER_EXTENSION = "register:extension"
    LEN_FRAMES = "len_frames"
    FRAME_UPDATE = "frame_update"
    INVALIDATE_GEOMETRY = "invalidate:geometry"
    INVALIDATE_FIGURE = "invalidate:figure"
    INVALIDATE_SELECTION = "invalidate:selection"
    INVALIDATE_SELECTION_GROUPS = "invalidate:selection_groups"
    INVALIDATE_BOOKMARK = "bookmarks:invalidate"
    FILESYSTEMS_UPDATE = "filesystems:update"


class LockConfig:
    """Lock configuration constants."""

    DEFAULT_TTL = 60  # Default lock TTL in seconds
    MAX_TTL = 300  # Maximum allowed lock TTL (5 minutes)
    DEFAULT_REFRESH_INTERVAL = 30  # Default refresh interval in seconds
