"""Constants for the ZnDraw application."""


class SocketEvents:
    """Socket.IO event names."""

    JOB_ASSIGNED = "job:assign"  # Server sends job to specific worker
    JOB_STATE_CHANGED = "job:update"  # Job state updates to room
    INVALIDATE = "invalidate"
    INVALIDATE_SCHEMA = "schema:invalidate"
    REGISTER_EXTENSION = "register:extension"
    LEN_FRAMES = "len_frames"
    FRAME_UPDATE = "frame:update"
    INVALIDATE_GEOMETRY = "geometry:invalidate"
    INVALIDATE_FIGURE = "figure:invalidate"
    INVALIDATE_SELECTION = "selection:invalidate"
    INVALIDATE_SELECTION_GROUPS = "selection-groups:invalidate"
    INVALIDATE_BOOKMARK = "bookmarks:invalidate"
    INVALIDATE_FRAMES = "frames:invalidate"
    FILESYSTEMS_UPDATE = "filesystems:update"
    FILESYSTEM_LIST = "filesystem:list"
    FILESYSTEM_METADATA = "filesystem:metadata"
    FILESYSTEM_LOAD = "filesystem:load"
    ROOM_DELETE = "room:delete"
    ROOM_UPDATE = "room:update"
    PROGRESS_UPDATED = "progress:update"
    PROGRESS_COMPLETED = "progress:complete"
    CHAT_MESSAGE_NEW = "chat:new"
    CHAT_MESSAGE_UPDATED = "chat:update"
    ACTIVE_CAMERA_UPDATE = "active-camera:update"
    SELECTION_UPDATE = "selection:update"
    FRAME_SELECTION_UPDATE = "frame-selection:update"


class LockConfig:
    """Lock configuration constants."""

    DEFAULT_TTL = 60  # Default lock TTL in seconds
    MAX_TTL = 300  # Maximum allowed lock TTL (5 minutes)
    DEFAULT_REFRESH_INTERVAL = 30  # Default refresh interval in seconds
