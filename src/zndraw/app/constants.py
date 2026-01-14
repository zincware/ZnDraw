"""Constants for the ZnDraw application."""


class SocketEvents:
    """Socket.IO event names."""

    JOB_ASSIGNED = "job:assigned"  # Server sends job to specific worker
    JOB_STATE_CHANGED = "job:state_changed"  # Job state updates to room
    INVALIDATE = "invalidate"
    INVALIDATE_SCHEMA = "invalidate:schema"
    REGISTER_EXTENSION = "register:extension"
    LEN_FRAMES = "len_frames"
    FRAME_UPDATE = "frame_update"
    INVALIDATE_GEOMETRY = "geometry:invalidate "
    INVALIDATE_FIGURE = "invalidate:figure"
    INVALIDATE_SELECTION = "invalidate:selection"
    INVALIDATE_SELECTION_GROUPS = "invalidate:selection_groups"
    INVALIDATE_BOOKMARK = "bookmarks:invalidate"
    INVALIDATE_FRAMES = "frames:invalidate"
    FILESYSTEMS_UPDATE = "filesystems:update"
    FILESYSTEM_LIST = "filesystem:list"
    FILESYSTEM_METADATA = "filesystem:metadata"
    FILESYSTEM_LOAD = "filesystem:load"
    ROOM_DELETE = "room:delete"
    ROOM_UPDATE = "room:update"
    PROGRESS_UPDATED = "progress:updated"
    PROGRESS_COMPLETED = "progress:completed"
    CHAT_MESSAGE_NEW = "chat:message:new"
    CHAT_MESSAGE_UPDATED = "chat:message:updated"
    ACTIVE_CAMERA_UPDATE = "active_camera:update"
    SELECTION_UPDATE = "selection:update"
    FRAME_SELECTION_UPDATE = "frame_selection:update"


class LockConfig:
    """Lock configuration constants."""

    DEFAULT_TTL = 60  # Default lock TTL in seconds
    MAX_TTL = 300  # Maximum allowed lock TTL (5 minutes)
    DEFAULT_REFRESH_INTERVAL = 30  # Default refresh interval in seconds
