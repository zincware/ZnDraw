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
