"""Queue management and event emission for ZnDraw."""

import logging
from typing import Any

from .constants import SocketEvents
from .redis_keys import ExtensionKeys
from .worker_stats import WorkerStats

log = logging.getLogger(__name__)


def emit_queue_update(
    redis_client: Any,
    room_id: str,
    category: str,
    extension: str,
    socketio_instance: Any,
) -> None:
    """Emit queue update event with current statistics.

    Centralizes the logic for emitting queue updates to avoid duplication.

    Args:
        redis_client: Redis client instance
        room_id: The room identifier
        category: The extension category
        extension: The extension name
        socketio_instance: The SocketIO instance to emit through
    """
    keys = ExtensionKeys.for_extension(room_id, category, extension)
    stats = WorkerStats.fetch(redis_client, keys)

    socketio_instance.emit(
        SocketEvents.QUEUE_UPDATE,
        {
            "roomId": room_id,
            "category": category,
            "extension": extension,
            **stats.to_dict(),
        },
        to=f"room:{room_id}",
    )

    log.debug(
        f"Queue update for {category}/{extension}: {stats.queue_length} queued, "
        f"{stats.idle_count} idle, {stats.progressing_count} progressing"
    )
