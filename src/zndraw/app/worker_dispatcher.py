"""Centralized worker queue dispatching logic.

This module consolidates queue dispatching functionality that was previously
duplicated in both events.py and routes.py.
"""

from typing import Any
import json
import logging

from .redis_keys import ExtensionKeys
from .constants import SocketEvents
from .queue_manager import emit_queue_update

log = logging.getLogger(__name__)


def dispatch_next_task(
    redis_client: Any,
    socketio_instance: Any,
    worker_id: str,
    room_id: str,
    category: str
) -> bool:
    """Check queues for all extensions a worker can handle and dispatch next task.

    Centralized function used by both:
    - Worker registration (events.py)
    - Job completion (routes.py)

    Args:
        redis_client: Redis client instance
        socketio_instance: SocketIO instance for emitting events
        worker_id: Worker client_id (UUID)
        room_id: Room identifier
        category: Extension category (e.g., 'modifiers', 'selections')

    Returns:
        bool: True if a task was dispatched, False otherwise
    """
    user_extensions_key = ExtensionKeys.user_extensions_key(room_id, category, worker_id)
    registered_extensions = redis_client.smembers(user_extensions_key)

    if not registered_extensions:
        log.debug(f"Worker {worker_id} has no registered extensions in category '{category}'")
        return False

    # Get the Socket.IO SID for this worker's client_id (needed for emit)
    worker_sid = redis_client.get(f"client_id:{worker_id}:sid")
    if not worker_sid:
        log.warning(f"Cannot dispatch: no Socket.IO SID found for worker {worker_id}")
        return False

    # Check each extension's queue
    for extension_name in registered_extensions:
        keys = ExtensionKeys.for_extension(room_id, category, extension_name)

        # Check if worker is still idle (might have been assigned already)
        is_idle = redis_client.sismember(keys.idle_workers, worker_id)
        if not is_idle:
            continue

        # Try to pop a task from queue (LPOP for FIFO)
        queued_task_json = redis_client.lpop(keys.queue)

        if queued_task_json:
            try:
                queued_task = json.loads(queued_task_json)
                task_data = queued_task.get("data")
                task_room = queued_task.get("room", room_id)
                job_id = queued_task.get("jobId")

                log.info(
                    f"Dispatching queued task (job {job_id}) for extension '{extension_name}' "
                    f"to worker {worker_id} (sid: {worker_sid}) in room '{room_id}'"
                )

                # Move worker from idle to progressing atomically
                moved = redis_client.smove(keys.idle_workers, keys.progressing_workers, worker_id)

                if moved:
                    # Emit task to worker using Socket.IO SID
                    socketio_instance.emit(
                        SocketEvents.TASK_RUN,
                        {
                            "data": task_data,
                            "extension": extension_name,
                            "category": category,
                            "room": task_room,
                            "jobId": job_id,
                        },
                        to=worker_sid,
                    )

                    # Notify all clients in room about queue update
                    emit_queue_update(redis_client, room_id, category, extension_name, socketio_instance)

                    log.info(
                        f"Successfully dispatched task from queue to worker {worker_id} "
                        f"for extension '{extension_name}'"
                    )

                    # Worker is now busy, stop checking other queues
                    return True
                else:
                    # Worker was already moved to progressing by another thread
                    # Put the task back in the queue
                    redis_client.lpush(keys.queue, queued_task_json)
                    log.warning(
                        f"Worker {worker_id} was already in progressing state, "
                        f"re-queued task for extension '{extension_name}'"
                    )
                    return False

            except json.JSONDecodeError as e:
                log.error(f"Failed to parse queued task: {e}")
                continue
            except Exception as e:
                log.error(f"Error dispatching task to worker {worker_id}: {e}", exc_info=True)
                continue

    log.debug(f"No queued tasks found for worker {worker_id} in category '{category}'")
    return False
