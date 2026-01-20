"""Utilities for cleaning up stale workers after server restart.

Stale workers are those that were registered with a previous server instance
that has since restarted. They are detected by comparing their registration
timestamp with the cluster heartbeat timestamp.
"""

import logging
import typing as t

from .cluster_heartbeat import ClusterHeartbeat
from .constants import SocketEvents
from .redis_keys import ExtensionKeys, FilesystemKeys, WorkerKeys

log = logging.getLogger(__name__)


def cleanup_stale_extension_workers(
    redis_client: t.Any,
    socketio: t.Any,
    room_id: str | None,
    category: str,
    extension: str,
) -> list[str]:
    """Clean up stale workers for a specific extension.

    Parameters
    ----------
    redis_client
        Redis client
    socketio
        SocketIO instance for emitting events
    room_id : str | None
        Room ID for room-scoped extensions, None for global
    category : str
        Extension category
    extension : str
        Extension name

    Returns
    -------
    list[str]
        List of cleaned up worker IDs
    """
    if room_id is None:
        keys = ExtensionKeys.for_global_extension(category, extension)
    else:
        keys = ExtensionKeys.for_extension(room_id, category, extension)

    # Get all workers with their registration timestamps
    workers = redis_client.hgetall(keys.workers)
    if not workers:
        return []

    stale_workers = []
    for worker_id, reg_timestamp in workers.items():
        if ClusterHeartbeat.is_worker_stale(redis_client, reg_timestamp):
            stale_workers.append(worker_id)

    if not stale_workers:
        return []

    log.info(
        f"Cleaning up {len(stale_workers)} stale worker(s) for "
        f"{category}/{extension} (room={room_id}): {stale_workers}"
    )

    # Clean up each stale worker
    for worker_id in stale_workers:
        _cleanup_single_extension_worker(
            redis_client, worker_id, room_id, category, extension, keys
        )

    # After cleanup, check if extension should be deleted
    remaining_workers = redis_client.hlen(keys.workers)
    if remaining_workers == 0:
        pending_jobs = redis_client.zcard(keys.pending_jobs)
        if pending_jobs == 0:
            # Delete the extension schema
            redis_client.hdel(keys.schema, extension)
            log.debug(
                f"Deleted orphaned extension schema: {category}/{extension} "
                f"(room={room_id})"
            )

            # Emit schema invalidation
            if room_id is None:
                socketio.emit(SocketEvents.INVALIDATE_SCHEMA, {"category": category})
            else:
                socketio.emit(
                    SocketEvents.INVALIDATE_SCHEMA,
                    {"category": category},
                    to=f"room:{room_id}",
                )

    return stale_workers


def _cleanup_single_extension_worker(
    redis_client: t.Any,
    worker_id: str,
    room_id: str | None,
    category: str,
    extension: str,
    keys: ExtensionKeys,
) -> None:
    """Clean up a single stale extension worker.

    Removes the worker from the extension's worker registry, cleans up
    reverse lookups, capacity keys, and active jobs tracking.

    Parameters
    ----------
    redis_client
        Redis client instance
    worker_id : str
        Socket ID of the stale worker
    room_id : str | None
        Room ID for room-scoped extensions, None for global
    category : str
        Extension category (modifiers, selections, etc.)
    extension : str
        Extension name
    keys : ExtensionKeys
        Pre-computed Redis keys for this extension
    """
    # Remove from workers hash
    redis_client.hdel(keys.workers, worker_id)

    # Remove from reverse lookup
    if room_id is None:
        worker_extensions_key = ExtensionKeys.global_user_extensions_key(
            category, worker_id
        )
    else:
        worker_extensions_key = ExtensionKeys.user_extensions_key(
            room_id, category, worker_id
        )
    redis_client.srem(worker_extensions_key, extension)

    # Clean up worker capacity
    capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
    redis_client.delete(capacity_key)

    # Clean up active jobs
    worker_keys = WorkerKeys(worker_id)
    redis_client.delete(worker_keys.active_jobs())

    log.debug(
        f"Cleaned up stale extension worker {worker_id} for "
        f"{category}/{extension} (room={room_id})"
    )


def cleanup_stale_filesystem_worker(
    redis_client: t.Any,
    socketio: t.Any,
    room_id: str | None,
    fs_name: str,
) -> str | None:
    """Clean up stale worker for a filesystem.

    Parameters
    ----------
    redis_client
        Redis client
    socketio
        SocketIO instance for emitting events
    room_id : str | None
        Room ID for room-scoped filesystems, None for global
    fs_name : str
        Filesystem name

    Returns
    -------
    str | None
        The cleaned up worker ID, or None if not stale.
    """
    if room_id is None:
        keys = FilesystemKeys.for_global_filesystem(fs_name)
    else:
        keys = FilesystemKeys.for_filesystem(room_id, fs_name)

    worker_id = redis_client.get(keys.worker)
    if not worker_id:
        return None

    # Get filesystem metadata to check registration time
    metadata = redis_client.hgetall(keys.metadata)
    if not metadata:
        return None

    # Check registration timestamp from metadata
    reg_timestamp = metadata.get("registration_timestamp")
    if reg_timestamp is None:
        # Legacy filesystem without timestamp - can't determine staleness
        log.debug(
            f"Filesystem {fs_name} has no registration_timestamp, skipping stale check"
        )
        return None

    if not ClusterHeartbeat.is_worker_stale(redis_client, reg_timestamp):
        return None

    # Worker is stale - clean up
    log.info(
        f"Cleaning up stale filesystem worker {worker_id} for "
        f"{fs_name} (room={room_id})"
    )

    # Clean up filesystem data
    redis_client.delete(keys.metadata)
    redis_client.delete(keys.worker)

    # Clean up reverse lookup
    if room_id is None:
        worker_filesystems_key = FilesystemKeys.global_user_filesystems_key(worker_id)
    else:
        worker_filesystems_key = FilesystemKeys.user_filesystems_key(room_id, worker_id)
    redis_client.srem(worker_filesystems_key, fs_name)

    # Emit update event
    if room_id is None:
        socketio.emit(SocketEvents.FILESYSTEMS_UPDATE, {"scope": "global"})
    else:
        socketio.emit(
            SocketEvents.FILESYSTEMS_UPDATE,
            {"scope": "room"},
            to=f"room:{room_id}",
        )

    return worker_id


def cleanup_all_stale_workers_for_category(
    redis_client: t.Any,
    socketio: t.Any,
    room_id: str | None,
    category: str,
) -> dict[str, list[str]]:
    """Clean up all stale workers for a category.

    Parameters
    ----------
    redis_client
        Redis client
    socketio
        SocketIO instance for emitting events
    room_id : str | None
        Room ID for room-scoped extensions, None for global
    category : str
        Extension category

    Returns
    -------
    dict[str, list[str]]
        Dictionary mapping extension names to lists of cleaned worker IDs
    """
    if room_id is None:
        schema_key = ExtensionKeys.global_schema_key(category)
    else:
        schema_key = ExtensionKeys.schema_key(room_id, category)

    # Get all extensions in this category
    extensions = redis_client.hkeys(schema_key)

    cleaned: dict[str, list[str]] = {}
    for ext_name in extensions:
        stale = cleanup_stale_extension_workers(
            redis_client, socketio, room_id, category, ext_name
        )
        if stale:
            cleaned[ext_name] = stale

    return cleaned
