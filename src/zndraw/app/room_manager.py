"""Helper functions for room metadata management and WebSocket event emission."""

import json
import logging
import typing as t

from flask_socketio import SocketIO
from redis import Redis

from zndraw.app.constants import SocketEvents
from zndraw.app.models import RoomMetadata
from zndraw.app.redis_keys import RoomKeys

log = logging.getLogger(__name__)


def get_room_metadata(redis_client: t.Any, room_id: str) -> RoomMetadata:
    """Fetch complete room metadata from Redis.

    Parameters
    ----------
    redis_client : Redis
        Redis client instance
    room_id : str
        Room identifier

    Returns
    -------
    RoomMetadata
        Complete room metadata
    """
    room_keys = RoomKeys(room_id)

    # Get basic fields
    description = redis_client.get(room_keys.description())
    locked = redis_client.get(room_keys.locked()) == "1"
    hidden = redis_client.get(room_keys.hidden()) == "1"

    # Get frame count
    frame_count = redis_client.zcard(room_keys.trajectory_indices())

    # Check if default room
    default_room = redis_client.get("default_room")
    is_default = default_room == room_id

    # Get presenter (if exists) - assuming stored as presenter:{room_id}
    presenter_sid = redis_client.get(f"presenter:{room_id}")

    return RoomMetadata(
        id=room_id,
        description=description,
        frameCount=frame_count,
        locked=locked,
        hidden=hidden,
        isDefault=is_default,
        presenterSid=presenter_sid,
    )


def emit_room_update(
    socketio: SocketIO, room_id: str, skip_sid: str | None = None, **changes
):
    """Emit room:update event to both overview:public and room:<room_id>.

    This function broadcasts room metadata changes to all relevant clients:
    - Clients in overview:public (room list page)
    - Clients in room:<room_id> (specific room page)

    Only sends changed fields for efficiency.

    Parameters
    ----------
    socketio : SocketIO
        Flask-SocketIO instance
    room_id : str
        Room identifier
    skip_sid : str | None
        Optional socket ID to skip when broadcasting (typically the client that
        initiated the change, to prevent them from receiving their own update)
    **changes
        Changed fields to broadcast (e.g., frameCount=120, locked=True)
        All fields must be valid RoomMetadata attributes

    Raises
    ------
    ValueError
        If any field in changes is not a valid RoomMetadata attribute

    Examples
    --------
    >>> emit_room_update(socketio, "my-room", frameCount=120)
    >>> emit_room_update(socketio, "my-room", locked=True, description="Updated")
    >>> # Skip the initiating client
    >>> emit_room_update(socketio, "my-room", skip_sid=request.sid, presenterSid="abc123")
    """
    # Validate that all changes are valid RoomMetadata fields
    valid_fields = set(RoomMetadata.model_fields.keys())

    # Allow special fields that are event-specific but not part of the model
    special_fields = {"created"}  # Flag indicating new room creation

    # roomId is added in payload, not part of changes
    invalid_fields = set(changes.keys()) - valid_fields - special_fields

    if invalid_fields:
        raise ValueError(
            f"Invalid fields for room:update event: {invalid_fields}. "
            f"Valid RoomMetadata fields: {valid_fields}. "
            f"Special event fields: {special_fields}"
        )

    payload = {"roomId": room_id, **changes}

    # Broadcast to overview:public (room list), optionally skipping initiating client
    socketio.emit(
        SocketEvents.ROOM_UPDATE, payload, to="overview:public", skip_sid=skip_sid, namespace="/"
    )

    # Broadcast to specific room, optionally skipping initiating client
    socketio.emit(
        SocketEvents.ROOM_UPDATE, payload, to=f"room:{room_id}", skip_sid=skip_sid, namespace="/"
    )

    log.debug(f"Emitted room:update for '{room_id}' (skip_sid={skip_sid}): {changes}")


def emit_room_delete(socketio: SocketIO, room_id: str):
    """Emit room:delete event to both overview:public and room:<room_id>.

    Parameters
    ----------
    socketio : SocketIO
        Flask-SocketIO instance
    room_id : str
        Room identifier that was deleted
    """
    payload = {"roomId": room_id}

    # Broadcast to overview:public (room list)
    socketio.emit(SocketEvents.ROOM_DELETE, payload, to="overview:public", namespace="/")

    # Broadcast to specific room (clients should navigate away)
    socketio.emit(SocketEvents.ROOM_DELETE, payload, to=f"room:{room_id}", namespace="/")

    log.info(f"Emitted room:delete for '{room_id}'")
