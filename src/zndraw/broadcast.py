"""Room-scoped broadcast helpers."""

from uuid import UUID

from zndraw_socketio import AsyncServerWrapper

from zndraw.models import Room
from zndraw.socket_events import RoomScopedEvent


def room_channel(room_id: str | UUID) -> str:
    return f"room:{room_id}"


async def broadcast_to_room(
    sio: AsyncServerWrapper,
    event: RoomScopedEvent,
    room: Room,
    *,
    also_notify_user: UUID | str | None = None,
    skip_sid: str | None = None,
) -> None:
    """Emit ``event`` on the room channel; optionally fan out to one user channel."""
    assert UUID(room.id) == event.room_id, (
        f"event.room_id {event.room_id} does not match room.id {room.id}"
    )
    if skip_sid is not None:
        await sio.emit(event, room=room_channel(room.id), skip_sid=skip_sid)
    else:
        await sio.emit(event, room=room_channel(room.id))
    if also_notify_user is not None:
        await sio.emit(event, room=f"user:{also_notify_user}")
