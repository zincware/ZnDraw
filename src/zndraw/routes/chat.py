"""Chat REST API endpoints for room messages."""

from datetime import UTC, datetime

from fastapi import APIRouter, Query, status
from sqlalchemy import func
from sqlmodel import col, select
from zndraw_auth import User

from zndraw.dependencies import (
    CurrentUserDep,
    SessionDep,
    SioDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    MessageNotFound,
    NotAuthenticated,
    NotMessageOwner,
    RoomNotFound,
    problem_responses,
)
from zndraw.models import Message
from zndraw.schemas import (
    MessageCreate,
    MessageEditRequest,
    MessageResponse,
    MessagesMetadata,
    MessagesResponse,
)
from zndraw.socket_events import MessageEdited, MessageNew

router = APIRouter(prefix="/v1/rooms/{room_id}/chat/messages", tags=["chat"])


def _datetime_to_unix_ms(dt: datetime) -> int:
    """Convert a datetime to unix milliseconds."""
    return int(dt.timestamp() * 1000)


def _message_to_response(msg: Message, email: str | None = None) -> MessageResponse:
    """Convert a Message model to MessageResponse."""
    return MessageResponse(
        id=msg.id,  # type: ignore[arg-type]
        room_id=msg.room_id,
        user_id=msg.user_id,
        content=msg.content,
        created_at=msg.created_at,
        updated_at=msg.updated_at,
        email=email,
    )


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_messages(
    session: SessionDep,
    _current_user: CurrentUserDep,
    room_id: str,
    limit: int = Query(default=30, ge=1, le=100),
    before: int | None = Query(default=None, description="Unix ms cursor"),
) -> MessagesResponse:
    """List messages with cursor pagination (newest first)."""
    await verify_room(session, room_id)

    stmt = select(Message).where(Message.room_id == room_id)

    if before is not None:
        before_dt = datetime.fromtimestamp(before / 1000, tz=UTC)
        stmt = stmt.where(col(Message.created_at) < before_dt)

    # Fetch one extra to determine has_more
    stmt = stmt.order_by(col(Message.created_at).desc()).limit(limit + 1)
    result = await session.execute(stmt)
    rows = list(result.scalars().all())

    has_more = len(rows) > limit
    rows = rows[:limit]

    # Get total count
    count_stmt = (
        select(func.count()).select_from(Message).where(Message.room_id == room_id)
    )
    total_count = (await session.execute(count_stmt)).scalar_one()

    # Look up emails for all user_ids
    user_ids = {row.user_id for row in rows}
    email_map: dict[str, str | None] = {}
    if user_ids:
        users_result = await session.execute(
            select(User).where(col(User.id).in_(user_ids))
        )
        for user in users_result.scalars().all():
            email_map[str(user.id)] = user.email

    messages = [
        _message_to_response(row, email_map.get(str(row.user_id))) for row in rows
    ]

    oldest_ts = _datetime_to_unix_ms(rows[-1].created_at) if rows else None
    newest_ts = _datetime_to_unix_ms(rows[0].created_at) if rows else None

    return MessagesResponse(
        items=messages,
        metadata=MessagesMetadata(
            has_more=has_more,
            total_count=total_count,
            oldest_timestamp=oldest_ts,
            newest_timestamp=newest_ts,
        ),
    )


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def create_message(
    session: SessionDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    room_id: str,
    request: MessageCreate,
) -> MessageResponse:
    """Create a new chat message."""
    await verify_room(session, room_id)

    msg = Message(
        room_id=room_id,
        user_id=current_user.id,  # type: ignore[arg-type]
        content=request.content,
    )
    session.add(msg)
    await session.commit()
    await session.refresh(msg)

    email = current_user.email

    await sio.emit(
        MessageNew(
            id=msg.id,  # type: ignore[arg-type]
            room_id=room_id,
            user_id=current_user.id,  # type: ignore[arg-type]
            content=msg.content,
            created_at=msg.created_at,
            email=email,
        ),
        room=room_channel(room_id),
    )

    return _message_to_response(msg, email)


@router.patch(
    "/{message_id}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, MessageNotFound, NotMessageOwner
    ),
)
async def edit_message(
    session: SessionDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    room_id: str,
    message_id: int,
    request: MessageEditRequest,
) -> MessageResponse:
    """Edit an existing chat message. Only the author can edit."""
    await verify_room(session, room_id)

    msg = await session.get(Message, message_id)
    if msg is None or msg.room_id != room_id:
        raise MessageNotFound.exception(f"Message {message_id} not found")

    if msg.user_id != current_user.id:
        raise NotMessageOwner.exception("Only the message author can edit")

    msg.content = request.content
    msg.updated_at = datetime.now(UTC)
    await session.commit()
    await session.refresh(msg)

    await sio.emit(
        MessageEdited(
            id=msg.id,  # type: ignore[arg-type]
            room_id=room_id,
            content=msg.content,
            updated_at=msg.updated_at,
        ),
        room=room_channel(room_id),
    )

    return _message_to_response(msg, current_user.email)
