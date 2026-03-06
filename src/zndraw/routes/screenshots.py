"""Screenshot REST API endpoints for room screenshots."""

import base64
from pathlib import Path

from fastapi import APIRouter, Form, Query, Response, UploadFile, status
from sqlalchemy import func
from sqlmodel import col, select

from zndraw.dependencies import (
    CurrentUserDep,
    MediaPathDep,
    RedisDep,
    SessionDep,
    SioDep,
    verify_room,
)
from zndraw.exceptions import (
    InvalidScreenshotFormat,
    NoFrontendSession,
    NotAuthenticated,
    RoomNotFound,
    ScreenshotNotFound,
    ScreenshotNotPending,
    ScreenshotTooLarge,
    problem_responses,
)
from zndraw.models import Screenshot
from zndraw.redis import RedisKey
from zndraw.schemas import (
    OffsetPage,
    ScreenshotCaptureCreate,
    ScreenshotListItem,
    ScreenshotResponse,
    StatusResponse,
)
from zndraw.socket_events import ScreenshotRequest

router = APIRouter(prefix="/v1/rooms/{room_id}/screenshots", tags=["screenshots"])

_ALLOWED_FORMATS = {"png", "jpeg", "webp"}
_MAX_FILE_SIZE = 10 * 1024 * 1024  # 10 MB


def _screenshot_dir(media_path: Path, room_id: str) -> Path:
    return media_path / room_id / "screenshots"


def _screenshot_path(media_path: Path, room_id: str, row: Screenshot) -> Path:
    return _screenshot_dir(media_path, room_id) / f"{row.id}.{row.format}"


def _validate_upload(format: str, file_bytes: bytes) -> bytes:
    """Validate screenshot format and size, returning the bytes on success."""
    if format not in _ALLOWED_FORMATS:
        raise InvalidScreenshotFormat.exception(
            f"Format '{format}' is not supported. "
            f"Use one of: {', '.join(sorted(_ALLOWED_FORMATS))}"
        )
    if len(file_bytes) > _MAX_FILE_SIZE:
        raise ScreenshotTooLarge.exception(
            f"File size {len(file_bytes)} exceeds maximum of {_MAX_FILE_SIZE} bytes"
        )
    return file_bytes


def _to_response(row: Screenshot, data: bytes | None = None) -> ScreenshotResponse:
    encoded = base64.b64encode(data).decode("ascii") if data is not None else None
    return ScreenshotResponse(
        id=row.id,  # type: ignore[arg-type]
        room_id=row.room_id,
        format=row.format,
        size=row.size,
        width=row.width,
        height=row.height,
        status=row.status,  # type: ignore[arg-type]
        created_by_id=row.created_by_id,
        created_at=row.created_at,
        data=encoded,
    )


# =============================================================================
# POST — Upload screenshot (multipart)
# =============================================================================


@router.post(
    "/upload",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, InvalidScreenshotFormat, ScreenshotTooLarge
    ),
)
async def upload_screenshot(
    session: SessionDep,
    current_user: CurrentUserDep,
    media_path: MediaPathDep,
    room_id: str,
    file: UploadFile,
    format: str = Form(default="png"),
    width: int | None = Form(default=None),
    height: int | None = Form(default=None),
) -> ScreenshotResponse:
    """Upload a screenshot file directly."""
    await verify_room(session, room_id)

    file_bytes = _validate_upload(format, await file.read())

    row = Screenshot(
        room_id=room_id,
        format=format,
        size=len(file_bytes),
        width=width,
        height=height,
        status="completed",
        created_by_id=current_user.id,  # type: ignore[arg-type]
    )
    session.add(row)
    await session.commit()
    await session.refresh(row)

    # Write file to disk
    dir_path = _screenshot_dir(media_path, room_id)
    dir_path.mkdir(parents=True, exist_ok=True)
    file_path = _screenshot_path(media_path, room_id, row)
    try:
        file_path.write_bytes(file_bytes)
    except OSError:
        await session.delete(row)
        await session.commit()
        raise

    return _to_response(row, data=file_bytes)


# =============================================================================
# POST — Request capture (JSON)
# =============================================================================


@router.post(
    "",
    status_code=status.HTTP_202_ACCEPTED,
    responses=problem_responses(NotAuthenticated, RoomNotFound, NoFrontendSession),
)
async def request_capture(
    session: SessionDep,
    current_user: CurrentUserDep,
    sio: SioDep,
    redis: RedisDep,
    room_id: str,
    request: ScreenshotCaptureCreate,
) -> Response:
    """Request a screenshot capture from a frontend session."""
    await verify_room(session, room_id)

    # Verify session_id is a live frontend session
    active_cam = await redis.hexists(  # type: ignore[misc]
        RedisKey.active_cameras(room_id), request.session_id
    )
    if not active_cam:
        raise NoFrontendSession.exception(
            f"Session '{request.session_id}' is not an active frontend session"
        )

    row = Screenshot(
        room_id=room_id,
        status="pending",
        size=0,
        created_by_id=current_user.id,  # type: ignore[arg-type]
    )
    session.add(row)
    await session.commit()
    await session.refresh(row)

    upload_url = f"/v1/rooms/{room_id}/screenshots/{row.id}"

    await sio.emit(
        ScreenshotRequest(
            screenshot_id=row.id,  # type: ignore[arg-type]
            upload_url=upload_url,
        ),
        to=request.session_id,
    )

    resp = _to_response(row)
    return Response(
        content=resp.model_dump_json(),
        status_code=status.HTTP_202_ACCEPTED,
        media_type="application/json",
        headers={"Location": upload_url},
    )


# =============================================================================
# GET — List screenshots
# =============================================================================


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_screenshots(
    session: SessionDep,
    _current_user: CurrentUserDep,
    room_id: str,
    limit: int = Query(default=20, ge=1, le=100),
    offset: int = Query(default=0, ge=0),
) -> OffsetPage[ScreenshotListItem]:
    """List completed screenshots with offset/limit pagination."""
    await verify_room(session, room_id)

    stmt = (
        select(Screenshot)
        .where(Screenshot.room_id == room_id, Screenshot.status == "completed")
        .order_by(col(Screenshot.created_at).desc())
        .offset(offset)
        .limit(limit)
    )
    result = await session.execute(stmt)
    rows = list(result.scalars().all())

    count_stmt = (
        select(func.count())
        .select_from(Screenshot)
        .where(Screenshot.room_id == room_id, Screenshot.status == "completed")
    )
    total = (await session.execute(count_stmt)).scalar_one()

    items = [
        ScreenshotListItem(
            id=r.id,  # type: ignore[arg-type]
            room_id=r.room_id,
            format=r.format,
            size=r.size,
            width=r.width,
            height=r.height,
            created_by_id=r.created_by_id,
            created_at=r.created_at,
        )
        for r in rows
    ]

    return OffsetPage(items=items, total=total, limit=limit, offset=offset)


# =============================================================================
# GET — Single screenshot
# =============================================================================


@router.get(
    "/{screenshot_id}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, ScreenshotNotFound),
)
async def get_screenshot(
    session: SessionDep,
    _current_user: CurrentUserDep,
    media_path: MediaPathDep,
    room_id: str,
    screenshot_id: int,
) -> ScreenshotResponse:
    """Get a single screenshot by ID."""
    await verify_room(session, room_id)

    row = await session.get(Screenshot, screenshot_id)
    if row is None or row.room_id != room_id:
        raise ScreenshotNotFound.exception(f"Screenshot {screenshot_id} not found")

    data: bytes | None = None
    if row.status == "completed":
        file_path = _screenshot_path(media_path, room_id, row)
        if file_path.exists():
            data = file_path.read_bytes()

    return _to_response(row, data=data)


# =============================================================================
# PATCH — Complete a pending screenshot
# =============================================================================


@router.patch(
    "/{screenshot_id}",
    responses=problem_responses(
        NotAuthenticated,
        RoomNotFound,
        ScreenshotNotFound,
        ScreenshotNotPending,
        InvalidScreenshotFormat,
        ScreenshotTooLarge,
    ),
)
async def complete_screenshot(
    session: SessionDep,
    _current_user: CurrentUserDep,
    media_path: MediaPathDep,
    room_id: str,
    screenshot_id: int,
    file: UploadFile,
    format: str = Form(default="png"),
    width: int | None = Form(default=None),
    height: int | None = Form(default=None),
) -> ScreenshotResponse:
    """Complete a pending screenshot by uploading the captured image."""
    await verify_room(session, room_id)

    row = await session.get(Screenshot, screenshot_id)
    if row is None or row.room_id != room_id:
        raise ScreenshotNotFound.exception(f"Screenshot {screenshot_id} not found")

    if row.status != "pending":
        raise ScreenshotNotPending.exception(
            f"Screenshot {screenshot_id} is already completed"
        )

    file_bytes = _validate_upload(format, await file.read())

    row.status = "completed"
    row.size = len(file_bytes)
    row.format = format
    row.width = width
    row.height = height

    dir_path = _screenshot_dir(media_path, room_id)
    dir_path.mkdir(parents=True, exist_ok=True)
    file_path = _screenshot_path(media_path, room_id, row)
    try:
        file_path.write_bytes(file_bytes)
    except OSError:
        await session.rollback()
        raise

    await session.commit()
    await session.refresh(row)

    return _to_response(row, data=file_bytes)


# =============================================================================
# DELETE — Delete a screenshot
# =============================================================================


@router.delete(
    "/{screenshot_id}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, ScreenshotNotFound),
)
async def delete_screenshot(
    session: SessionDep,
    _current_user: CurrentUserDep,
    media_path: MediaPathDep,
    room_id: str,
    screenshot_id: int,
) -> StatusResponse:
    """Delete a screenshot and its file."""
    await verify_room(session, room_id)

    row = await session.get(Screenshot, screenshot_id)
    if row is None or row.room_id != room_id:
        raise ScreenshotNotFound.exception(f"Screenshot {screenshot_id} not found")

    # Delete file from disk
    file_path = _screenshot_path(media_path, room_id, row)
    file_path.unlink(missing_ok=True)

    await session.delete(row)
    await session.commit()

    return StatusResponse()


# =============================================================================
# Cleanup helper
# =============================================================================


def delete_room_screenshots(media_path: Path, room_id: str) -> None:
    """Remove the screenshots directory for a room.

    Call this when deleting a room. DB rows should cascade via FK.
    """
    import shutil

    dir_path = _screenshot_dir(media_path, room_id)
    if dir_path.exists():
        shutil.rmtree(dir_path)
