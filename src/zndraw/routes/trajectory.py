"""Trajectory REST API endpoints for file download and upload.

Provides browser-downloadable trajectory files with format selection,
drag-and-drop upload from the frontend, and temporary download tokens
for browser/wget downloads without Authorization headers.
"""

import io
import re
import uuid
from typing import Annotated

import ase.io
from asebytes import decode, encode
from fastapi import APIRouter, Query, Request, UploadFile, status
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field

from zndraw.connectivity import add_connectivity
from zndraw.dependencies import (
    CurrentUserDep,
    OptionalUserDep,
    RedisDep,
    SessionDep,
    SioDep,
    StorageDep,
    WritableRoomDep,
    room_channel,
    verify_room,
)
from zndraw.enrichment import add_colors, add_radii
from zndraw.exceptions import (
    InvalidPayload,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    problem_responses,
)
from zndraw.redis import RedisKey
from zndraw.routes.rooms import build_room_update
from zndraw.schemas import FrameBulkResponse
from zndraw.socket_events import FramesInvalidate
from zndraw.storage.router import StorageRouter

router = APIRouter(prefix="/v1/rooms/{room_id}/trajectory", tags=["trajectory"])

# Supported formats and their MIME types / file extensions
# MIME type and file extension per format.
# extxyz and xyz share "chemical/x-xyz" — no registered MIME for extxyz.
_FORMAT_INFO: dict[str, tuple[str, str]] = {
    "extxyz": ("chemical/x-xyz", "extxyz"),
    "xyz": ("chemical/x-xyz", "xyz"),
    "cif": ("chemical/x-cif", "cif"),
    "pdb": ("chemical/x-pdb", "pdb"),
}

_MAX_UPLOAD_SIZE = 100 * 1024 * 1024  # 100 MB
_DEFAULT_TOKEN_TTL = 300  # 5 minutes
_MAX_TOKEN_TTL = 3600  # 1 hour


class DownloadTokenRequest(BaseModel):
    """Request body for creating a download token."""

    ttl: int = Field(
        default=_DEFAULT_TOKEN_TTL,
        ge=1,
        le=_MAX_TOKEN_TTL,
        description="Token TTL in seconds",
    )


class DownloadTokenResponse(BaseModel):
    """Response body for a created download token."""

    token: str
    url: str
    expires_in: int


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound, InvalidPayload),
)
async def download_trajectory(
    session: SessionDep,
    storage: StorageDep,
    redis: RedisDep,
    user: OptionalUserDep,
    room_id: str,
    format: Annotated[str, Query(description="Output format")] = "extxyz",  # noqa: A002
    indices: Annotated[
        str | None, Query(description="Comma-separated frame indices")
    ] = None,
    selection: Annotated[
        str | None, Query(description="Comma-separated atom indices to filter")
    ] = None,
    filename: Annotated[
        str | None, Query(description="Custom filename for download")
    ] = None,
    token: Annotated[
        str | None, Query(description="Download token (alternative to JWT auth)")
    ] = None,
) -> StreamingResponse:
    """Download room frames as a trajectory file.

    Authenticate via JWT header or a temporary download token.
    Supported formats: extxyz, xyz, cif, pdb.
    """
    if user is None:
        if token is None:
            raise NotAuthenticated.exception()
        stored_room = await redis.getdel(RedisKey.download_token(token))
        if stored_room is None or stored_room != room_id:
            raise NotAuthenticated.exception()

    await verify_room(session, room_id)

    if format not in _FORMAT_INFO:
        raise InvalidPayload.exception(
            f"Unsupported format '{format}'. "
            f"Use one of: {', '.join(sorted(_FORMAT_INFO))}"
        )

    total = await storage.get_length(room_id)
    if total == 0:
        raise InvalidPayload.exception("Room has no frames to download")

    if isinstance(storage, StorageRouter) and await storage.has_mount(room_id):
        raise InvalidPayload.exception(
            "Trajectory download not supported for provider-backed rooms"
        )

    # Determine which indices to export
    if indices is not None:
        try:
            index_list = [int(i.strip()) for i in indices.split(",")]
        except ValueError as err:
            raise InvalidPayload.exception("Invalid indices format") from err
        invalid = [i for i in index_list if i < 0 or i >= total]
        if invalid:
            raise InvalidPayload.exception(
                f"Frame indices {invalid} out of range (0-{total - 1})"
            )
    else:
        index_list = list(range(total))

    # Parse atom selection filter
    atom_indices: list[int] | None = None
    if selection is not None:
        try:
            atom_indices = [int(i.strip()) for i in selection.split(",")]
        except ValueError as err:
            raise InvalidPayload.exception("Invalid selection format") from err

    mime_type, ext = _FORMAT_INFO[format]

    # Validate atom selection against the first frame (before streaming starts)
    if atom_indices is not None:
        first_frame = await storage.get(room_id, index_list[0])
        if first_frame is not None:
            n_atoms = len(decode(first_frame))
            invalid_atoms = [i for i in atom_indices if i < 0 or i >= n_atoms]
            if invalid_atoms:
                raise InvalidPayload.exception(
                    f"Atom indices {invalid_atoms} out of range (0-{n_atoms - 1})"
                )

    # Build filename (sanitize user input for Content-Disposition)
    if filename is None:
        filename = f"{room_id}_{len(index_list)}_frames.{ext}"
    else:
        filename = re.sub(r"[^\w\-.]", "_", filename)

    _DOWNLOAD_BATCH = 100

    async def _generate():
        for batch_start in range(0, len(index_list), _DOWNLOAD_BATCH):
            batch_indices = index_list[batch_start : batch_start + _DOWNLOAD_BATCH]
            raw_frames = await storage.get_many(room_id, batch_indices)
            for raw_frame in raw_frames:
                if raw_frame is None:
                    continue
                atoms = decode(raw_frame)
                if atom_indices is not None:
                    atoms = atoms[atom_indices]
                buf = io.StringIO()
                ase.io.write(buf, atoms, format=format)
                yield buf.getvalue().encode("utf-8")

    return StreamingResponse(
        _generate(),
        media_type=mime_type,
        headers={
            "Content-Disposition": f'attachment; filename="{filename}"',
        },
    )


@router.post(
    "/download-tokens",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def create_download_token(
    session: SessionDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    request: Request,
    body: DownloadTokenRequest | None = None,
) -> DownloadTokenResponse:
    """Create a temporary download token for this room's trajectory.

    The token can be used as a query parameter on the GET endpoint,
    enabling downloads via browser navigation or wget without auth headers.
    """
    await verify_room(session, room_id)

    ttl = body.ttl if body else _DEFAULT_TOKEN_TTL
    token_value = uuid.uuid4().hex
    await redis.set(RedisKey.download_token(token_value), room_id, ex=ttl)

    base_url = str(request.base_url).rstrip("/")
    url = f"{base_url}/v1/rooms/{room_id}/trajectory?token={token_value}"

    return DownloadTokenResponse(token=token_value, url=url, expires_in=ttl)


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, InvalidPayload
    ),
)
async def upload_trajectory(
    session: SessionDep,
    storage: StorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    _current_user: CurrentUserDep,
    file: UploadFile,
    format: Annotated[  # noqa: A002
        str | None,
        Query(description="Format hint (inferred from extension if omitted)"),
    ] = None,
) -> FrameBulkResponse:
    """Upload a trajectory file and append frames to the room.

    Accepts multipart/form-data with a trajectory file.
    Parses the file using ASE and appends frames to storage.
    """
    content = await file.read(_MAX_UPLOAD_SIZE + 1)

    if len(content) > _MAX_UPLOAD_SIZE:
        raise InvalidPayload.exception(
            f"File exceeds maximum upload size of {_MAX_UPLOAD_SIZE} bytes"
        )

    if len(content) == 0:
        raise InvalidPayload.exception("Uploaded file is empty")

    # Determine format: explicit > extension > error
    fmt = format
    if fmt is None and file.filename:
        suffix = (
            file.filename.rsplit(".", 1)[-1].lower() if "." in file.filename else None
        )
        if suffix in _FORMAT_INFO:
            fmt = suffix
    if fmt is None:
        fmt = "extxyz"  # Default fallback

    if fmt not in _FORMAT_INFO:
        raise InvalidPayload.exception(
            f"Unsupported format '{fmt}'. Use one of: {', '.join(sorted(_FORMAT_INFO))}"
        )

    # Parse trajectory (ASE text-based formats need StringIO)
    buf = io.StringIO(content.decode("utf-8"))
    try:
        atoms_list = ase.io.read(buf, index=":", format=fmt)
    except Exception as exc:
        raise InvalidPayload.exception(
            f"Failed to parse trajectory file: {exc}"
        ) from exc

    if not isinstance(atoms_list, list):
        atoms_list = [atoms_list]

    if len(atoms_list) == 0:
        raise InvalidPayload.exception("Trajectory file contains no frames")

    # Convert to storage format and store in batches
    old_total = await storage.get_length(room_id)
    batch_size = 500
    new_total = old_total

    for i in range(0, len(atoms_list), batch_size):
        batch = atoms_list[i : i + batch_size]
        for atoms in batch:
            add_colors(atoms)
            add_radii(atoms)
            if len(atoms) < 100 and "connectivity" not in atoms.info:
                add_connectivity(atoms)
        frames = [encode(atoms) for atoms in batch]
        new_total = await storage.extend(room_id, frames)

    # Broadcast invalidation
    await sio.emit(
        FramesInvalidate(room_id=room_id, action="add", count=new_total),
        room=room_channel(room_id),
    )
    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")

    return FrameBulkResponse(
        frames=[],  # Don't echo back all frames
        total=new_total,
        start=old_total,
        stop=new_total,
    )
