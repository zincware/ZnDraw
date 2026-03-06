"""Frame REST API endpoints for room frame CRUD operations.

Returns MessagePack-encoded binary data for efficient transfer:
- Success (200 OK): application/x-msgpack with raw frame data
- Error (4xx/5xx): application/problem+json per RFC 9457

Frame data format: list[dict[bytes, bytes]] where:
- Keys: b"arrays.positions", b"cell", etc.
- Values: msgpack-numpy encoded bytes
"""

import asyncio
from typing import Annotated, NoReturn

import msgpack
from fastapi import APIRouter, Query, Request, Response, status
from sqlalchemy import select as sa_select
from sqlalchemy.ext.asyncio import AsyncSession
from zndraw_joblib.dependencies import ResultBackend, request_hash
from zndraw_joblib.events import Emission, ProviderRequest, emit as joblib_emit
from zndraw_joblib.exceptions import ProviderTimeout
from zndraw_joblib.models import ProviderRecord
from zndraw_socketio import AsyncServerWrapper

from zndraw.dependencies import (
    CurrentUserFactoryDep,
    JobLibSettingsDep,
    ResultBackendDep,
    SessionDep,
    SessionMakerDep,
    SioDep,
    StorageDep,
    WritableRoomDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    FrameNotFound,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    UnprocessableContent,
    problem_responses,
)
from zndraw.redis import RedisKey
from zndraw.responses import MessagePackResponse
from zndraw.routes.rooms import build_room_update
from zndraw.schemas import (
    FrameBulkResponse,
    FrameCreateRequest,
    FrameMergeResponse,
    FrameMetadataResponse,
    FrameResponse,
    FrameUpdateRequest,
    PropertyMeta,
    StatusResponse,
)
from zndraw.socket_events import FramesInvalidate
from zndraw.storage.base import RawFrame, to_raw_frame

router = APIRouter(prefix="/v1/rooms/{room_id}/frames", tags=["frames"])

_REQUIRED_FRAME_KEYS = frozenset({b"arrays.colors", b"arrays.radii"})


def _validate_frame_keys(frame: RawFrame) -> None:
    """Raise 422 if required rendering keys are missing from a frame."""
    missing = _REQUIRED_FRAME_KEYS - frame.keys()
    if missing:
        names = ", ".join(sorted(k.decode() for k in missing))
        raise UnprocessableContent.exception(f"Frame missing required keys: {names}")


def _raise_frame_not_found(index: int | list[int], total: int) -> NoReturn:
    """Raise FrameNotFound with a consistent message."""
    range_str = f"0-{total - 1}" if total > 0 else "no frames"
    label = (
        f"Frame index {index}" if isinstance(index, int) else f"Frame indices {index}"
    )
    raise FrameNotFound.exception(f"{label} out of range ({range_str})")


# =============================================================================
# Provider dispatch helpers
# =============================================================================


async def _find_frames_provider(
    session: AsyncSession, room_id: str
) -> ProviderRecord | None:
    """Find a frames provider registered for this room."""
    result = await session.execute(
        sa_select(ProviderRecord).where(
            ProviderRecord.room_id == room_id,
            ProviderRecord.category == "frames",
        )
    )
    return result.scalar_one_or_none()


async def _dispatch_provider_frame(
    result_backend: ResultBackend,
    sio: AsyncServerWrapper,
    provider: ProviderRecord,
    index: int,
    timeout: float = 5.0,
    inflight_ttl: int = 30,
) -> RawFrame:
    """Check provider cache, dispatch if needed, and long-poll for the result.

    Returns
    -------
    RawFrame
        The unpacked frame.

    Raises
    ------
    ProblemException
        ProviderTimeout (504) with Retry-After when the provider does not
        deliver within *timeout* seconds.
    """
    params = {"index": str(index)}
    rhash = request_hash(params)
    cache_key = RedisKey.provider_result(provider.full_name, rhash)

    # Fast path: cache hit
    cached = await result_backend.get(cache_key)
    if cached is not None:
        return msgpack.unpackb(cached, raw=True)

    # Dispatch if not already inflight
    inflight_key = RedisKey.provider_inflight(provider.full_name, rhash)
    acquired = await result_backend.acquire_inflight(inflight_key, inflight_ttl)
    if acquired:
        provider_room = f"providers:{provider.full_name}"
        await joblib_emit(
            sio,
            {
                Emission(
                    ProviderRequest.from_dict_params(
                        request_id=rhash,
                        provider_name=provider.full_name,
                        params=params,
                    ),
                    provider_room,
                )
            },
        )

    # Long-poll: wait for provider to deliver
    result = await result_backend.wait_for_key(cache_key, timeout)
    if result is not None:
        return msgpack.unpackb(result, raw=True)

    raise ProviderTimeout.exception(
        detail=f"Provider '{provider.full_name}' did not respond within {timeout}s",
        headers={"Retry-After": "2"},
    )


# =============================================================================
# Frame CRUD Operations
# =============================================================================


def _filter_frames_by_keys(
    frames: list[RawFrame], key_bytes: frozenset[bytes]
) -> list[RawFrame]:
    """Filter frames to only include specified keys.

    Parameters
    ----------
    frames
        List of raw frame data (dict[bytes, bytes]).
    key_bytes
        Pre-computed set of byte keys to retain.
    """
    return [{k: v for k, v in f.items() if k in key_bytes} for f in frames]


@router.get(
    "",
    response_class=MessagePackResponse,
    responses=problem_responses(NotAuthenticated, RoomNotFound, FrameNotFound),
)
async def list_frames(
    session_maker: SessionMakerDep,
    storage: StorageDep,
    current_user: CurrentUserFactoryDep,
    sio: SioDep,
    result_backend: ResultBackendDep,
    joblib_settings: JobLibSettingsDep,
    room_id: str,
    start: Annotated[int, Query(ge=0, description="Start index (inclusive)")] = 0,
    stop: Annotated[
        int | None, Query(ge=0, description="Stop index (exclusive)")
    ] = None,
    indices: Annotated[
        str | None, Query(description="Comma-separated frame indices")
    ] = None,
    keys: Annotated[
        str | None, Query(description="Comma-separated frame keys to include")
    ] = None,
) -> Response:
    """List frames with optional range or specific indices.

    Returns MessagePack-encoded list[dict[bytes, bytes]] for efficient binary transfer.
    Each frame has bytes keys (e.g., b"arrays.positions") and msgpack-numpy encoded values.

    On storage miss with a registered frames provider, dispatches a provider
    request and returns 504 with Retry-After if the result is not yet cached.

    - Use start/stop for range queries
    - Use indices for specific frame indices (comma-separated)
    - Use keys to filter which keys are included in each frame
    """
    async with session_maker() as session:
        await verify_room(session, room_id)
        total = await storage.get_length(room_id)

        requested_indices: list[int]
        if indices:
            requested_indices = [int(i) for i in indices.split(",")]
            invalid = [i for i in requested_indices if i < 0 or i >= total]
            if invalid:
                _raise_frame_not_found(invalid, total)
            frames_or_none = await storage.get_many(room_id, requested_indices)
        else:
            effective_stop = stop if stop is not None else total
            effective_start = min(start, total)
            effective_stop = min(effective_stop, total)
            requested_indices = list(range(effective_start, effective_stop))
            frames_or_none = await storage.get_range(
                room_id, effective_start, effective_stop
            )

        has_missing = any(f is None for f in frames_or_none)
        provider = (
            await _find_frames_provider(session, room_id) if has_missing else None
        )
    # Session closed, lock released ^

    # Resolve missing frames — dispatch concurrently for O(1×timeout)
    frames: list[RawFrame]
    if has_missing:
        if provider is None:
            missing = [
                idx for idx, f in zip(requested_indices, frames_or_none) if f is None
            ]
            _raise_frame_not_found(missing, total)
        # Dispatch all missing frames concurrently (all-or-nothing:
        # if any dispatch times out, gather propagates the 504 so the
        # client retries the entire batch via Retry-After).
        missing_tasks: dict[int, asyncio.Task[RawFrame]] = {}
        async with asyncio.TaskGroup() as tg:
            for i, (idx, f) in enumerate(zip(requested_indices, frames_or_none)):
                if f is None:
                    missing_tasks[i] = tg.create_task(
                        _dispatch_provider_frame(
                            result_backend,
                            sio,
                            provider,
                            idx,
                            timeout=joblib_settings.provider_long_poll_default_seconds,
                            inflight_ttl=joblib_settings.provider_inflight_ttl_seconds,
                        )
                    )
        frames = [
            f if f is not None else missing_tasks[i].result()
            for i, f in enumerate(frames_or_none)
        ]
    else:
        frames = [f for f in frames_or_none if f is not None]

    # Filter keys if specified
    if keys:
        key_bytes = frozenset(k.strip().encode() for k in keys.split(","))
        frames = _filter_frames_by_keys(frames, key_bytes)

    return MessagePackResponse(content=frames)


@router.get(
    "/{index}",
    response_class=MessagePackResponse,
    responses=problem_responses(NotAuthenticated, RoomNotFound, FrameNotFound),
)
async def get_frame(
    session_maker: SessionMakerDep,
    storage: StorageDep,
    current_user: CurrentUserFactoryDep,
    sio: SioDep,
    result_backend: ResultBackendDep,
    joblib_settings: JobLibSettingsDep,
    room_id: str,
    index: int,
) -> Response:
    """Get a single frame by index.

    Returns MessagePack-encoded list containing the single frame (for consistency).
    Format: [dict[bytes, bytes]] - same as list_frames but with one element.

    On storage miss with a registered frames provider, dispatches a provider
    request and returns 504 with Retry-After if the result is not yet cached.
    """
    async with session_maker() as session:
        await verify_room(session, room_id)
        total = await storage.get_length(room_id)
        if index < 0 or index >= total:
            _raise_frame_not_found(index, total)

        frame = await storage.get(room_id, index)
        if frame is not None:
            return MessagePackResponse(content=[frame])

        provider = await _find_frames_provider(session, room_id)
    # Session closed, lock released ^

    if provider is None:
        _raise_frame_not_found(index, total)

    frame = await _dispatch_provider_frame(
        result_backend,
        sio,
        provider,
        index,
        timeout=joblib_settings.provider_long_poll_default_seconds,
        inflight_ttl=joblib_settings.provider_inflight_ttl_seconds,
    )
    return MessagePackResponse(content=[frame])


def _extract_property_meta(value_bytes: bytes) -> PropertyMeta:
    """Extract dtype/shape metadata from a msgpack-numpy encoded value.

    Parameters
    ----------
    value_bytes
        Raw msgpack-numpy encoded bytes for a single frame property.
    """
    decoded = msgpack.unpackb(value_bytes, raw=True)

    # msgpack-numpy arrays have an "nd" marker
    if isinstance(decoded, dict) and b"nd" in decoded:
        dtype_raw = decoded.get(b"type", b"unknown")
        dtype = dtype_raw.decode() if isinstance(dtype_raw, bytes) else str(dtype_raw)
        shape_raw = decoded.get(b"shape", ())
        shape = list(shape_raw) if isinstance(shape_raw, (list, tuple)) else []
        return PropertyMeta(dtype=dtype, shape=shape, type="array")

    # Scalar types (bool before int — bool is a subclass of int in Python)
    if isinstance(decoded, bool):
        dtype = "bool"
    elif isinstance(decoded, float):
        dtype = "float64"
    elif isinstance(decoded, int):
        dtype = "int64"
    elif isinstance(decoded, (str, bytes)):
        dtype = "str"
    else:
        dtype = "unknown"
    return PropertyMeta(dtype=dtype, type="scalar")


@router.get(
    "/{index}/metadata",
    responses=problem_responses(NotAuthenticated, RoomNotFound, FrameNotFound),
)
async def get_frame_metadata(
    session_maker: SessionMakerDep,
    storage: StorageDep,
    current_user: CurrentUserFactoryDep,
    sio: SioDep,
    result_backend: ResultBackendDep,
    joblib_settings: JobLibSettingsDep,
    room_id: str,
    index: int,
) -> FrameMetadataResponse:
    """Get metadata (keys with dtype/shape) for a specific frame.

    Returns a mapping of key names to their dtype and shape,
    without transferring the actual array data.

    On storage miss with a registered frames provider, dispatches a provider
    request and returns 504 with Retry-After if the result is not yet cached.
    """
    async with session_maker() as session:
        await verify_room(session, room_id)
        total = await storage.get_length(room_id)
        if index < 0 or index >= total:
            _raise_frame_not_found(index, total)

        frame = await storage.get(room_id, index)
        provider = (
            await _find_frames_provider(session, room_id) if frame is None else None
        )
    # Session closed, lock released ^

    if frame is None:
        if provider is None:
            _raise_frame_not_found(index, total)
        frame = await _dispatch_provider_frame(
            result_backend,
            sio,
            provider,
            index,
            timeout=joblib_settings.provider_long_poll_default_seconds,
            inflight_ttl=joblib_settings.provider_inflight_ttl_seconds,
        )

    metadata: dict[str, PropertyMeta] = {}
    for key_bytes, value_bytes in frame.items():
        key = key_bytes.decode()
        metadata[key] = _extract_property_meta(value_bytes)

    return FrameMetadataResponse(frame_id=index, metadata=metadata, source_room=room_id)


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, UnprocessableContent
    ),
)
async def append_frames(
    session: SessionDep,
    storage: StorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    request: FrameCreateRequest,
) -> FrameBulkResponse:
    """Append frames to the room's frame storage.

    Returns the updated frame count and range of appended frames.
    Broadcasts FramesInvalidateBroadcast(action="add") to the room.
    """

    raw_frames = [to_raw_frame(f) for f in request.frames]
    for frame in raw_frames:
        _validate_frame_keys(frame)

    old_total = await storage.get_length(room_id)
    new_total = await storage.extend(room_id, raw_frames)

    # Broadcast invalidation to room with new total frame count
    await sio.emit(
        FramesInvalidate(room_id=room_id, action="add", count=new_total),
        room=room_channel(room_id),
    )
    # Notify @overview of updated frame count
    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")

    return FrameBulkResponse(
        frames=request.frames,
        total=new_total,
        start=old_total,
        stop=new_total,
    )


@router.put(
    "/{index}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, FrameNotFound, RoomLocked, UnprocessableContent
    ),
)
async def update_frame(
    storage: StorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    index: int,
    request: FrameUpdateRequest,
) -> FrameResponse:
    """Update a single frame at the specified index.

    Broadcasts FramesInvalidate(action="modify") to the room.
    """

    total = await storage.get_length(room_id)
    if index < 0 or index >= total:
        _raise_frame_not_found(index, total)

    raw_frame = to_raw_frame(request.data)
    _validate_frame_keys(raw_frame)
    await storage.set_item(room_id, index, raw_frame)

    await sio.emit(
        FramesInvalidate(room_id=room_id, action="modify", indices=[index]),
        room=room_channel(room_id),
    )

    return FrameResponse(index=index, data=request.data)


@router.patch(
    "/{index}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, FrameNotFound, RoomLocked
    ),
)
async def merge_frame(
    storage: StorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    index: int,
    request: Request,
) -> FrameMergeResponse:
    """Partially update a frame by merging new data into the existing frame.

    Accepts a msgpack-encoded body with only the keys to update.
    Existing keys not present in the body are preserved.
    """
    total = await storage.get_length(room_id)
    if index < 0 or index >= total:
        _raise_frame_not_found(index, total)

    body = await request.body()
    # Use raw=False to preserve msgpack str/bin distinction.
    # The frontend's msgpack-numpy format uses str type for keys and dtype
    # (e.g., "type": "<f4"). With raw=True these become bytes, get re-packed
    # as msgpack bin, and the frontend decoder fails (checks typeof === "string").
    unpacked = msgpack.unpackb(body, raw=False)

    # Re-pack each value to match storage format (dict[bytes, bytes]).
    # Keys need bytes encoding; values keep their str/bin types intact.
    partial: dict[bytes, bytes] = {}
    for k, v in unpacked.items():
        key = k.encode() if isinstance(k, str) else k
        packed = msgpack.packb(v, use_bin_type=True)
        partial[key] = packed  # type: ignore[assignment]  # packb returns bytes, stubs say bytes | None

    await storage.merge_item(room_id, index, partial)

    await sio.emit(
        FramesInvalidate(room_id=room_id, action="modify", indices=[index]),
        room=room_channel(room_id),
    )

    updated_keys = [k.decode() for k in partial]
    return FrameMergeResponse(index=index, updated_keys=updated_keys)


@router.delete(
    "/{index}",
    status_code=status.HTTP_200_OK,
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, FrameNotFound, RoomLocked
    ),
)
async def delete_frame(
    session: SessionDep,
    storage: StorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    index: int,
) -> StatusResponse:
    """Delete a single frame at the specified index.

    Frames after the deleted index are shifted to fill the gap.
    Broadcasts FramesInvalidateBroadcast(action="delete") to the room.
    """

    total = await storage.get_length(room_id)
    if index < 0 or index >= total:
        _raise_frame_not_found(index, total)

    await storage.delete_range(room_id, index, index + 1)

    # Get new total after deletion for broadcast
    new_total = await storage.get_length(room_id)

    # Broadcast invalidation to room with new frame count
    await sio.emit(
        FramesInvalidate(
            room_id=room_id, action="delete", indices=[index], count=new_total
        ),
        room=room_channel(room_id),
    )
    # Notify @overview of updated frame count
    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")

    return StatusResponse()
