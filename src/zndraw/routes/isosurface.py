"""Isosurface REST endpoint — marching cubes extraction for volumetric data."""

from __future__ import annotations

from typing import Annotated

import msgpack
import msgpack_numpy
import numpy as np
from fastapi import APIRouter, Query, Response

from zndraw.dependencies import (
    CurrentUserFactoryDep,
    JobLibSettingsDep,
    ResultBackendDep,
    SessionMakerDep,
    SioDep,
    StorageDep,
    verify_room,
)
from zndraw.exceptions import (
    FrameNotFound,
    RoomNotFound,
    UnprocessableContent,
    problem_responses,
)
from zndraw.routes.frames import _dispatch_provider_frame, _find_frames_provider

router = APIRouter(
    prefix="/v1/rooms/{room_id}/frames/{index}/isosurface",
    tags=["isosurface"],
)


def _extract_mesh(
    grid: np.ndarray,
    origin: np.ndarray,
    cell: np.ndarray,
    isovalue: float,
    step_size: int,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Run marching cubes and transform vertices to world space.

    Parameters
    ----------
    grid : np.ndarray
        3D scalar field (Nx, Ny, Nz).
    origin : np.ndarray
        World-space origin (3,).
    cell : np.ndarray
        Axis vectors spanning the grid (3, 3).
    isovalue : float
        Scalar threshold for surface extraction.
    step_size : int
        Marching cubes step size (1 = full resolution).

    Returns
    -------
    tuple[np.ndarray, np.ndarray] | None
        (vertices, faces) as float32/uint32 arrays, or None if no surface.
    """
    from skimage.measure import marching_cubes

    try:
        verts, faces, _, _ = marching_cubes(grid, level=isovalue, step_size=step_size)
    except (ValueError, RuntimeError):
        return None

    if len(verts) == 0:
        return None

    shape = np.array(grid.shape, dtype=np.float64)
    fractional = verts / (shape - 1)
    world = origin + fractional @ cell

    return world.astype(np.float32), faces.astype(np.uint32)


@router.get(
    "",
    responses={
        **problem_responses(RoomNotFound, FrameNotFound),
        **problem_responses(UnprocessableContent),
    },
)
async def get_isosurface(
    session_maker: SessionMakerDep,
    storage: StorageDep,
    sio: SioDep,
    result_backend: ResultBackendDep,
    joblib_settings: JobLibSettingsDep,
    _current_user: CurrentUserFactoryDep,
    room_id: str,
    index: int,
    cube_key: Annotated[str, Query(description="Frame key for volumetric data dict")],
    isovalue: Annotated[float, Query(description="Scalar threshold")],
    resolution: Annotated[
        float, Query(ge=0.0, le=1.0, description="Mesh resolution (0=coarse, 1=fine)")
    ] = 1.0,
    sigma: Annotated[
        float, Query(ge=0.0, description="Gaussian smoothing sigma (0=disabled)")
    ] = 0.0,
) -> Response:
    """Extract an isosurface mesh from volumetric frame data."""
    async with session_maker() as session:
        await verify_room(session, room_id)
        total = await storage.get_length(room_id)
        if index < 0 or index >= total:
            FrameNotFound.raise_out_of_range(index, total)

        frame = await storage.get(room_id, index)
        provider = (
            await _find_frames_provider(session, room_id) if frame is None else None
        )

    if frame is None:
        if provider is None:
            raise FrameNotFound.exception(f"Frame at index {index} is empty")
        frame = await _dispatch_provider_frame(
            result_backend,
            sio,
            provider,
            index,
            timeout=joblib_settings.provider_long_poll_default_seconds,
            inflight_ttl=joblib_settings.provider_inflight_ttl_seconds,
        )

    key_bytes = cube_key.encode()
    if key_bytes not in frame:
        raise UnprocessableContent.exception(
            f"Key '{cube_key}' not found in frame {index}"
        )

    cube_dict = msgpack.unpackb(
        frame[key_bytes], object_hook=msgpack_numpy.decode, raw=False
    )

    required_keys = {"grid", "origin", "cell"}
    if not isinstance(cube_dict, dict):
        raise UnprocessableContent.exception(f"Key '{cube_key}' is not a dict")
    missing = required_keys - cube_dict.keys()
    if missing:
        raise UnprocessableContent.exception(
            f"Cube data missing keys: {', '.join(sorted(missing))}"
        )

    grid = np.asarray(cube_dict["grid"])
    if grid.ndim != 3:
        raise UnprocessableContent.exception(f"Grid must be 3D, got {grid.ndim}D")

    origin = np.asarray(cube_dict["origin"], dtype=np.float64)
    if origin.shape != (3,):
        raise UnprocessableContent.exception(
            f"Origin must have shape (3,), got {origin.shape}"
        )
    cell = np.asarray(cube_dict["cell"], dtype=np.float64)
    if cell.shape != (3, 3):
        raise UnprocessableContent.exception(
            f"Cell must have shape (3, 3), got {cell.shape}"
        )

    if sigma > 0:
        from scipy.ndimage import gaussian_filter

        grid = gaussian_filter(grid, sigma=sigma)

    # resolution 1.0 → step_size 1, resolution 0.0 → step_size 8
    step_size = max(1, round(1 + (1 - resolution) * 7))

    result = _extract_mesh(grid, origin, cell, isovalue, step_size)

    if result is None:
        packed = msgpack.packb({"vertices": [], "faces": []})
    else:
        vertices, faces = result
        packed = msgpack.packb(
            {"vertices": vertices, "faces": faces},
            default=msgpack_numpy.encode,
        )

    return Response(content=packed, media_type="application/x-msgpack")
