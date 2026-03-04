"""CLI commands for GIF animation export.

Two composable axes:
- **Trajectory**: iterate over frame indices (``--traj-start/stop/step``)
- **Camera curve**: sweep camera along a curve (``--orbit`` or ``--curve KEY``)

If an axis is omitted it stays fixed at its current value throughout.
"""

from __future__ import annotations

import math
import sys
import time
import uuid
from pathlib import Path
from typing import Annotated

import typer

from .connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    cli_error_handler,
    die,
    EXIT_CLIENT_ERROR,
    get_zndraw,
    resolve_room,
)

gif_app = typer.Typer(name="gif", help="GIF animation export")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_session(vis):
    """Return the first active session or exit with an error."""
    sids = list(vis.sessions)
    if not sids:
        die(
            "No Active Sessions",
            "No browser sessions found. Open the room in a browser first.",
            400,
            EXIT_CLIENT_ERROR,
        )
    return vis.sessions[sids[0]]


def _build_schedule(
    traj_values: list[int] | None,
    curve_values: list[float] | None,
) -> list[tuple[int | None, float | None]]:
    """Zip both axes, stretching the shorter one via linear interpolation.

    When one axis has fewer values than the other, its values are evenly
    distributed across the total frame count (instead of stalling at the
    last value).

    Parameters
    ----------
    traj_values
        Frame indices or ``None`` for fixed trajectory.
    curve_values
        Curve progress values or ``None`` for fixed camera.

    Returns
    -------
    list[tuple[int | None, float | None]]
        Each entry is ``(frame_index_or_None, progress_or_None)``.
    """
    if traj_values is None and curve_values is None:
        return [(None, None)]

    import numpy as np

    n_traj = len(traj_values) if traj_values else 0
    n_curve = len(curve_values) if curve_values else 0
    total = max(n_traj, n_curve, 1)

    def _stretch(values: list, total: int) -> list:
        """Resample *values* to length *total* via nearest-index mapping."""
        indices = np.linspace(0, len(values) - 1, total).astype(int)
        return [values[i] for i in indices]

    stretched_traj = _stretch(traj_values, total) if traj_values else [None] * total
    stretched_curve = _stretch(curve_values, total) if curve_values else [None] * total
    return list(zip(stretched_traj, stretched_curve))


def _assemble_gif(frames: list[bytes], output: Path, fps: int) -> None:
    """Convert PNG byte frames into an animated GIF.

    Parameters
    ----------
    frames
        List of raw PNG bytes.
    output
        Destination file path.
    fps
        Frames per second.
    """
    try:
        from PIL import Image
    except ImportError:
        die(
            "Missing Dependency",
            "Pillow is required for GIF export. Install with: uv add 'zndraw[gif]'",
            400,
            EXIT_CLIENT_ERROR,
        )

    import io

    images = [Image.open(io.BytesIO(f)) for f in frames]
    duration_ms = max(1, round(1000 / fps))
    images[0].save(
        output,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        disposal=2,
    )


def _parse_center(center_str: str) -> tuple[float, float, float]:
    """Parse ``x,y,z`` string into a 3-tuple."""
    parts = center_str.split(",")
    if len(parts) != 3:
        msg = f"Expected x,y,z format, got: {center_str!r}"
        raise typer.BadParameter(msg)
    return (float(parts[0]), float(parts[1]), float(parts[2]))


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------


@gif_app.command("capture")
def capture(
    output: Annotated[Path, typer.Option("-o", "--output", help="Output GIF path")] = Path("output.gif"),
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    # Camera mode
    orbit: Annotated[bool, typer.Option("--orbit", help="Auto-create orbit circle")] = False,
    curve: Annotated[str | None, typer.Option("--curve", help="Existing curve geometry key")] = None,
    radius: Annotated[float, typer.Option(help="Orbit radius (with --orbit)")] = 15.0,
    center: Annotated[str | None, typer.Option(help="Orbit center as x,y,z (default: camera target)")] = None,
    # Curve axis
    curve_start: Annotated[float, typer.Option(help="Curve start progress (0.0-1.0)")] = 0.0,
    curve_stop: Annotated[float, typer.Option(help="Curve stop progress (0.0-1.0)")] = 1.0,
    curve_step: Annotated[float, typer.Option(help="Curve step size (0 = no curve)")] = 0.0,
    # Trajectory axis
    traj_start: Annotated[int | None, typer.Option(help="Trajectory start frame")] = None,
    traj_stop: Annotated[int | None, typer.Option(help="Trajectory stop frame")] = None,
    traj_step: Annotated[int, typer.Option(help="Trajectory step")] = 1,
    # Output options
    fps: Annotated[int, typer.Option(help="Frames per second")] = 20,
    delay: Annotated[float, typer.Option(help="Delay between captures (seconds)")] = 0.02,
) -> None:
    """Capture an animated GIF from a ZnDraw room.

    Composes two independent axes: trajectory playback and camera motion.
    If neither axis is specified, captures a single still frame.
    """
    with cli_error_handler():
        # Validate mutually exclusive camera modes
        if orbit and curve is not None:
            msg = "--orbit and --curve are mutually exclusive"
            raise typer.BadParameter(msg)

        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        session = _get_session(vis)

        # Build trajectory range
        has_traj = traj_start is not None or traj_stop is not None
        traj_values: list[int] | None = None
        if has_traj:
            t_start = traj_start if traj_start is not None else 0
            t_stop = traj_stop if traj_stop is not None else len(vis)
            traj_values = list(range(t_start, t_stop, traj_step))
            if not traj_values:
                msg = "Trajectory range is empty"
                raise typer.BadParameter(msg)

        # Build curve range
        has_curve = orbit or curve is not None or curve_step > 0
        curve_values: list[float] | None = None
        if has_curve:
            if curve_step <= 0:
                msg = "--curve-step must be > 0 when using camera curves"
                raise typer.BadParameter(msg)
            import numpy as np

            curve_values = np.arange(
                curve_start, curve_stop + curve_step / 2, curve_step
            ).clip(0.0, 1.0).tolist()
            if not curve_values:
                msg = "Curve range is empty"
                raise typer.BadParameter(msg)

        schedule = _build_schedule(traj_values, curve_values)

        # Save originals for cleanup
        original_step = vis.step
        original_camera_key = session.active_camera
        temp_keys: list[str] = []

        try:
            # Set up camera curve if needed
            camera_curve_key: str | None = None
            if orbit:
                from zndraw.geometries.circle_curve import CircleCurve
                from zndraw.geometries.camera import Camera
                from zndraw.transformations import CurveAttachment

                # Determine orbit center
                current_cam = session.camera
                if center is not None:
                    orbit_center = _parse_center(center)
                elif isinstance(current_cam.target, tuple):
                    orbit_center = current_cam.target
                else:
                    orbit_center = (0.0, 0.0, 0.0)

                # Create temporary circle curve (XZ plane for horizontal orbit)
                curve_key = f"_gif_orbit_{uuid.uuid4().hex[:8]}"
                vis.geometries[curve_key] = CircleCurve(
                    position=[orbit_center],
                    radius=radius,
                    rotation=[(math.pi / 2, 0.0, 0.0)],
                    active=False,
                )
                temp_keys.append(curve_key)
                camera_curve_key = curve_key

            elif curve is not None:
                camera_curve_key = curve

            # Create temp camera if curve mode
            temp_camera_key: str | None = None
            if camera_curve_key is not None:
                from zndraw.geometries.camera import Camera
                from zndraw.transformations import CurveAttachment

                current_cam = session.camera
                cam_key = f"_gif_cam_{uuid.uuid4().hex[:8]}"
                target = current_cam.target
                if isinstance(target, CurveAttachment):
                    target = (0.0, 0.0, 0.0)

                vis.geometries[cam_key] = Camera(
                    position=CurveAttachment(
                        geometry_key=camera_curve_key,
                        progress=0.0,
                    ),
                    target=target,
                    up=current_cam.up,
                    fov=current_cam.fov,
                    near=current_cam.near,
                    far=current_cam.far,
                    helper_visible=False,
                    active=True,
                )
                temp_keys.append(cam_key)
                temp_camera_key = cam_key
                session.active_camera = cam_key

            # Capture frames
            frames: list[bytes] = []
            total = len(schedule)
            prev_traj: int | None = None
            prev_curve: float | None = None

            for idx, (t_val, c_val) in enumerate(schedule):
                # Update trajectory
                if t_val is not None and t_val != prev_traj:
                    vis.step = t_val
                    prev_traj = t_val

                # Update camera progress
                if c_val is not None and c_val != prev_curve and temp_camera_key is not None:
                    from zndraw.geometries.camera import Camera
                    from zndraw.transformations import CurveAttachment

                    cam = vis.geometries[temp_camera_key]
                    assert isinstance(cam, Camera)
                    # Camera is frozen, so create a new instance
                    new_pos = CurveAttachment(
                        geometry_key=camera_curve_key,
                        progress=c_val,
                    )
                    vis.geometries[temp_camera_key] = cam.model_copy(
                        update={"position": new_pos}
                    )
                    prev_curve = c_val

                time.sleep(delay)
                img = session.screenshot()
                frames.append(img.data)
                sys.stderr.write(f"\rCapturing frame {idx + 1}/{total}")
                sys.stderr.flush()

            sys.stderr.write("\n")

            if not frames:
                die("No Frames", "No frames were captured", 400, EXIT_CLIENT_ERROR)

            _assemble_gif(frames, output, fps)
            sys.stderr.write(f"Saved {len(frames)} frames to {output}\n")

        finally:
            # Restore original state
            session.active_camera = original_camera_key
            vis.step = original_step

            # Clean up temp geometries
            for key in temp_keys:
                try:
                    del vis.geometries[key]
                except Exception:
                    pass

            vis.disconnect()
