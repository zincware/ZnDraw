"""Unit tests for GIF export CLI module."""

from __future__ import annotations

import io
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, PropertyMock, patch

import pytest

from zndraw.cli_agent.gif import _assemble_gif, _build_schedule, _get_session

if TYPE_CHECKING:
    from pathlib import Path

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _minimal_png(color: tuple[int, int, int] = (255, 0, 0)) -> bytes:
    """Generate a minimal valid 2x2 PNG image."""
    from PIL import Image

    buf = io.BytesIO()
    img = Image.new("RGB", (2, 2), color=color)
    img.save(buf, format="PNG")
    return buf.getvalue()


def _make_vis_and_session(*, step_val: int = 0):
    """Create a mock vis + session pair for capture tests.

    Returns
    -------
    vis, session, step_box, geom_store, geom_writes, disconnect_calls
    """
    from zndraw.geometries.camera import Camera

    vis = MagicMock()
    vis.__len__ = MagicMock(return_value=100)

    session = MagicMock()
    session.active_camera = "default-cam"
    session.camera = Camera(position=(10.0, 10.0, 10.0), target=(0.0, 0.0, 0.0))

    vis.sessions.__iter__ = MagicMock(return_value=iter(["sid1"]))
    vis.sessions.__getitem__ = MagicMock(return_value=session)

    # Track step changes
    step_box = [step_val]
    type(vis).step = PropertyMock(
        fget=lambda _self: step_box[0],
        fset=lambda _self, v: step_box.__setitem__(0, v),
    )

    # Track geometry writes/deletes (MagicMock dunder tracking is unreliable)
    geom_store: dict[str, object] = {}
    geom_writes: list[tuple[str, object]] = []

    def set_geom(key, value):
        geom_store[key] = value
        geom_writes.append((key, value))

    def get_geom(key):
        from zndraw.geometries.camera import Camera as Cam

        return geom_store.get(key, Cam())

    def del_geom(key):
        geom_store.pop(key, None)

    vis.geometries.__setitem__ = MagicMock(side_effect=set_geom)
    vis.geometries.__getitem__ = MagicMock(side_effect=get_geom)
    vis.geometries.__delitem__ = MagicMock(side_effect=del_geom)

    # Explicit disconnect tracking
    disconnect_calls: list[bool] = []
    vis.disconnect = MagicMock(side_effect=lambda: disconnect_calls.append(True))

    screenshot_mock = MagicMock()
    screenshot_mock.data = _minimal_png()
    session.screenshot.return_value = screenshot_mock

    return vis, session, step_box, geom_store, geom_writes, disconnect_calls


# ---------------------------------------------------------------------------
# _build_schedule
# ---------------------------------------------------------------------------


def test_build_schedule_traj_only():
    """Traj range, no curve -> schedule has (frame, None) tuples."""
    schedule = _build_schedule([0, 1, 2, 3], None)
    assert len(schedule) == 4
    assert schedule[0] == (0, None)
    assert schedule[3] == (3, None)
    assert all(c is None for _, c in schedule)


def test_build_schedule_curve_only():
    """Curve range, no traj -> schedule has (None, progress) tuples."""
    schedule = _build_schedule(None, [0.0, 0.5, 1.0])
    assert len(schedule) == 3
    assert schedule[0] == (None, 0.0)
    assert schedule[2] == (None, 1.0)
    assert all(t is None for t, _ in schedule)


def test_build_schedule_traj_shorter_stretches():
    """Traj shorter than curve -> traj values stretched via linspace."""
    schedule = _build_schedule([0, 1], [0.0, 0.25, 0.5, 0.75, 1.0])
    assert len(schedule) == 5
    # np.linspace(0, 1, 5).astype(int) = [0, 0, 0, 0, 1]
    assert schedule[0] == (0, 0.0)
    assert schedule[1] == (0, 0.25)
    assert schedule[2] == (0, 0.5)
    assert schedule[3] == (0, 0.75)
    assert schedule[4] == (1, 1.0)


def test_build_schedule_both_axes_equal():
    """Same length -> 1:1 zip."""
    schedule = _build_schedule([10, 20, 30], [0.0, 0.5, 1.0])
    assert len(schedule) == 3
    assert schedule[0] == (10, 0.0)
    assert schedule[1] == (20, 0.5)
    assert schedule[2] == (30, 1.0)


def test_build_schedule_neither_axis():
    """No axes -> single still frame."""
    schedule = _build_schedule(None, None)
    assert schedule == [(None, None)]


def test_build_schedule_curve_shorter_stretches():
    """Curve shorter than traj -> curve values stretched via linspace."""
    schedule = _build_schedule([0, 1, 2, 3], [0.0, 1.0])
    assert len(schedule) == 4
    # np.linspace(0, 1, 4).astype(int) = [0, 0, 0, 1]
    assert schedule[0] == (0, 0.0)
    assert schedule[1] == (1, 0.0)
    assert schedule[2] == (2, 0.0)
    assert schedule[3] == (3, 1.0)


# ---------------------------------------------------------------------------
# _assemble_gif
# ---------------------------------------------------------------------------


def test_assemble_gif_creates_valid_file(tmp_path: Path):
    """Minimal PNG bytes -> _assemble_gif -> verify valid multi-frame GIF."""
    from PIL import Image

    # Use different colors so PIL doesn't optimize away identical frames
    frames = [
        _minimal_png((255, 0, 0)),
        _minimal_png((0, 255, 0)),
        _minimal_png((0, 0, 255)),
    ]
    output = tmp_path / "test.gif"
    _assemble_gif(frames, output, fps=10)

    assert output.exists()
    assert output.read_bytes()[:3] == b"GIF"

    img = Image.open(output)
    assert img.format == "GIF"
    assert getattr(img, "n_frames", 1) == 3


def test_assemble_gif_single_frame(tmp_path: Path):
    """Single frame produces a valid GIF."""
    output = tmp_path / "single.gif"
    _assemble_gif([_minimal_png()], output, fps=20)
    assert output.exists()
    assert output.read_bytes()[:3] == b"GIF"


def test_assemble_gif_fps_affects_duration(tmp_path: Path):
    """Different FPS should produce different frame durations."""
    from PIL import Image

    frames = [_minimal_png() for _ in range(2)]

    out_slow = tmp_path / "slow.gif"
    _assemble_gif(frames, out_slow, fps=5)

    out_fast = tmp_path / "fast.gif"
    _assemble_gif(frames, out_fast, fps=50)

    slow_duration = Image.open(out_slow).info.get("duration", 0)
    fast_duration = Image.open(out_fast).info.get("duration", 0)
    assert slow_duration > fast_duration


# ---------------------------------------------------------------------------
# _get_session
# ---------------------------------------------------------------------------


def test_get_session_returns_first():
    """Returns the first active session."""
    vis = MagicMock()
    session_mock = MagicMock()
    vis.sessions.__iter__ = MagicMock(return_value=iter(["sid1", "sid2"]))
    vis.sessions.__getitem__ = MagicMock(return_value=session_mock)

    result = _get_session(vis)
    assert result is session_mock
    vis.sessions.__getitem__.assert_called_once_with("sid1")


def test_get_session_exits_when_no_sessions():
    """Exits with error when vis.sessions is empty."""
    vis = MagicMock()
    vis.sessions.__iter__ = MagicMock(return_value=iter([]))

    with pytest.raises(SystemExit):
        _get_session(vis)


# ---------------------------------------------------------------------------
# capture command — via CliRunner
# ---------------------------------------------------------------------------


def test_capture_orbit_and_curve_mutually_exclusive():
    """Error when both --orbit and --curve specified."""
    from typer.testing import CliRunner

    from zndraw.cli_agent.gif import gif_app

    runner = CliRunner()
    result = runner.invoke(
        gif_app,
        [
            "--orbit",
            "--curve",
            "some-key",
            "--curve-step",
            "0.1",
            "--room",
            "test-room",
        ],
    )
    assert result.exit_code != 0


@patch("zndraw.cli_agent.gif.get_zndraw")  # why: orchestration test of geometry creation/cleanup, not connectivity
@patch("zndraw.cli_agent.gif.resolve_room", return_value="test-room")  # why: orchestration test of geometry creation/cleanup, not connectivity
def test_capture_orbit_creates_temp_geometries(
    mock_resolve, mock_get_zndraw, tmp_path: Path
):
    """Verify CircleCurve + Camera created in geom_store when --orbit."""
    from typer.testing import CliRunner

    from zndraw.cli_agent.gif import gif_app

    vis, _session, _step_box, geom_store, geom_writes, _disconnect_calls = (
        _make_vis_and_session()
    )
    mock_get_zndraw.return_value = vis

    output = tmp_path / "orbit.gif"
    runner = CliRunner()
    result = runner.invoke(
        gif_app,
        [
            "-o",
            str(output),
            "--orbit",
            "--radius",
            "20",
            "--curve-step",
            "0.5",
            "--delay",
            "0.0",
            "--room",
            "test-room",
        ],
    )

    assert result.exit_code == 0, result.output

    # Check that geometries were written (tracked via side_effect)
    written_types = {type(v).__name__ for _, v in geom_writes}
    assert "CircleCurve" in written_types
    assert "Camera" in written_types

    # Cleanup should have deleted temp keys (geom_store empty)
    assert len(geom_store) == 0


@patch("zndraw.cli_agent.gif.get_zndraw")  # why: orchestration test of geometry creation/cleanup, not connectivity
@patch("zndraw.cli_agent.gif.resolve_room", return_value="test-room")  # why: orchestration test of geometry creation/cleanup, not connectivity
def test_capture_restores_step(mock_resolve, mock_get_zndraw, tmp_path: Path):
    """Verify step is restored to its original value after capture."""
    from typer.testing import CliRunner

    from zndraw.cli_agent.gif import gif_app

    vis, _session, step_box, _geom_store, _geom_writes, _disconnect_calls = (
        _make_vis_and_session(step_val=7)
    )
    mock_get_zndraw.return_value = vis

    output = tmp_path / "restore.gif"
    runner = CliRunner()
    result = runner.invoke(
        gif_app,
        [
            "-o",
            str(output),
            "--traj-start",
            "0",
            "--traj-stop",
            "3",
            "--delay",
            "0.0",
            "--room",
            "test-room",
        ],
    )

    assert result.exit_code == 0, result.output
    # Step should be restored to original value 7
    assert step_box[0] == 7


@patch("zndraw.cli_agent.gif.get_zndraw")  # why: orchestration test of geometry creation/cleanup, not connectivity
@patch("zndraw.cli_agent.gif.resolve_room", return_value="test-room")  # why: orchestration test of geometry creation/cleanup, not connectivity
def test_capture_cleans_up_on_error(mock_resolve, mock_get_zndraw, tmp_path: Path):
    """Verify cleanup runs even when screenshot raises."""
    from typer.testing import CliRunner

    from zndraw.cli_agent.gif import gif_app

    vis, session, step_box, _geom_store, _geom_writes, disconnect_calls = (
        _make_vis_and_session(step_val=5)
    )
    mock_get_zndraw.return_value = vis

    # Screenshot raises
    session.screenshot.side_effect = TimeoutError("test error")

    output = tmp_path / "error.gif"
    runner = CliRunner()
    result = runner.invoke(
        gif_app,
        [
            "-o",
            str(output),
            "--traj-start",
            "0",
            "--traj-stop",
            "3",
            "--delay",
            "0.0",
            "--room",
            "test-room",
        ],
    )

    assert result.exit_code != 0

    # Step should be restored to original value despite error
    assert step_box[0] == 5

    # disconnect should have been called (tracked via side_effect)
    assert len(disconnect_calls) == 1


def test_pillow_missing_gives_clear_error(tmp_path: Path):
    """Verify exit when Pillow import fails."""
    with (
        patch.dict("sys.modules", {"PIL": None, "PIL.Image": None}),  # why: tests error message when PIL is not installed
        pytest.raises(SystemExit),
    ):
        _assemble_gif([b"fake"], tmp_path / "test.gif", 20)
