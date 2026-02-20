"""Tests for client-side source mount/unmount and fetch handler."""

import time

import ase
import numpy as np
import pytest
from ase.data import covalent_radii

from zndraw.client import ZnDraw

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class FakeSource:
    """Minimal FrameSource for testing."""

    def __init__(self, frames: list[ase.Atoms]) -> None:
        self._frames = frames

    def __len__(self) -> int:
        return len(self._frames)

    def __getitem__(self, index: int) -> ase.Atoms:
        return self._frames[index]


def _make_atoms(n: int = 1) -> ase.Atoms:
    return ase.Atoms("H" * n, positions=np.zeros((n, 3)))


# ---------------------------------------------------------------------------
# vis.mount()
# ---------------------------------------------------------------------------


def test_mount_stores_source(server: str) -> None:
    """mount() stores the source and updates cached length."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms() for _ in range(10)])
        vis.mount(source)

        assert vis._mount is source
        assert len(vis) == 10
    finally:
        vis.disconnect()


def test_mount_on_non_empty_room_raises(server: str) -> None:
    """mount() raises when room already has frames."""
    vis = ZnDraw(url=server)
    try:
        vis.append(_make_atoms())
        source = FakeSource([_make_atoms()])

        with pytest.raises(RuntimeError, match="empty"):
            vis.mount(source)
    finally:
        vis.disconnect()


def test_mount_twice_raises(server: str) -> None:
    """mount() raises when a source is already mounted."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms() for _ in range(5)])
        vis.mount(source)

        with pytest.raises(RuntimeError, match="already"):
            vis.mount(source)
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# vis.unmount()
# ---------------------------------------------------------------------------


def test_unmount_clears_mount(server: str) -> None:
    """unmount() clears the stored source and resets room."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms() for _ in range(10)])
        vis.mount(source)
        vis.unmount()

        assert vis._mount is None
        assert len(vis) == 0
    finally:
        vis.disconnect()


def test_unmount_without_mount_raises(server: str) -> None:
    """unmount() raises when no source is mounted."""
    vis = ZnDraw(url=server)
    try:
        with pytest.raises(RuntimeError, match="No mount"):
            vis.unmount()
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Full fetch round-trip: mount → GET frame → SIO fetch → upload → retry
# ---------------------------------------------------------------------------


def test_fetch_round_trip(server: str) -> None:
    """After mount, reading a frame triggers fetch and returns data."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms(i + 1) for i in range(5)])
        vis.mount(source)

        result = _wait_for_frame(vis, 0)
        assert isinstance(result, ase.Atoms)
        assert len(result) == 1
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Helper: wait for frame fetch round-trip
# ---------------------------------------------------------------------------


def _wait_for_frame(vis: ZnDraw, index: int, timeout: float = 5.0) -> ase.Atoms:
    """Retry vis[index] until it succeeds or timeout."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            return vis[index]
        except IndexError:
            time.sleep(0.2)
    raise TimeoutError(f"Frame {index} not available within {timeout}s")


# ---------------------------------------------------------------------------
# Write rejection on mounted room
# ---------------------------------------------------------------------------


def test_append_on_mounted_room_raises(server: str) -> None:
    """Appending to a mounted room raises (room is read-only)."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms() for _ in range(5)])
        vis.mount(source)

        with pytest.raises(Exception):
            vis.append(_make_atoms())
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Fetch multiple frames with correct data
# ---------------------------------------------------------------------------


def test_fetch_returns_correct_atom_counts(server: str) -> None:
    """Different frames have correct number of atoms after fetch."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms(i + 1) for i in range(3)])
        vis.mount(source)

        frame0 = _wait_for_frame(vis, 0)
        assert len(frame0) == 1

        frame1 = _wait_for_frame(vis, 1)
        assert len(frame1) == 2

        frame2 = _wait_for_frame(vis, 2)
        assert len(frame2) == 3
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Mount lifecycle: mount → read → unmount → empty
# ---------------------------------------------------------------------------


def test_mount_read_unmount_lifecycle(server: str) -> None:
    """Full lifecycle: mount, read frames, unmount, verify empty."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms(2) for _ in range(3)])
        vis.mount(source)
        assert len(vis) == 3

        frame = _wait_for_frame(vis, 0)
        assert len(frame) == 2

        vis.unmount()
        assert len(vis) == 0
        assert vis._mount is None
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Length update via PATCH /rooms/{room_id}
# ---------------------------------------------------------------------------


def test_update_frame_count_via_room_patch(server: str) -> None:
    """PATCH room with frame_count updates len(vis)."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms() for _ in range(5)])
        vis.mount(source)
        assert len(vis) == 5

        vis.api.update_room({"frame_count": 10})
        # Server stores new count in Redis; client re-queries length
        vis._cached_length = None  # force re-query
        assert len(vis) == 10
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Radii and connectivity computed on fetch
# ---------------------------------------------------------------------------


def test_copy_from_mounted_room_raises(server: str) -> None:
    """Copying from a room with a mounted source returns 409."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms() for _ in range(5)])
        vis.mount(source)

        import uuid

        import httpx

        response = httpx.post(
            f"{server}/v1/rooms",
            headers={"Authorization": f"Bearer {vis.api.token}"},
            json={"room_id": str(uuid.uuid4()), "copy_from": vis.room},
        )
        assert response.status_code == 409
        body = response.json()
        assert body["type"] == "/v1/problems/room-read-only"
    finally:
        vis.disconnect()


def test_fetched_frame_has_radii(server: str) -> None:
    """Fetched frames from source have radii even when source omits them."""
    vis = ZnDraw(url=server)
    try:
        # Bare atoms — no radii, no colors
        atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        assert "radii" not in atoms.arrays

        source = FakeSource([atoms])
        vis.mount(source)

        frame = _wait_for_frame(vis, 0)
        assert "radii" in frame.arrays
        expected = covalent_radii[frame.numbers].astype(np.float32)
        np.testing.assert_array_equal(frame.arrays["radii"], expected)
    finally:
        vis.disconnect()


def test_fetched_frame_has_colors(server: str) -> None:
    """Fetched frames from source have colors even when source omits them."""
    vis = ZnDraw(url=server)
    try:
        atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        assert "colors" not in atoms.arrays

        source = FakeSource([atoms])
        vis.mount(source)

        frame = _wait_for_frame(vis, 0)
        assert "colors" in frame.arrays
        assert len(frame.arrays["colors"]) == 3
    finally:
        vis.disconnect()


def test_fetched_frame_has_connectivity(server: str) -> None:
    """Fetched frames from source have connectivity computed."""
    vis = ZnDraw(url=server)
    try:
        # Two H atoms within bonding distance
        atoms = ase.Atoms("H2", positions=[[0, 0, 0], [0.5, 0, 0]])
        assert "connectivity" not in atoms.info

        source = FakeSource([atoms])
        vis.mount(source)

        frame = _wait_for_frame(vis, 0)
        assert "connectivity" in frame.info
        assert len(frame.info["connectivity"]) >= 1
    finally:
        vis.disconnect()


def test_fetched_frame_preserves_existing_radii(server: str) -> None:
    """Source-provided radii are preserved (not overwritten)."""
    vis = ZnDraw(url=server)
    try:
        atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
        custom_radii = np.array([0.5, 0.5], dtype=np.float32)
        atoms.set_array("radii", custom_radii)

        source = FakeSource([atoms])
        vis.mount(source)

        frame = _wait_for_frame(vis, 0)
        np.testing.assert_array_equal(frame.arrays["radii"], custom_radii)
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Deep copy via copyFrom
# ---------------------------------------------------------------------------


def test_copy_from_room_copies_bookmarks(server: str) -> None:
    """Creating a room with copyFrom deep-copies bookmarks."""
    import uuid

    import httpx

    vis = ZnDraw(url=server)
    try:
        vis.append(_make_atoms())
        vis.bookmarks[0] = "start"

        new_room_id = str(uuid.uuid4())
        response = httpx.post(
            f"{server}/v1/rooms",
            headers={"Authorization": f"Bearer {vis.api.token}"},
            json={"room_id": new_room_id, "copy_from": vis.room},
        )
        assert response.status_code == 201

        copy = ZnDraw(url=server, room=new_room_id)
        try:
            assert len(copy.bookmarks) == 1
            assert copy.bookmarks[0] == "start"
        finally:
            copy.disconnect()
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Remount returns fresh data (not stale cache)
# ---------------------------------------------------------------------------


def test_provider_disconnect_clears_frame_count(server: str) -> None:
    """When a mounted client disconnects, the room's frame count drops to 0."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms() for _ in range(5)])
        vis.mount(source)
        assert len(vis) == 5
        room_id = vis.room
    finally:
        vis.disconnect()

    # After provider disconnects, an observer should see 0 frames
    observer = ZnDraw(url=server, room=room_id)
    try:
        deadline = time.monotonic() + 5.0
        while time.monotonic() < deadline:
            observer._cached_length = None  # force re-query
            if len(observer) == 0:
                break
            time.sleep(0.3)
        assert len(observer) == 0, (
            f"Expected 0 frames after provider disconnect, got {len(observer)}"
        )
    finally:
        observer.disconnect()


def test_remount_serves_new_source_data(server: str) -> None:
    """After unmount + remount with different source, new data is served."""
    vis = ZnDraw(url=server)
    try:
        # Mount source A: 3 atoms per frame
        source_a = FakeSource([_make_atoms(3)])
        vis.mount(source_a)

        frame_a = _wait_for_frame(vis, 0)
        assert len(frame_a) == 3

        vis.unmount()

        # Mount source B: 7 atoms per frame — different data
        source_b = FakeSource([_make_atoms(7)])
        vis.mount(source_b)

        frame_b = _wait_for_frame(vis, 0)
        assert len(frame_b) == 7, (
            f"Expected 7 atoms from new mount, got {len(frame_b)} (stale cache?)"
        )
    finally:
        vis.disconnect()


# ---------------------------------------------------------------------------
# Slicing on mounted rooms
# ---------------------------------------------------------------------------


def test_slice_on_mounted_room(server: str) -> None:
    """vis[:5] on a mounted room returns 5 frames, not []."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms(i + 1) for i in range(10)])
        vis.mount(source)
        assert len(vis) == 10

        frames = vis[:5]
        assert len(frames) == 5
        for i, frame in enumerate(frames):
            assert isinstance(frame, ase.Atoms)
            assert len(frame) == i + 1
    finally:
        vis.disconnect()


def test_slice_middle_range_on_mounted_room(server: str) -> None:
    """vis[3:7] on a mounted room returns 4 frames."""
    vis = ZnDraw(url=server)
    try:
        source = FakeSource([_make_atoms(i + 1) for i in range(10)])
        vis.mount(source)

        frames = vis[3:7]
        assert len(frames) == 4
        for i, frame in enumerate(frames):
            assert len(frame) == (3 + i) + 1
    finally:
        vis.disconnect()
