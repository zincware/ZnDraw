"""Tests for mount functionality (lazy frame serving)."""

from __future__ import annotations

import tempfile
from pathlib import Path

import ase
import httpx

from zndraw.client import ZnDraw


def test_mount_serves_frames(server_url: str) -> None:
    """mount should create a room and serve the correct frame count."""
    import asebytes

    with tempfile.TemporaryDirectory() as tmpdir:
        lmdb_path = str(Path(tmpdir) / "test.lmdb")
        db = asebytes.ASEIO(lmdb_path)
        db.append(ase.Atoms("H", positions=[[0.0, 0.0, 0.0]]))
        db.append(ase.Atoms("HO", positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))

        # Reopen read-only for mounting
        db = asebytes.ASEIO(lmdb_path)
        assert len(db) == 2

        vis = ZnDraw(url=server_url)
        vis.mount(db)

        assert vis.room is not None

        # Verify the room is accessible via REST
        with httpx.Client(base_url=server_url, timeout=10.0) as client:
            resp = client.post("/v1/auth/guest")
            token = resp.json()["access_token"]
            headers = {"Authorization": f"Bearer {token}"}

            resp = client.get(f"/v1/rooms/{vis.room}/step", headers=headers)
            assert resp.status_code == 200
            assert resp.json()["total_frames"] == 2

        vis.disconnect()
