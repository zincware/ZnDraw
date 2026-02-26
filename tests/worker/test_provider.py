"""E2E tests for filesystem provider registration, listing, and read flow.

Tests cover: register_fs → list providers → read (long-poll → 200).
"""

import time
from urllib.parse import urlencode

import fsspec
import pytest
from pydantic import BaseModel
from zndraw_joblib.schemas import JobSummary, ProviderResponse

from zndraw import ZnDraw


class FileItem(BaseModel):
    """Expected shape of a single filesystem provider read result."""

    name: str
    path: str
    size: int
    type: str


# =============================================================================
# Helpers
# =============================================================================


def _list_providers(vis: ZnDraw) -> list[ProviderResponse]:
    """GET /rooms/{room}/providers and validate with Pydantic."""
    resp = vis.api.http.get(
        f"{vis.api.base_url}/v1/joblib/rooms/{vis.room}/providers",
        headers=vis.api.get_headers(),
    )
    resp.raise_for_status()
    return [ProviderResponse.model_validate(p) for p in resp.json()["items"]]


def _list_jobs(vis: ZnDraw) -> list[JobSummary]:
    """GET /rooms/{room}/jobs and validate with Pydantic."""
    resp = vis.api.http.get(
        f"{vis.api.base_url}/v1/joblib/rooms/{vis.room}/jobs",
        headers=vis.api.get_headers(),
    )
    resp.raise_for_status()
    return [JobSummary.model_validate(j) for j in resp.json()["items"]]


def _read_provider(
    vis: ZnDraw, provider_full_name: str, params: dict[str, str] | None = None
) -> tuple[int, list[FileItem]]:
    """GET the provider read endpoint, return (status_code, items).

    Items are empty for non-200 responses.
    """
    qs = "?" + urlencode(params) if params else ""
    resp = vis.api.http.get(
        f"{vis.api.base_url}/v1/joblib/rooms/{vis.room}/providers/"
        f"{provider_full_name}{qs}",
        headers=vis.api.get_headers(),
    )
    if resp.status_code == 200:
        return 200, [FileItem.model_validate(item) for item in resp.json()]
    return resp.status_code, []


def _poll_until_200(
    vis: ZnDraw,
    full_name: str,
    params: dict[str, str],
    *,
    timeout: float = 15,
    interval: float = 0.5,
) -> list[FileItem]:
    """Poll a provider read until 200, return validated items or fail."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        status, items = _read_provider(vis, full_name, params)
        if status == 200:
            return items
        time.sleep(interval)
    pytest.fail(f"Provider read for {params} did not return 200 within {timeout}s")


# =============================================================================
# Registration & Listing
# =============================================================================


def test_register_fs_creates_provider(server):
    """register_fs registers a FilesystemRead provider visible in list."""
    vis = ZnDraw(url=server)
    try:
        vis.register_fs(fsspec.filesystem("memory"), name="mem")

        providers = _list_providers(vis)
        fs_providers = [p for p in providers if p.category == "filesystem"]
        assert len(fs_providers) == 1
        assert fs_providers[0].name == "mem"
        assert fs_providers[0].full_name == f"{vis.room}:filesystem:mem"
    finally:
        vis.jobs.disconnect()
        vis.disconnect()


def test_register_fs_also_registers_load_file_job(server):
    """register_fs registers both the provider and a LoadFile modifier job."""
    vis = ZnDraw(url=server)
    try:
        vis.register_fs(fsspec.filesystem("memory"), name="mem")

        jobs = _list_jobs(vis)
        job_names = {j.name for j in jobs}
        assert "LoadFile" in job_names
    finally:
        vis.jobs.disconnect()
        vis.disconnect()


def test_provider_not_visible_in_other_room(server):
    """A room-scoped provider is not visible from a different room."""
    room_a = ZnDraw(url=server)
    room_b = ZnDraw(url=server)
    try:
        room_a.register_fs(fsspec.filesystem("memory"), name="mem")

        resp = room_b.api.http.get(
            f"{room_b.api.base_url}/v1/joblib/rooms/{room_b.room}/providers",
            headers=room_b.api.get_headers(),
        )
        resp.raise_for_status()
        assert len(resp.json()["items"]) == 0
    finally:
        room_a.jobs.disconnect()
        room_a.disconnect()
        room_b.disconnect()


# =============================================================================
# Read Flow: long-poll → 200
# =============================================================================


def test_provider_read_returns_200(server):
    """Read long-polls until provider delivers, then returns 200."""
    fs = fsspec.filesystem("memory")
    fs.mkdir("/testdir")
    fs.pipe("/testdir/a.txt", b"hello")
    fs.pipe("/testdir/b.txt", b"world")

    vis = ZnDraw(url=server)
    try:
        vis.register_fs(fs, name="local")
        full_name = f"{vis.room}:filesystem:local"

        # Long-poll: server holds connection until provider delivers
        items = _poll_until_200(vis, full_name, {"path": "/testdir"})

        names = {item.name for item in items}
        assert "a.txt" in names
        assert "b.txt" in names
        for item in items:
            assert item.type in ("file", "directory")
    finally:
        vis.jobs.disconnect()
        vis.disconnect()


def test_provider_read_cached_returns_200_immediately(server):
    """Second identical read should return 200 from cache without re-dispatch."""
    fs = fsspec.filesystem("memory")
    fs.mkdir("/cached")
    fs.pipe("/cached/x.txt", b"data")

    vis = ZnDraw(url=server)
    try:
        vis.register_fs(fs, name="cached")
        full_name = f"{vis.room}:filesystem:cached"

        # First read → wait for 200
        _poll_until_200(vis, full_name, {"path": "/cached"})

        # Second read with same params should be an immediate cache hit
        status, items = _read_provider(vis, full_name, {"path": "/cached"})
        assert status == 200
        assert any(item.name == "x.txt" for item in items)
    finally:
        vis.jobs.disconnect()
        vis.disconnect()


def test_provider_read_different_params_dispatches_again(server):
    """Different query params produce different cache keys → separate requests."""
    fs = fsspec.filesystem("memory")
    fs.mkdir("/dir1")
    fs.pipe("/dir1/one.txt", b"1")
    fs.mkdir("/dir2")
    fs.pipe("/dir2/two.txt", b"2")

    vis = ZnDraw(url=server)
    try:
        vis.register_fs(fs, name="multi")
        full_name = f"{vis.room}:filesystem:multi"

        # Read /dir1
        _poll_until_200(vis, full_name, {"path": "/dir1"})

        # Read /dir2 (different cache key)
        items = _poll_until_200(vis, full_name, {"path": "/dir2"})
        assert any(item.name == "two.txt" for item in items)
    finally:
        vis.jobs.disconnect()
        vis.disconnect()


# =============================================================================
# Disconnect & Cleanup
# =============================================================================


def test_provider_disconnect_removes_provider(server):
    """Disconnecting the worker removes the provider from listings."""
    vis = ZnDraw(url=server)
    try:
        vis.register_fs(fsspec.filesystem("memory"), name="temp")
        assert len(_list_providers(vis)) > 0

        vis.jobs.disconnect()

        # After disconnect, provider should be gone
        time.sleep(0.5)
        assert len(_list_providers(vis)) == 0
    finally:
        vis.disconnect()


# =============================================================================
# Auth Mode
# =============================================================================


def test_guest_can_register_fs_provider(server_auth):
    """Guest user can register a filesystem provider (no admin required)."""
    vis = ZnDraw(url=server_auth)
    try:
        vis.register_fs(fsspec.filesystem("memory"), name="guest-fs")

        providers = _list_providers(vis)
        fs_providers = [p for p in providers if p.category == "filesystem"]
        assert len(fs_providers) == 1
        assert fs_providers[0].name == "guest-fs"
    finally:
        vis.jobs.disconnect()
        vis.disconnect()
