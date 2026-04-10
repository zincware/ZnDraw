"""E2E tests for #906 — global extension cleanup after worker disconnect.

Covers the two acceptance criteria:

1. After ``jobs.disconnect()``, ``GET /v1/joblib/rooms/@global/jobs`` no
   longer includes the extension (graceful path).
2. After SIGKILL, the background sweeper removes the extension within
   ``worker_timeout_seconds`` + one sweeper interval.

The SIGKILL test overrides the sweeper timings via ``server_factory`` so
it runs in a few seconds instead of the ~90s it would take with
production defaults. The subprocess helper lives in
``_sigkill_worker_child.py`` and is owned by the
``spawn_sigkill_worker`` fixture.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import time
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import ClassVar

import pytest

from zndraw import ZnDraw
from zndraw_joblib.client import Category, Extension

_SIGKILL_HELPER = Path(__file__).parent / "_sigkill_worker_child.py"


class NoopCleanup(Extension):
    """No-op modifier registered by the graceful-disconnect test."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self, vis, **_kwargs):  # pragma: no cover - never executed
        pass


def _list_global_names(client: ZnDraw) -> list[str]:
    resp = client.api.http.get(
        f"{client.api.base_url}/v1/joblib/rooms/@global/jobs",
        headers=client.api.get_headers(),
    )
    resp.raise_for_status()
    return [item["name"] for item in resp.json()["items"]]


def _poll_until_gone(client: ZnDraw, ext_name: str, *, max_seconds: float) -> float:
    """Poll the ``@global`` jobs endpoint until *ext_name* disappears.

    Returns elapsed seconds. Raises ``AssertionError`` on timeout.
    """
    start = time.monotonic()
    while time.monotonic() - start < max_seconds:
        if ext_name not in _list_global_names(client):
            return time.monotonic() - start
        time.sleep(0.2)
    raise AssertionError(
        f"Extension {ext_name!r} still listed after {max_seconds:.1f}s"
    )


# =============================================================================
# Fixtures
# =============================================================================


SigkillWorkerSpawn = Callable[[str], "tuple[subprocess.Popen[str], str]"]


@pytest.fixture
def spawn_sigkill_worker() -> Iterator[SigkillWorkerSpawn]:
    """Spawn the ``_sigkill_worker_child.py`` helper as a subprocess.

    Yields a factory ``(server_url) -> (proc, worker_id)`` that launches
    the helper and blocks until it prints ``READY <worker_id>``. Any
    processes spawned during the test are killed at teardown as a
    belt-and-braces cleanup in case the test itself fails before
    reaching its own ``kill``.
    """
    procs: list[subprocess.Popen[str]] = []

    def _spawn(server_url: str) -> tuple[subprocess.Popen[str], str]:
        # sys.executable + constant helper path + test-fixture URL — S603 is safe
        proc = subprocess.Popen(  # noqa: S603
            [sys.executable, str(_SIGKILL_HELPER), server_url],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        procs.append(proc)
        assert proc.stdout is not None

        deadline = time.time() + 20.0
        while time.time() < deadline:
            line = proc.stdout.readline()
            if not line:
                if proc.poll() is not None:
                    err = proc.stderr.read() if proc.stderr else ""
                    pytest.fail(f"sigkill worker child exited before READY:\n{err}")
                continue
            if line.startswith("READY "):
                return proc, line.strip().split(" ", 1)[1]
        pytest.fail("sigkill worker child did not print READY within 20s")

    yield _spawn

    for p in procs:
        if p.poll() is None:
            p.kill()
            p.wait(timeout=5)


# =============================================================================
# Test 1: graceful disconnect clears the listing (happy path)
# =============================================================================


def test_graceful_disconnect_removes_global_extension(server):
    """Registering then disconnecting a worker removes the extension.

    The extension must no longer appear in
    ``GET /v1/joblib/rooms/@global/jobs`` within 2 seconds of
    ``jobs.disconnect()``. Matches the first acceptance bullet of #906.
    """
    worker = ZnDraw(url=server)
    observer = ZnDraw(url=server)
    try:
        worker.jobs.register(NoopCleanup)
        assert "NoopCleanup" in _list_global_names(observer)

        worker.jobs.disconnect()
        elapsed = _poll_until_gone(observer, "NoopCleanup", max_seconds=2.0)
        assert elapsed < 2.0
    finally:
        worker.disconnect()
        observer.disconnect()


# =============================================================================
# Test 2: SIGKILL cleanup via background sweeper
# =============================================================================


def test_sigkill_worker_cleared_by_sweeper(server_factory, spawn_sigkill_worker):
    """SIGKILL'd worker is cleaned up by the sweeper within the configured window.

    Uses overridden sweeper settings so the test runs in ~5s instead of
    ~90s. Matches the second acceptance bullet of #906.
    """
    instance = server_factory(
        {
            "ZNDRAW_JOBLIB_WORKER_TIMEOUT_SECONDS": "3",
            "ZNDRAW_JOBLIB_SWEEPER_INTERVAL_SECONDS": "1",
        }
    )
    observer = ZnDraw(url=instance.url)

    try:
        proc, _worker_id = spawn_sigkill_worker(instance.url)
        assert "SigkillCleanup" in _list_global_names(observer), (
            "extension not registered by child"
        )

        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(timeout=5)

        # Budget: worker_timeout (3s) + sweeper_interval (1s) + generous slack
        _poll_until_gone(observer, "SigkillCleanup", max_seconds=15.0)
    finally:
        observer.disconnect()
