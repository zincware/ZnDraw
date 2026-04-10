"""Subprocess helper for ``test_global_cleanup_e2e``.

Registers a ``@global`` extension named ``SigkillCleanup`` and sleeps
forever. Prints ``READY <worker_id>`` to stdout once registration
completes; the parent test then sends ``SIGKILL`` to exercise the
sweeper cleanup path.

The leading underscore in the filename prevents pytest from collecting
this module as a test. It is only ever executed as a subprocess via
``python _sigkill_worker_child.py <server_url>``.
"""

from __future__ import annotations

import sys
from typing import ClassVar

from zndraw import ZnDraw
from zndraw_joblib.client import Category, Extension


class SigkillCleanup(Extension):
    """No-op modifier used to register a ``@global`` job."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self, vis, **_kwargs):  # pragma: no cover - never executed
        pass


if __name__ == "__main__":
    server_url = sys.argv[1]
    worker = ZnDraw(url=server_url)
    worker.jobs.register(SigkillCleanup)
    # stdout is the IPC channel back to the parent test process
    print(f"READY {worker.jobs.worker_id}", flush=True)  # noqa: T201
    # vis.wait() blocks on the Socket.IO transport. The parent test will
    # SIGKILL this process, which bypasses Python entirely — no cleanup,
    # no atexit, no try/except needed for the test's purpose.
    worker.wait()
