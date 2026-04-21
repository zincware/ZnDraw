"""End-to-end test for the @internal:modifiers:LoadFile dispatch path.

REGRESSION PIN: the original e2e bug was that ``InternalExtensionExecutor``
calls ``instance.run(vis)`` with no kwargs, while ``LoadFile.run`` looks up
its filesystem handler in ``kwargs.get("providers")``. No code path ever
injected the server-side @internal filesystem handlers into that dict, so
the task failed at runtime with::

    ValueError: Provider '@internal:filesystem:FilesystemRead' not found.
    Available providers: []

This test submits a real LoadFile task against a real uvicorn server that
has ``FILEBROWSER_PATH`` pointed at a temp dir, then waits for the task to
reach a terminal status. A passing test requires:

1. ``LoadFile`` is registered as ``@internal:modifiers:LoadFile`` (so the
   POST to ``/tasks/@internal:modifiers:LoadFile`` is accepted).
2. ``InternalExtensionExecutor`` resolves the @internal filesystem handler
   and passes it via ``run(vis, providers=...)`` (so ``LoadFile.run`` can
   read the file).
3. The room ends up with frames from the file.
"""

from __future__ import annotations

from zndraw import ZnDraw
from zndraw.extensions.filesystem import LoadFile


def test_load_file_e2e_via_internal_dispatch(server_factory, water_xyz):
    instance = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": str(water_xyz.parent.resolve()),
        }
    )

    vis = ZnDraw(url=instance.url)
    try:
        task = vis.run(
            LoadFile(
                provider_name="@internal:filesystem:FilesystemRead",
                path="/water.xyz",
            )
        )
        task.wait(timeout=30)

        assert task.status == "completed", f"expected completed, got {task.status!r}"

        # LoadFile extends the room with atoms from the file.
        assert len(vis) >= 1
        loaded = vis[-1]
        assert len(loaded) == 3
        assert loaded.get_chemical_formula() == "H2O"
    finally:
        vis.disconnect()
