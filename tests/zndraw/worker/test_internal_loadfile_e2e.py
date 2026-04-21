"""End-to-end test for the @internal:modifiers:LoadFile dispatch path.

REGRESSION PIN: the original e2e bug was that ``InternalExtensionExecutor``
calls ``instance.run(vis)`` with no kwargs, while ``LoadFile.run`` looks up
its filesystem handler in ``kwargs.get("providers")``. No code path ever
injected the server-side @internal filesystem handlers into that dict, so
the task failed at runtime with::

    ValueError: Provider '@internal:filesystem:FilesystemRead' not found.
    Available providers: []

ADDITIONAL PIN (#923): ``LoadFile`` used to open the provider file as a
text stream and call ``ase.io.read`` on the handle, which can't read the
random-access formats (``.h5`` / ``.h5md`` / ``.lmdb``) that the
``zndraw <file>`` CLI supports. Every format is parametrised so the full
format matrix runs through the real server + dispatch path.
"""

from __future__ import annotations

from zndraw import ZnDraw
from zndraw.extensions.filesystem import LoadFile


def test_load_file_e2e_via_internal_dispatch(server_factory, water_file):
    instance = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": str(water_file.parent.resolve()),
        }
    )

    vis = ZnDraw(url=instance.url)
    try:
        task = vis.run(
            LoadFile(
                provider_name="@internal:filesystem:FilesystemRead",
                path=f"/{water_file.name}",
            )
        )
        task.wait(timeout=30)

        assert task.status == "completed", f"expected completed, got {task.status!r}"

        assert len(vis) >= 1
        loaded = vis[-1]
        assert len(loaded) == 3
        assert loaded.get_chemical_formula() == "H2O"
    finally:
        vis.disconnect()


def test_load_file_honours_slice_e2e(server_factory, water_trajectory_h5):
    """LoadFile must apply ``start``/``stop``/``step`` on asebytes formats too."""
    instance = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": str(water_trajectory_h5.parent.resolve()),
        }
    )

    vis = ZnDraw(url=instance.url)
    try:
        task = vis.run(
            LoadFile(
                provider_name="@internal:filesystem:FilesystemRead",
                path=f"/{water_trajectory_h5.name}",
                start=1,
                stop=4,
                step=2,
            )
        )
        task.wait(timeout=30)

        assert task.status == "completed", f"expected completed, got {task.status!r}"
        # slice(1, 4, 2) over 5 frames -> 2 frames
        assert len(vis) == 2
    finally:
        vis.disconnect()
