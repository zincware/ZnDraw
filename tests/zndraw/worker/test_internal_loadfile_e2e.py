"""End-to-end test for the @internal:modifiers:LoadFile dispatch path.

REGRESSION PIN (#923): ``LoadFile`` used to open the provider file as a
text stream and call ``ase.io.read`` on the handle, which can't read the
random-access formats (``.h5`` / ``.h5md`` / ``.lmdb``) that the
``zndraw <file>`` CLI supports. Every format is parametrised so the full
format matrix runs through the real server + dispatch path.
"""

from zndraw import ZnDraw
from zndraw.extensions.filesystem import LoadFile


def test_load_file_e2e_via_internal_dispatch(server_factory, water_file):
    path, frames = water_file
    instance = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": str(path.parent.resolve()),
        }
    )

    vis = ZnDraw(url=instance.url)
    try:
        task = vis.run(
            LoadFile(
                provider_name="@internal:filesystem:FilesystemRead",
                path=f"/{path.name}",
            )
        )
        task.wait(timeout=30)

        assert task.status == "completed", f"expected completed, got {task.status!r}"
        assert list(vis) == frames
    finally:
        vis.disconnect()


def test_load_file_honours_slice_e2e(server_factory, water_trajectory_h5):
    """LoadFile must apply ``start``/``stop``/``step`` on asebytes formats too."""
    path, frames = water_trajectory_h5
    instance = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": str(path.parent.resolve()),
        }
    )

    vis = ZnDraw(url=instance.url)
    try:
        task = vis.run(
            LoadFile(
                provider_name="@internal:filesystem:FilesystemRead",
                path=f"/{path.name}",
                start=1,
                stop=4,
                step=2,
            )
        )
        task.wait(timeout=30)

        assert task.status == "completed", f"expected completed, got {task.status!r}"
        assert list(vis) == frames[1:4:2]
    finally:
        vis.disconnect()
