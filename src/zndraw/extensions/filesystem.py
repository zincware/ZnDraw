"""Filesystem extension for loading files from registered providers."""

import typing as t
from pathlib import Path

from fsspec.implementations.dirfs import DirFileSystem
from fsspec.implementations.local import LocalFileSystem

from zndraw.extensions.abc import Category, Extension
from zndraw.io import open_frames


class LoadFile(Extension):
    """Load a file from a registered filesystem provider into the room."""

    category: t.ClassVar[Category] = Category.MODIFIER

    provider_name: str
    path: str
    start: int | None = None
    stop: int | None = None
    step: int | None = None

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        """Read frames via ``open_frames`` and extend the room.

        Uses the same loader as the ``zndraw <file>`` CLI so random-access
        formats (``.h5``, ``.h5md``, ``.lmdb``) are supported alongside
        everything ``ase.io.iread`` can stream.
        """
        providers: dict[str, t.Any] = kwargs.get("providers") or {}
        fs = providers.get(self.provider_name)
        if fs is None:
            available = list(providers)
            raise ValueError(
                f"Provider '{self.provider_name}' not found. "
                f"Available providers: {available}"
            )

        local_path = _require_local_path(fs, self.path)
        frames = open_frames(
            local_path,
            start=self.start,
            stop=self.stop,
            step=self.step,
        )
        vis.extend(frames)


def _require_local_path(fs: t.Any, path: str) -> str:
    """Resolve an fsspec filesystem + path to a local filesystem path.

    Random-access formats like ``.h5`` can't be opened through an
    ``fs.open`` handle — ``asebytes.ASEIO`` needs a real on-disk path.
    The only handler the server materialises is a ``DirFileSystem`` over
    ``LocalFileSystem`` (for the ``@internal`` provider), so we cover
    that and bare ``LocalFileSystem``.
    """
    inner = fs.fs if isinstance(fs, DirFileSystem) else fs
    if not isinstance(inner, LocalFileSystem):
        raise NotImplementedError(
            f"LoadFile requires a local filesystem provider; "
            f"got backend {type(inner).__name__}"
        )
    if isinstance(fs, DirFileSystem):
        return str(Path(fs.path) / path.lstrip("/"))
    return path
