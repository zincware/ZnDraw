"""Unified file/URI reading for atomistic data."""

from __future__ import annotations

import typing as t

if t.TYPE_CHECKING:
    from collections.abc import Iterator

    import ase
    import asebytes

_ASEBYTES_EXTENSIONS = frozenset({".h5", ".h5md"})


def open_frames(
    path_or_uri: str,
    *,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
) -> asebytes.ASEIO | Iterator[ase.Atoms]:
    """Open a file for frame access.

    Routes ``.h5`` / ``.h5md`` to ``asebytes.ASEIO`` (random-access).
    Everything else streams via ``ase.io.iread``.

    Parameters
    ----------
    path_or_uri
        File path or URI string.
    start, stop, step
        Optional slice parameters for frame selection.

    Returns
    -------
    asebytes.ASEIO | Iterator[ase.Atoms]
        Random-access source for asebytes formats, streaming iterator
        for everything else.
    """
    from pathlib import Path

    suffix = Path(path_or_uri).suffix.lower()

    if suffix in _ASEBYTES_EXTENSIONS:
        import asebytes as _asebytes

        db = _asebytes.ASEIO(path_or_uri)
        if any(x is not None for x in [start, stop, step]):
            return db[slice(start, stop, step)]
        return db

    return _iread(path_or_uri, start=start, stop=stop, step=step)


def _iread(
    path: str,
    *,
    start: int | None,
    stop: int | None,
    step: int | None,
) -> Iterator[ase.Atoms]:
    """Stream via ``ase.io.iread`` — never loads full file."""
    import ase.io

    index: str | slice = (
        slice(start, stop, step)
        if any(x is not None for x in [start, stop, step])
        else ":"
    )
    yield from ase.io.iread(path, index=index)
