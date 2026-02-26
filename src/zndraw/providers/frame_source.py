"""Providers for on-demand frame reads from a FrameSource.

FrameSourceRead: binary provider for individual frame data (msgpack).
FrameSourceLength: JSON provider returning the source length.

Results are cached by the result backend. For frame data this MUST be an
LMDB-backed backend (not Redis) — frame payloads are too large for Redis.
See ``CompositeResultBackend`` in Task 3.2.
"""

import typing as t

import ase
import msgpack
from zndraw_joblib import Provider


@t.runtime_checkable
class FrameSource(t.Protocol):
    """Protocol for objects that can serve as a frame source.

    Any object implementing ``__len__`` and ``__getitem__(int)`` satisfies
    this protocol (e.g. ``asebytes.ASEIO``, ``list[ase.Atoms]``).
    """

    def __len__(self) -> int: ...
    def __getitem__(self, index: int, /) -> ase.Atoms: ...


class FrameSourceLength(Provider):
    """Return the length of a mounted FrameSource.

    The handler must satisfy ``__len__``. Returns a JSON dict with
    a single ``"length"`` key.
    """

    category: t.ClassVar[str] = "frames_meta"

    def read(self, handler: FrameSource) -> dict[str, int]:
        """Return ``{"length": len(handler)}``."""
        return {"length": len(handler)}


class FrameSourceRead(Provider):
    """Read a single frame from a mounted FrameSource.

    The handler must satisfy the ``FrameSource`` protocol (``__len__``,
    ``__getitem__``).  Each ``read()`` returns one msgpack-encoded
    ``RawFrame`` (``dict[bytes, bytes]``) with colors, radii, and
    connectivity added on the fly.
    """

    category: t.ClassVar[str] = "frames"
    content_type: t.ClassVar[str] = "application/x-msgpack"

    index: int = 0

    def read(self, handler: FrameSource) -> bytes:
        """Read frame at *index*, enrich it, and return msgpack bytes."""
        from asebytes import encode

        from zndraw.connectivity import add_connectivity
        from zndraw.enrichment import add_colors, add_radii

        atoms = handler[self.index]

        add_colors(atoms)
        add_radii(atoms)

        if len(atoms) < 100 and "connectivity" not in atoms.info:
            add_connectivity(atoms)

        raw_frame: dict[bytes, bytes] = encode(atoms)
        return msgpack.packb(raw_frame)  # type: ignore[return-value]
