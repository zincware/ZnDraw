"""Unit tests for LoadFile extension."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import ase
import ase.io
import fsspec
import pytest
from fsspec.implementations.dirfs import DirFileSystem

from zndraw.extensions.filesystem import LoadFile


class _FakeVis:
    """Minimal ``vis`` stand-in: records frames ``extend`` is called with."""

    def __init__(self) -> None:
        self.frames: list[ase.Atoms] = []

    def extend(self, frames: Any) -> None:
        self.frames.extend(frames)


def _water() -> ase.Atoms:
    return ase.Atoms("H2O", positions=[[0, 0, 0], [0, 0, 1], [1, 0, 0]])


def _dir_fs(root: Path) -> DirFileSystem:
    return DirFileSystem(path=str(root), fs=fsspec.filesystem("file"))


def test_load_file_schema() -> None:
    schema = LoadFile.model_json_schema()
    assert "provider_name" in schema["properties"]
    assert "path" in schema["properties"]
    assert "start" in schema["properties"]
    assert "stop" in schema["properties"]
    assert "step" in schema["properties"]


def test_load_file_requires_provider() -> None:
    ext = LoadFile(provider_name="missing", path="/data/file.xyz")
    with pytest.raises(ValueError, match="not found"):
        ext.run(None, providers={})


def test_load_file_category() -> None:
    from zndraw_joblib.client import Category

    assert LoadFile.category == Category.MODIFIER


def test_load_file_reads_xyz(tmp_path: Path) -> None:
    ase.io.write(tmp_path / "water.xyz", _water())
    vis = _FakeVis()

    LoadFile(provider_name="fs", path="/water.xyz").run(
        vis, providers={"fs": _dir_fs(tmp_path)}
    )

    assert len(vis.frames) == 1
    assert vis.frames[0].get_chemical_formula() == "H2O"


@pytest.mark.parametrize("suffix", [".h5", ".h5md", ".lmdb"])
def test_load_file_reads_asebytes_formats(tmp_path: Path, suffix: str) -> None:
    """Regression for #923 — asebytes formats must round-trip through LoadFile.

    Before the fix, ``LoadFile`` opened the file as a text stream and ran it
    through ``ase.io.read``, which can't read the random-access formats
    (``.h5`` / ``.h5md`` / ``.lmdb``) that ``zndraw <file>`` accepts.
    """
    import asebytes

    with asebytes.ASEIO(str(tmp_path / f"trj{suffix}")) as db:
        db.extend([_water(), _water()])

    vis = _FakeVis()
    LoadFile(provider_name="fs", path=f"/trj{suffix}").run(
        vis, providers={"fs": _dir_fs(tmp_path)}
    )

    assert len(vis.frames) == 2
    assert all(a.get_chemical_formula() == "H2O" for a in vis.frames)


def test_load_file_honours_slice(tmp_path: Path) -> None:
    ase.io.write(tmp_path / "trj.xyz", [_water() for _ in range(5)])
    vis = _FakeVis()

    LoadFile(provider_name="fs", path="/trj.xyz", start=1, stop=4, step=2).run(
        vis, providers={"fs": _dir_fs(tmp_path)}
    )

    assert len(vis.frames) == 2


def test_load_file_rejects_non_local_backend() -> None:
    """Remote fsspec backends aren't supported yet — asebytes needs a real path."""
    memfs = fsspec.filesystem("memory")
    memfs.pipe_file("/data/file.xyz", b"ignored")

    with pytest.raises(NotImplementedError, match="local filesystem"):
        LoadFile(provider_name="fs", path="/data/file.xyz").run(
            _FakeVis(), providers={"fs": memfs}
        )
