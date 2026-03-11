"""Tests for zndraw.io — unified file reading."""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING

import ase
import pytest

if TYPE_CHECKING:
    from pathlib import Path

from zndraw.io import open_frames


@pytest.fixture
def xyz_file(tmp_path: Path) -> Path:
    """Create a small .xyz file with 3 frames."""
    path = tmp_path / "test.xyz"
    frames = [ase.Atoms("H2", positions=[[0, 0, 0], [0, 0, i]]) for i in range(3)]
    ase.io.write(str(path), frames, format="extxyz")
    return path


@pytest.fixture
def extxyz_file(tmp_path: Path) -> Path:
    """Create a small .extxyz file with 3 frames."""
    path = tmp_path / "test.extxyz"
    frames = [ase.Atoms("H2", positions=[[0, 0, 0], [0, 0, i]]) for i in range(3)]
    ase.io.write(str(path), frames, format="extxyz")
    return path


@pytest.fixture
def pdb_file(tmp_path: Path) -> Path:
    """Create a small .pdb file (unregistered in asebytes)."""
    path = tmp_path / "test.pdb"
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    ase.io.write(str(path), atoms, format="proteindatabank")
    return path


def test_open_frames_xyz(xyz_file: Path) -> None:
    """xyz returns a streaming iterator."""
    result = open_frames(str(xyz_file))
    assert isinstance(result, Iterator)
    frames = list(result)
    assert len(frames) == 3
    assert all(isinstance(f, ase.Atoms) for f in frames)


def test_open_frames_extxyz(extxyz_file: Path) -> None:
    """extxyz returns a streaming iterator."""
    result = open_frames(str(extxyz_file))
    assert isinstance(result, Iterator)
    frames = list(result)
    assert len(frames) == 3


def test_open_frames_pdb(pdb_file: Path) -> None:
    """pdb (unregistered in asebytes) falls back to ase.io.iread."""
    result = open_frames(str(pdb_file))
    assert isinstance(result, Iterator)
    frames = list(result)
    assert len(frames) == 1
    assert frames[0].get_chemical_formula() == "H2O"


def test_open_frames_slice(xyz_file: Path) -> None:
    """start/stop/step params slice the frames."""
    result = open_frames(str(xyz_file), start=0, stop=2)
    frames = list(result)
    assert len(frames) == 2


def test_open_frames_step(xyz_file: Path) -> None:
    """step parameter works."""
    result = open_frames(str(xyz_file), step=2)
    frames = list(result)
    assert len(frames) == 2  # frames 0, 2


@pytest.mark.parametrize("suffix", [".h5", ".h5md", ".lmdb"])
def test_asebytes_extensions_recognized(suffix: str) -> None:
    """h5/h5md/lmdb suffixes are in _ASEBYTES_EXTENSIONS."""
    from zndraw.io import _ASEBYTES_EXTENSIONS

    assert suffix in _ASEBYTES_EXTENSIONS
