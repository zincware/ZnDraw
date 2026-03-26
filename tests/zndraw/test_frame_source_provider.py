"""Unit tests for FrameSourceRead provider."""

import ase
import msgpack
import numpy as np
import pytest

from zndraw.providers.frame_source import FrameSourceRead


@pytest.fixture
def water():
    """Simple water molecule for testing."""
    return ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])


@pytest.fixture
def source(water):
    """Minimal FrameSource (list of Atoms)."""
    return [water, ase.Atoms("He"), ase.Atoms("Ar")]


def test_read_returns_bytes(source):
    """read() returns bytes (msgpack-encoded RawFrame)."""
    provider = FrameSourceRead(index=0)
    result = provider.read(source)
    assert isinstance(result, bytes)


def test_read_roundtrips_to_raw_frame(source):
    """Unpacking the result yields a dict[bytes, bytes]."""
    provider = FrameSourceRead(index=0)
    raw = msgpack.unpackb(provider.read(source), raw=True)
    assert isinstance(raw, dict)
    assert all(isinstance(k, bytes) for k in raw)
    assert all(isinstance(v, bytes) for v in raw.values())


def test_read_contains_positions(source):
    """Encoded frame contains arrays.positions key."""
    provider = FrameSourceRead(index=0)
    raw = msgpack.unpackb(provider.read(source), raw=True)
    assert b"arrays.positions" in raw


def test_read_adds_colors(source):
    """Colors array is added on the fly."""
    provider = FrameSourceRead(index=0)
    raw = msgpack.unpackb(provider.read(source), raw=True)
    assert b"arrays.colors" in raw


def test_read_adds_radii(source):
    """Radii array is added on the fly."""
    provider = FrameSourceRead(index=0)
    raw = msgpack.unpackb(provider.read(source), raw=True)
    assert b"arrays.radii" in raw


def test_read_adds_connectivity(water):
    """Connectivity is added for small molecules."""
    provider = FrameSourceRead(index=0)
    raw = msgpack.unpackb(provider.read([water]), raw=True)
    assert b"info.connectivity" in raw


def test_read_respects_index(source):
    """Different indices read different frames."""
    r0 = msgpack.unpackb(FrameSourceRead(index=0).read(source), raw=True)
    r1 = msgpack.unpackb(FrameSourceRead(index=1).read(source), raw=True)
    # He (index=1) has 1 atom, H2O (index=0) has 3
    assert r0[b"arrays.positions"] != r1[b"arrays.positions"]


def test_read_preserves_existing_colors():
    """If atoms already have colors, they are not overwritten."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms.set_array("colors", np.array(["#ff0000", "#00ff00"], dtype="U7"))
    raw = msgpack.unpackb(FrameSourceRead(index=0).read([atoms]), raw=True)
    assert b"arrays.colors" in raw


@pytest.mark.parametrize("index", [0, 1, 2])
def test_read_all_indices(source, index):
    """Every valid index produces a valid msgpack result."""
    result = FrameSourceRead(index=index).read(source)
    raw = msgpack.unpackb(result, raw=True)
    assert b"arrays.numbers" in raw


def test_category_and_content_type():
    """Class vars are set correctly."""
    assert FrameSourceRead.category == "frames"
    assert FrameSourceRead.content_type == "application/x-msgpack"
