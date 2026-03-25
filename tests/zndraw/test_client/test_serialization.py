"""Behavioral round-trip tests for client serialization helpers."""

import pytest
from molify import smiles2atoms

from zndraw.client.serialization import (
    _estimate_frame_size,
    atoms_to_json_dict,
    json_dict_to_atoms,
)


@pytest.fixture(params=["O", "CCO", "c1ccccc1"], ids=["water", "ethanol", "benzene"])
def atoms(request: pytest.FixtureRequest):
    """Create real ASE Atoms from SMILES."""
    return smiles2atoms(request.param)


def test_round_trip_preserves_atom_count(atoms):
    """Round-trip through JSON dict preserves number of atoms."""
    data = atoms_to_json_dict(atoms)
    restored = json_dict_to_atoms(data)
    assert len(restored) == len(atoms)


def test_round_trip_preserves_positions(atoms):
    """Round-trip through JSON dict preserves atomic positions."""
    data = atoms_to_json_dict(atoms)
    restored = json_dict_to_atoms(data)
    assert restored.positions == pytest.approx(atoms.positions, abs=1e-6)


def test_round_trip_preserves_symbols(atoms):
    """Round-trip through JSON dict preserves chemical symbols."""
    data = atoms_to_json_dict(atoms)
    restored = json_dict_to_atoms(data)
    assert restored.get_chemical_symbols() == atoms.get_chemical_symbols()


def test_atoms_to_json_dict_returns_string_keys(atoms):
    """atoms_to_json_dict returns a dict with string keys (b64-prefixed)."""
    data = atoms_to_json_dict(atoms)
    assert isinstance(data, dict)
    for key in data:
        assert isinstance(key, str)
        assert key.startswith("b64:")


def test_estimate_frame_size_positive_for_nonempty():
    """_estimate_frame_size returns a positive integer for a non-empty frame."""
    atoms = smiles2atoms("O")
    frame = atoms_to_json_dict(atoms)
    size = _estimate_frame_size(frame)
    assert isinstance(size, int)
    assert size > 0


def test_estimate_frame_size_zero_for_empty():
    """_estimate_frame_size returns 0 for an empty frame dict."""
    assert _estimate_frame_size({}) == 0


@pytest.mark.parametrize("threshold", [1, 1000])
def test_connectivity_threshold(threshold):
    """Round-trip works regardless of connectivity_threshold setting."""
    atoms = smiles2atoms("O")  # 3 atoms
    data = atoms_to_json_dict(atoms, connectivity_threshold=threshold)
    restored = json_dict_to_atoms(data)
    assert len(restored) == len(atoms)
    assert restored.get_chemical_symbols() == atoms.get_chemical_symbols()
