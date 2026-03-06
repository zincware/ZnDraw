"""Tests for automatic connectivity calculation."""

import ase
import numpy as np

from zndraw.connectivity import add_connectivity


def test_auto_connectivity_small_molecule():
    """Small molecule gets connectivity computed."""
    # Water molecule — O-H bonds should be detected
    atoms = ase.Atoms(
        "OHH",
        positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
    )
    add_connectivity(atoms)

    conn = atoms.info["connectivity"]
    assert len(conn) == 2  # Two O-H bonds


def test_no_bonds_far_apart():
    """Atoms too far apart have no bonds."""
    atoms = ase.Atoms(
        "HH",
        positions=[[0, 0, 0], [100, 0, 0]],
    )
    add_connectivity(atoms)

    conn = atoms.info["connectivity"]
    assert len(conn) == 0


def test_existing_connectivity_overwritten():
    """add_connectivity overwrites existing connectivity info."""
    atoms = ase.Atoms(
        "HH",
        positions=[[0, 0, 0], [0.74, 0, 0]],
    )
    atoms.info["connectivity"] = np.array([[0, 1, 99]], dtype=np.int32)

    add_connectivity(atoms)

    conn = atoms.info["connectivity"]
    # Should be recalculated, bond order should be 1, not 99
    assert len(conn) == 1
    assert conn[0][2] == 1


def test_empty_atoms():
    """Empty atoms don't crash."""
    atoms = ase.Atoms()
    add_connectivity(atoms)
    assert "connectivity" not in atoms.info


def test_single_atom():
    """Single atom has no bonds."""
    atoms = ase.Atoms("H", positions=[[0, 0, 0]])
    add_connectivity(atoms)

    conn = atoms.info["connectivity"]
    assert len(conn) == 0


def test_connectivity_is_numpy_array():
    """Connectivity is stored as numpy int32 array."""
    atoms = ase.Atoms(
        "HH",
        positions=[[0, 0, 0], [0.74, 0, 0]],
    )
    add_connectivity(atoms)

    conn = atoms.info["connectivity"]
    assert isinstance(conn, np.ndarray)
    assert conn.dtype == np.int32


def test_scale_factor_increases_bonds():
    """Larger scale factor detects more bonds."""
    # Place atoms at a distance that's borderline
    atoms = ase.Atoms(
        "HH",
        positions=[[0, 0, 0], [1.5, 0, 0]],
    )

    add_connectivity(atoms, scale=1.0)
    conn_narrow = len(atoms.info["connectivity"])

    add_connectivity(atoms, scale=2.0)
    conn_wide = len(atoms.info["connectivity"])

    assert conn_wide >= conn_narrow


def test_no_self_bonds():
    """No atom is bonded to itself."""
    atoms = ase.Atoms(
        "H3",
        positions=[[0, 0, 0], [0.74, 0, 0], [0.37, 0.64, 0]],
    )
    add_connectivity(atoms)

    conn = atoms.info["connectivity"]
    for row in conn:
        assert row[0] != row[1]


def test_no_duplicate_bonds():
    """Each bond appears exactly once (i < j)."""
    atoms = ase.Atoms(
        "H4",
        positions=[[0, 0, 0], [0.74, 0, 0], [0.37, 0.64, 0], [1.1, 0.64, 0]],
    )
    add_connectivity(atoms)

    conn = atoms.info["connectivity"]
    pairs = {(row[0], row[1]) for row in conn}
    # All i < j
    for i, j in pairs:
        assert i < j
    # No duplicates
    assert len(pairs) == len(conn)
