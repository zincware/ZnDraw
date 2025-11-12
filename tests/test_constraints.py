"""Tests for constraint serialization and storage.

This test suite verifies that ASE constraints can be properly serialized,
stored in ASEBytes storage, and reconstructed. All constraints that inherit from
ase.constraints.FixConstraint should be supported via ASE's todict()
and dict2constraint() methods.
"""

import ase
import ase.constraints
import numpy as np
import numpy.testing as npt
import pytest

from asebytes import encode, decode
from zndraw.storage import ASEBytesStorageBackend


def test_fixatoms_serialization():
    """Test FixAtoms constraint serialization and deserialization."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixAtoms(indices=[0, 2])
    atoms.set_constraint(constraint)

    # Serialize with asebytes (natively supports constraints)
    data = encode(atoms)

    # Deserialize
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    npt.assert_array_equal(atoms2.constraints[0].index, [0, 2])


def test_fixbondlength_serialization():
    """Test FixBondLength constraint serialization and deserialization."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    constraint = ase.constraints.FixBondLength(0, 1)
    atoms.set_constraint(constraint)

    # Serialize with asebytes (natively supports constraints)
    data = encode(atoms)

    # Deserialize
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 1
    # Restored as FixBondLengths
    assert atoms2.constraints[0].__class__.__name__ == "FixBondLengths"


def test_fixbondlengths_serialization():
    """Test FixBondLengths constraint serialization and deserialization."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixBondLengths(pairs=[[0, 1], [0, 2]])
    atoms.set_constraint(constraint)

    # Serialize with asebytes (natively supports constraints)
    data = encode(atoms)

    # Deserialize
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 1
    assert atoms2.constraints[0].__class__.__name__ == "FixBondLengths"
    npt.assert_array_equal(atoms2.constraints[0].pairs, [[0, 1], [0, 2]])


def test_fixedline_serialization():
    """Test FixedLine constraint serialization and deserialization."""
    atoms = ase.Atoms("H", positions=[[0, 0, 0]])
    constraint = ase.constraints.FixedLine(0, [1, 0, 0])
    atoms.set_constraint(constraint)

    # Serialize with asebytes (natively supports constraints)
    data = encode(atoms)

    # Deserialize
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 1
    assert atoms2.constraints[0].__class__.__name__ == "FixedLine"
    npt.assert_array_equal(atoms2.constraints[0].index, [0])
    npt.assert_array_almost_equal(atoms2.constraints[0].dir, [1, 0, 0])


def test_fixedplane_serialization():
    """Test FixedPlane constraint serialization and deserialization."""
    atoms = ase.Atoms("H", positions=[[0, 0, 0]])
    constraint = ase.constraints.FixedPlane(0, [0, 0, 1])
    atoms.set_constraint(constraint)

    # Serialize with asebytes (natively supports constraints)
    data = encode(atoms)

    # Deserialize
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 1
    assert atoms2.constraints[0].__class__.__name__ == "FixedPlane"
    npt.assert_array_equal(atoms2.constraints[0].index, [0])
    npt.assert_array_almost_equal(atoms2.constraints[0].dir, [0, 0, 1])


def test_multiple_constraints_serialization():
    """Test serialization of multiple constraints on the same atoms object."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraints = [
        ase.constraints.FixAtoms(indices=[0]),
        ase.constraints.FixBondLength(1, 2),
    ]
    atoms.set_constraint(constraints)

    # Serialize with asebytes (natively supports constraints)
    data = encode(atoms)

    # Deserialize
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 2
    assert atoms2.constraints[0].__class__.__name__ == "FixAtoms"
    assert atoms2.constraints[1].__class__.__name__ == "FixBondLengths"


def test_no_constraints_serialization():
    """Test that atoms without constraints don't have a constraints key."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])

    # Serialize with asebytes
    data = encode(atoms)

    # Deserialize should work fine
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 0


def test_constraints_with_storage(tmp_path):
    """Integration test: constraints should persist through ASEBytes storage."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixAtoms(indices=[0, 2])
    atoms.set_constraint(constraint)

    # Serialize with asebytes
    data = encode(atoms)

    # Store in ASEBytes
    db_path = str(tmp_path / "test.lmdb")
    store = ASEBytesStorageBackend(db_path, map_size=1_000_000_000)
    store.append(data)

    # Retrieve from storage
    retrieved = store[0]

    # Convert back to atoms
    atoms2 = decode(retrieved)
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    npt.assert_array_equal(atoms2.constraints[0].index, [0, 2])


def test_constraints_with_variable_sized_atoms(tmp_path):
    """Test that constraints work with variable-sized atoms in ASEBytes storage."""
    # Create multiple frames with different sizes and constraints
    atoms1 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms1.set_constraint(ase.constraints.FixAtoms(indices=[0]))

    atoms2 = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms2.set_constraint(ase.constraints.FixAtoms(indices=[0, 2]))

    atoms3 = ase.Atoms(
        "CH4", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    )
    atoms3.set_constraint(ase.constraints.FixAtoms(indices=[0]))

    # Serialize with asebytes
    data1 = encode(atoms1)
    data2 = encode(atoms2)
    data3 = encode(atoms3)

    # Store in ASEBytes
    db_path = str(tmp_path / "test.lmdb")
    store = ASEBytesStorageBackend(db_path, map_size=1_000_000_000)
    store.extend([data1, data2, data3])

    # Retrieve and verify
    retrieved1 = store[0]
    retrieved2 = store[1]
    retrieved3 = store[2]

    # Check frame 1
    atoms1_restored = decode(retrieved1)
    npt.assert_array_equal(atoms1_restored.constraints[0].index, [0])

    # Check frame 2
    atoms2_restored = decode(retrieved2)
    npt.assert_array_equal(atoms2_restored.constraints[0].index, [0, 2])

    # Check frame 3
    atoms3_restored = decode(retrieved3)
    npt.assert_array_equal(atoms3_restored.constraints[0].index, [0])


def test_constraints_with_calculator():
    """Test that constraints work alongside calculator data."""
    from ase.calculators.singlepoint import SinglePointCalculator

    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms.set_constraint(ase.constraints.FixAtoms(indices=[0]))

    # Add calculator
    calc = SinglePointCalculator(atoms, energy=1.5, forces=[[0.1, 0, 0], [0.2, 0, 0]])
    atoms.calc = calc

    # Serialize with asebytes (supports both constraints and calculator)
    data = encode(atoms)

    # Deserialize
    atoms2 = decode(data)
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    assert atoms2.calc.results["energy"] == 1.5
    npt.assert_array_almost_equal(
        atoms2.calc.results["forces"], [[0.1, 0, 0], [0.2, 0, 0]]
    )


def test_constraint_roundtrip_fixatoms():
    """Test FixAtoms constraint serialization roundtrip."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixAtoms(indices=[0, 1])
    atoms.set_constraint(constraint)

    # Roundtrip: atoms -> bytes -> atoms
    data = encode(atoms)
    atoms2 = decode(data)

    # Verify constraint was restored
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    npt.assert_array_equal(atoms2.constraints[0].index, [0, 1])


def test_constraint_roundtrip_fixbondlength():
    """Test FixBondLength constraint serialization roundtrip."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    constraint = ase.constraints.FixBondLength(0, 1)
    atoms.set_constraint(constraint)

    # Roundtrip: atoms -> bytes -> atoms
    data = encode(atoms)
    atoms2 = decode(data)

    # Verify constraint was restored (as FixBondLengths)
    assert len(atoms2.constraints) == 1
    assert atoms2.constraints[0].__class__.__name__ == "FixBondLengths"


def test_constraints_preserved_across_multiple_frames(tmp_path):
    """Test that different constraints are preserved independently across frames."""
    db_path = str(tmp_path / "test.lmdb")
    store = ASEBytesStorageBackend(db_path, map_size=1_000_000_000)

    # Frame 1: FixAtoms
    atoms1 = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms1.set_constraint(ase.constraints.FixAtoms(indices=[0]))
    store.append(encode(atoms1))

    # Frame 2: FixBondLength
    atoms2 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms2.set_constraint(ase.constraints.FixBondLength(0, 1))
    store.append(encode(atoms2))

    # Frame 3: Multiple constraints
    atoms3 = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms3.set_constraint(
        [
            ase.constraints.FixAtoms(indices=[0]),
            ase.constraints.FixBondLength(1, 2),
        ]
    )
    store.append(encode(atoms3))

    # Retrieve and verify each frame independently
    retrieved1 = decode(store[0])
    assert len(retrieved1.constraints) == 1
    assert isinstance(retrieved1.constraints[0], ase.constraints.FixAtoms)

    retrieved2 = decode(store[1])
    assert len(retrieved2.constraints) == 1
    assert retrieved2.constraints[0].__class__.__name__ == "FixBondLengths"

    retrieved3 = decode(store[2])
    assert len(retrieved3.constraints) == 2
    assert isinstance(retrieved3.constraints[0], ase.constraints.FixAtoms)
    assert retrieved3.constraints[1].__class__.__name__ == "FixBondLengths"
