"""Tests for constraint serialization and storage.

This test suite verifies that ASE constraints can be properly serialized,
stored in Zarr, and reconstructed. All constraints that inherit from
ase.constraints.FixConstraint should be supported via ASE's todict()
and dict2constraint() methods.
"""
import ase
import ase.constraints
import numpy as np
import numpy.testing as npt
import pytest
import zarr
from zarr.storage import MemoryStore

from zndraw.storage import ZarrStorageSequence
from zndraw.utils import atoms_from_dict, atoms_to_dict


def test_fixatoms_serialization():
    """Test FixAtoms constraint serialization and deserialization."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixAtoms(indices=[0, 2])
    atoms.set_constraint(constraint)

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" in data
    assert len(data["constraints"]) == 1
    assert data["constraints"][0]["name"] == "FixAtoms"
    assert data["constraints"][0]["kwargs"]["indices"] == [0, 2]

    # Deserialize
    atoms2 = atoms_from_dict(data)
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    npt.assert_array_equal(atoms2.constraints[0].index, [0, 2])


def test_fixbondlength_serialization():
    """Test FixBondLength constraint serialization and deserialization."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    constraint = ase.constraints.FixBondLength(0, 1)
    atoms.set_constraint(constraint)

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" in data
    assert len(data["constraints"]) == 1
    # FixBondLength creates a FixBondLengths constraint
    assert data["constraints"][0]["name"] == "FixBondLengths"

    # Deserialize
    atoms2 = atoms_from_dict(data)
    assert len(atoms2.constraints) == 1
    # Restored as FixBondLengths
    assert atoms2.constraints[0].__class__.__name__ == "FixBondLengths"


def test_fixbondlengths_serialization():
    """Test FixBondLengths constraint serialization and deserialization."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixBondLengths(pairs=[[0, 1], [0, 2]])
    atoms.set_constraint(constraint)

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" in data
    assert len(data["constraints"]) == 1
    assert data["constraints"][0]["name"] == "FixBondLengths"

    # Deserialize
    atoms2 = atoms_from_dict(data)
    assert len(atoms2.constraints) == 1
    assert atoms2.constraints[0].__class__.__name__ == "FixBondLengths"
    npt.assert_array_equal(atoms2.constraints[0].pairs, [[0, 1], [0, 2]])


def test_fixedline_serialization():
    """Test FixedLine constraint serialization and deserialization."""
    atoms = ase.Atoms("H", positions=[[0, 0, 0]])
    constraint = ase.constraints.FixedLine(0, [1, 0, 0])
    atoms.set_constraint(constraint)

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" in data
    assert len(data["constraints"]) == 1
    assert data["constraints"][0]["name"] == "FixedLine"

    # Deserialize
    atoms2 = atoms_from_dict(data)
    assert len(atoms2.constraints) == 1
    assert atoms2.constraints[0].__class__.__name__ == "FixedLine"
    npt.assert_array_equal(atoms2.constraints[0].index, [0])
    npt.assert_array_almost_equal(atoms2.constraints[0].dir, [1, 0, 0])


def test_fixedplane_serialization():
    """Test FixedPlane constraint serialization and deserialization."""
    atoms = ase.Atoms("H", positions=[[0, 0, 0]])
    constraint = ase.constraints.FixedPlane(0, [0, 0, 1])
    atoms.set_constraint(constraint)

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" in data
    assert len(data["constraints"]) == 1
    assert data["constraints"][0]["name"] == "FixedPlane"

    # Deserialize
    atoms2 = atoms_from_dict(data)
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

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" in data
    assert len(data["constraints"]) == 2
    assert data["constraints"][0]["name"] == "FixAtoms"
    assert data["constraints"][1]["name"] == "FixBondLengths"

    # Deserialize
    atoms2 = atoms_from_dict(data)
    assert len(atoms2.constraints) == 2
    assert atoms2.constraints[0].__class__.__name__ == "FixAtoms"
    assert atoms2.constraints[1].__class__.__name__ == "FixBondLengths"


def test_no_constraints_serialization():
    """Test that atoms without constraints don't have a constraints key."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" not in data

    # Deserialize should work fine
    atoms2 = atoms_from_dict(data)
    assert len(atoms2.constraints) == 0


def test_constraints_with_zarr_storage():
    """Integration test: constraints should persist through zarr storage."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixAtoms(indices=[0, 2])
    atoms.set_constraint(constraint)

    # Convert to dict
    data = atoms_to_dict(atoms)

    # Store in Zarr
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)
    store.append(data)

    # Retrieve from Zarr
    retrieved = store[0]

    # Verify constraints are preserved
    assert "constraints" in retrieved
    assert len(retrieved["constraints"]) == 1
    assert retrieved["constraints"][0]["name"] == "FixAtoms"
    assert retrieved["constraints"][0]["kwargs"]["indices"] == [0, 2]

    # Convert back to atoms
    atoms2 = atoms_from_dict(retrieved)
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    npt.assert_array_equal(atoms2.constraints[0].index, [0, 2])


def test_constraints_with_variable_sized_atoms():
    """Test that constraints work with variable-sized atoms in zarr storage."""
    # Create multiple frames with different sizes and constraints
    atoms1 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms1.set_constraint(ase.constraints.FixAtoms(indices=[0]))

    atoms2 = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms2.set_constraint(ase.constraints.FixAtoms(indices=[0, 2]))

    atoms3 = ase.Atoms("CH4", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]])
    atoms3.set_constraint(ase.constraints.FixAtoms(indices=[0]))

    # Convert to dicts
    data1 = atoms_to_dict(atoms1)
    data2 = atoms_to_dict(atoms2)
    data3 = atoms_to_dict(atoms3)

    # Store in Zarr
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)
    store.extend([data1, data2, data3])

    # Retrieve and verify
    retrieved1 = store[0]
    retrieved2 = store[1]
    retrieved3 = store[2]

    # Check frame 1
    assert retrieved1["constraints"][0]["kwargs"]["indices"] == [0]
    atoms1_restored = atoms_from_dict(retrieved1)
    npt.assert_array_equal(atoms1_restored.constraints[0].index, [0])

    # Check frame 2
    assert retrieved2["constraints"][0]["kwargs"]["indices"] == [0, 2]
    atoms2_restored = atoms_from_dict(retrieved2)
    npt.assert_array_equal(atoms2_restored.constraints[0].index, [0, 2])

    # Check frame 3
    assert retrieved3["constraints"][0]["kwargs"]["indices"] == [0]
    atoms3_restored = atoms_from_dict(retrieved3)
    npt.assert_array_equal(atoms3_restored.constraints[0].index, [0])


def test_constraints_with_calculator():
    """Test that constraints work alongside calculator data."""
    from ase.calculators.singlepoint import SinglePointCalculator

    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms.set_constraint(ase.constraints.FixAtoms(indices=[0]))

    # Add calculator
    calc = SinglePointCalculator(atoms, energy=1.5, forces=[[0.1, 0, 0], [0.2, 0, 0]])
    atoms.calc = calc

    # Serialize
    data = atoms_to_dict(atoms)
    assert "constraints" in data
    assert "calc.energy" in data
    assert "calc.forces" in data

    # Deserialize
    atoms2 = atoms_from_dict(data)
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    assert atoms2.calc.results["energy"] == 1.5
    npt.assert_array_almost_equal(atoms2.calc.results["forces"], [[0.1, 0, 0], [0.2, 0, 0]])


def test_constraint_roundtrip_fixatoms():
    """Test FixAtoms constraint serialization roundtrip."""
    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixAtoms(indices=[0, 1])
    atoms.set_constraint(constraint)

    # Roundtrip: atoms -> dict -> atoms
    data = atoms_to_dict(atoms)
    atoms2 = atoms_from_dict(data)

    # Verify constraint was restored
    assert len(atoms2.constraints) == 1
    assert isinstance(atoms2.constraints[0], ase.constraints.FixAtoms)
    npt.assert_array_equal(atoms2.constraints[0].index, [0, 1])


def test_constraint_roundtrip_fixbondlength():
    """Test FixBondLength constraint serialization roundtrip."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    constraint = ase.constraints.FixBondLength(0, 1)
    atoms.set_constraint(constraint)

    # Roundtrip: atoms -> dict -> atoms
    data = atoms_to_dict(atoms)
    atoms2 = atoms_from_dict(data)

    # Verify constraint was restored (as FixBondLengths)
    assert len(atoms2.constraints) == 1
    assert atoms2.constraints[0].__class__.__name__ == "FixBondLengths"


def test_constraints_preserved_across_multiple_frames():
    """Test that different constraints are preserved independently across frames."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # Frame 1: FixAtoms
    atoms1 = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms1.set_constraint(ase.constraints.FixAtoms(indices=[0]))
    store.append(atoms_to_dict(atoms1))

    # Frame 2: FixBondLength
    atoms2 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms2.set_constraint(ase.constraints.FixBondLength(0, 1))
    store.append(atoms_to_dict(atoms2))

    # Frame 3: Multiple constraints
    atoms3 = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms3.set_constraint(
        [
            ase.constraints.FixAtoms(indices=[0]),
            ase.constraints.FixBondLength(1, 2),
        ]
    )
    store.append(atoms_to_dict(atoms3))

    # Retrieve and verify each frame independently
    retrieved1 = atoms_from_dict(store[0])
    assert len(retrieved1.constraints) == 1
    assert isinstance(retrieved1.constraints[0], ase.constraints.FixAtoms)

    retrieved2 = atoms_from_dict(store[1])
    assert len(retrieved2.constraints) == 1
    assert retrieved2.constraints[0].__class__.__name__ == "FixBondLengths"

    retrieved3 = atoms_from_dict(store[2])
    assert len(retrieved3.constraints) == 2
    assert isinstance(retrieved3.constraints[0], ase.constraints.FixAtoms)
    assert retrieved3.constraints[1].__class__.__name__ == "FixBondLengths"
