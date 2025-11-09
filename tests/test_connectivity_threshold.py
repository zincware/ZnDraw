import numpy as np
import pytest

import ase.build
from zndraw.zndraw import ZnDraw


def test_connectivity_added_on_append_small_structure(server):
    """Test that connectivity is automatically added for small structures on append."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Create a small structure (< 1000 atoms)
    atoms = ase.build.molecule("H2O")
    assert len(atoms) < 1000
    assert "connectivity" not in atoms.info

    vis.append(atoms)

    # Retrieve and check connectivity was added
    retrieved = vis[0]
    assert "connectivity" in retrieved.info
    assert isinstance(retrieved.info["connectivity"], np.ndarray)
    assert retrieved.info["connectivity"].dtype == np.int32


def test_connectivity_not_added_on_append_large_structure(server):
    """Test that connectivity is NOT automatically added for large structures on append."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Create a large structure (>= 1000 atoms)
    atoms = ase.build.bulk("Cu", cubic=True) * (10, 10, 10)  # 4000 atoms
    assert len(atoms) >= 1000
    assert "connectivity" not in atoms.info

    vis.append(atoms)

    # Retrieve and check connectivity was NOT added
    retrieved = vis[0]
    assert "connectivity" not in retrieved.info


def test_connectivity_preserved_when_already_present(server):
    """Test that existing connectivity is preserved and not recomputed."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Create atoms with custom connectivity
    atoms = ase.build.molecule("H2O")
    custom_connectivity = np.array([[0, 1, 1], [1, 2, 1]], dtype=np.int32)
    atoms.info["connectivity"] = custom_connectivity

    vis.append(atoms)

    # Retrieve and check custom connectivity was preserved
    retrieved = vis[0]
    assert "connectivity" in retrieved.info
    np.testing.assert_array_equal(retrieved.info["connectivity"], custom_connectivity)


def test_connectivity_added_on_extend(server):
    """Test that connectivity is automatically added for small structures on extend."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Create multiple small structures
    atoms_list = [ase.build.molecule("H2O") for _ in range(5)]
    for atoms in atoms_list:
        assert "connectivity" not in atoms.info

    vis.extend(atoms_list)

    # Check all frames have connectivity
    for i in range(len(atoms_list)):
        retrieved = vis[i]
        assert "connectivity" in retrieved.info
        assert isinstance(retrieved.info["connectivity"], np.ndarray)


def test_connectivity_added_on_insert(server):
    """Test that connectivity is automatically added for small structures on insert."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Add initial atom
    vis.append(ase.build.molecule("H2O"))

    # Insert a new atom
    atoms = ase.build.molecule("NH3")
    assert "connectivity" not in atoms.info

    vis.insert(0, atoms)

    # Check inserted frame has connectivity
    retrieved = vis[0]
    assert "connectivity" in retrieved.info


def test_connectivity_added_on_setitem_single(server):
    """Test that connectivity is automatically added for small structures on single setitem."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Add initial atom
    vis.append(ase.build.molecule("H2O"))

    # Replace with new atom
    atoms = ase.build.molecule("NH3")
    assert "connectivity" not in atoms.info

    vis[0] = atoms

    # Check replaced frame has connectivity
    retrieved = vis[0]
    assert "connectivity" in retrieved.info


def test_connectivity_added_on_setitem_list(server):
    """Test that connectivity is automatically added for small structures on list setitem."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Add initial atoms
    vis.extend([ase.build.molecule("H2O") for _ in range(3)])

    # Replace with new atoms
    atoms_list = [ase.build.molecule("NH3") for _ in range(3)]
    for atoms in atoms_list:
        assert "connectivity" not in atoms.info

    vis[0:3] = atoms_list

    # Check all replaced frames have connectivity
    for i in range(3):
        retrieved = vis[i]
        assert "connectivity" in retrieved.info


def test_custom_connectivity_threshold(server):
    """Test that custom connectivity_threshold parameter works."""
    # Create ZnDraw with custom threshold of 5 atoms
    vis = ZnDraw(url=server, room="testroom", user="testuser", connectivity_threshold=5)

    # Create a 3-atom structure (< 5, should get connectivity)
    small_atoms = ase.build.molecule("H2O")
    assert len(small_atoms) == 3
    assert "connectivity" not in small_atoms.info

    vis.append(small_atoms)
    retrieved_small = vis[0]
    assert "connectivity" in retrieved_small.info

    # Create a 10-atom structure (>= 5, should NOT get connectivity)
    large_atoms = ase.build.molecule("C6H6")  # benzene
    assert len(large_atoms) >= 5
    assert "connectivity" not in large_atoms.info

    vis.append(large_atoms)
    retrieved_large = vis[1]
    assert "connectivity" not in retrieved_large.info


def test_connectivity_threshold_zero_disables_auto_connectivity(server):
    """Test that connectivity_threshold=0 disables automatic connectivity."""
    vis = ZnDraw(url=server, room="testroom", user="testuser", connectivity_threshold=0)

    # Even small structures should not get connectivity
    atoms = ase.build.molecule("H2O")
    assert "connectivity" not in atoms.info

    vis.append(atoms)

    retrieved = vis[0]
    assert "connectivity" not in retrieved.info


def test_connectivity_mixed_sizes_in_extend(server):
    """Test that connectivity is selectively added based on size in extend."""
    vis = ZnDraw(url=server, room="testroom", user="testuser", connectivity_threshold=10)

    # Mix of small and large structures
    small_atoms = ase.build.molecule("H2O")  # 3 atoms
    large_atoms = ase.build.bulk("Cu", cubic=True) * (2, 2, 2)  # 32 atoms

    atoms_list = [small_atoms, large_atoms, small_atoms]
    vis.extend(atoms_list)

    # Small structures should have connectivity
    assert "connectivity" in vis[0].info
    assert "connectivity" in vis[2].info

    # Large structure should NOT have connectivity
    assert "connectivity" not in vis[1].info
