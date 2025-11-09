import numpy as np

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


def test_empty_atoms_append(server):
    """Test that empty Atoms objects can be appended without errors."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Create empty Atoms
    empty_atoms = ase.Atoms()
    assert len(empty_atoms) == 0
    assert "connectivity" not in empty_atoms.info

    # Should not raise an error
    vis.append(empty_atoms)

    # Retrieve and verify
    retrieved = vis[0]
    assert len(retrieved) == 0
    assert "connectivity" not in retrieved.info
    # Colors and radii should still be added
    assert "colors" in retrieved.arrays
    assert "radii" in retrieved.arrays


def test_empty_atoms_extend(server):
    """Test that empty Atoms objects can be extended without errors."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Create list with empty and non-empty Atoms
    empty_atoms = ase.Atoms()
    water = ase.build.molecule("H2O")

    atoms_list = [empty_atoms, water, empty_atoms]
    vis.extend(atoms_list)

    # Check empty atoms don't have connectivity
    assert len(vis[0]) == 0
    assert "connectivity" not in vis[0].info

    # Check non-empty atoms do have connectivity
    assert len(vis[1]) == 3
    assert "connectivity" in vis[1].info

    # Check second empty atoms
    assert len(vis[2]) == 0
    assert "connectivity" not in vis[2].info


def test_empty_atoms_insert(server):
    """Test that empty Atoms objects can be inserted without errors."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Add initial atoms
    vis.append(ase.build.molecule("H2O"))

    # Insert empty atoms
    empty_atoms = ase.Atoms()
    vis.insert(0, empty_atoms)

    # Verify
    assert len(vis[0]) == 0
    assert "connectivity" not in vis[0].info


def test_empty_atoms_setitem(server):
    """Test that empty Atoms objects can replace existing frames without errors."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Add initial atoms
    vis.append(ase.build.molecule("H2O"))

    # Replace with empty atoms
    empty_atoms = ase.Atoms()
    vis[0] = empty_atoms

    # Verify
    retrieved = vis[0]
    assert len(retrieved) == 0
    assert "connectivity" not in retrieved.info


def test_empty_atoms_with_custom_info(server):
    """Test that empty Atoms with custom info are handled correctly."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")

    # Create empty Atoms with custom info
    empty_atoms = ase.Atoms()
    empty_atoms.info["custom_key"] = "custom_value"

    vis.append(empty_atoms)

    # Verify custom info is preserved
    retrieved = vis[0]
    assert retrieved.info["custom_key"] == "custom_value"
    assert "connectivity" not in retrieved.info
