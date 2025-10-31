"""Test for connectivity shape mismatch when mixing numpy and list connectivity."""

import ase
import numpy as np
import pytest
import zarr
from zarr.storage import MemoryStore

from zndraw.connectivity import add_connectivity
from zndraw.storage import extend_zarr
from zndraw.utils import atoms_to_dict


def test_connectivity_list_vs_numpy_should_work():
    """Test that appending atoms with list connectivity to store with numpy connectivity works.

    This tests the fix for the issue where:
    1. First atoms have connectivity as numpy array with shape (n_bonds, 3)
    2. Second atoms have connectivity as list (from rdkit2ase)
    3. Storage should automatically convert the list to numpy array
    """
    # Create first atoms with numpy connectivity (using add_connectivity)
    atoms1 = ase.Atoms(
        "CH4", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]]
    )
    add_connectivity(atoms1)  # This creates numpy array connectivity
    dict1 = atoms_to_dict(atoms1)

    # Verify first atoms have numpy array connectivity
    assert isinstance(dict1["info.connectivity"], np.ndarray)
    assert dict1["info.connectivity"].shape == (4, 3)  # 4 bonds, 3 values each

    # Create storage and add first atoms
    store = MemoryStore()
    root = zarr.open_group(store, mode="w")
    extend_zarr(root, [dict1])

    # Create second atoms with list connectivity (simulating rdkit2ase)
    atoms2 = ase.Atoms("C2OH6", positions=np.random.rand(9, 3))
    atoms2.info["connectivity"] = [
        (0, 1, 1.0),
        (1, 2, 1.0),
        (0, 3, 1.0),
        (0, 4, 1.0),
        (0, 5, 1.0),
        (1, 6, 1.0),
        (1, 7, 1.0),
        (2, 8, 1.0),
    ]  # List, not numpy array
    dict2 = atoms_to_dict(atoms2)

    # Verify second atoms have list connectivity
    assert isinstance(dict2["info.connectivity"], list)

    # This should work - list should be automatically converted to numpy array
    extend_zarr(root, [dict2])

    # Verify both frames are stored
    assert root["info.connectivity"].shape[0] == 2
    # Both should have valid connectivity data
    assert root["info.connectivity"][0].shape[1] == 3  # First frame: 4 bonds x 3 values
    assert (
        root["info.connectivity"][1].shape[1] == 3
    )  # Second frame: 8 bonds x 3 values


def test_connectivity_numpy_to_numpy_works():
    """Test that appending atoms with numpy connectivity to store with numpy connectivity works."""
    # Create first atoms with numpy connectivity
    atoms1 = ase.Atoms(
        "CH4", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]]
    )
    add_connectivity(atoms1)
    dict1 = atoms_to_dict(atoms1)

    # Create storage and add first atoms
    store = MemoryStore()
    root = zarr.open_group(store, mode="w")
    extend_zarr(root, [dict1])

    # Create second atoms with numpy connectivity (different number of bonds)
    atoms2 = ase.Atoms("C2H6", positions=np.random.rand(8, 3))
    add_connectivity(atoms2)
    dict2 = atoms_to_dict(atoms2)

    # Verify both have numpy arrays
    assert isinstance(dict1["info.connectivity"], np.ndarray)
    assert isinstance(dict2["info.connectivity"], np.ndarray)

    # This should work - variable-length arrays are supported
    extend_zarr(root, [dict2])

    # Verify both frames are stored
    assert root["info.connectivity"].shape[0] == 2
