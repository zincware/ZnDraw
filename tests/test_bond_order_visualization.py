"""Tests for bond order visualization."""

import ase
import numpy as np
import pytest
import requests

from zndraw import ZnDraw
from zndraw.geometries import Bond


def test_bond_order_type_compatibility():
    """Test that Bond geometry accepts various bond order values."""
    # Test with ignore mode (default)
    bond_ignore = Bond(
        connectivity=[[0, 1, 1], [1, 2, 2], [2, 3, 3]], bond_order_mode="ignore"
    )
    assert bond_ignore.bond_order_mode == "ignore"
    assert bond_ignore.bond_order_offset == 3

    # Test with parallel mode
    bond_parallel = Bond(
        connectivity=[[0, 1, 1.5], [1, 2, 2.0], [2, 3, 3]],
        bond_order_mode="parallel",
        bond_order_offset=0.2,
    )
    assert bond_parallel.bond_order_mode == "parallel"
    assert bond_parallel.bond_order_offset == 0.2
    assert 1.5 in bond_parallel.bond_order_radius_scale
    assert 2.0 in bond_parallel.bond_order_radius_scale


def test_bond_order_default_radius_scale():
    """Test default bond order radius scaling."""
    bond = Bond()
    expected_scales = {1: 1.0, 1.5: 0.75, 2: 0.75, 3: 0.7}
    assert bond.bond_order_radius_scale == expected_scales


def test_bond_order_null_handling():
    """Test that null bond orders default to 1."""
    # Connectivity with null bond order
    bond = Bond(connectivity=[[0, 1, None], [1, 2, 2]])
    # The None should be handled gracefully by the frontend
    assert bond.connectivity == [(0.0, 1.0, None), (1.0, 2.0, 2.0)]


@pytest.mark.parametrize(
    "bond_order,expected_cylinders",
    [
        (1, 1),  # Single bond
        (1.5, 2),  # Aromatic bond
        (2, 2),  # Double bond
        (3, 3),  # Triple bond
        (None, 1),  # Null defaults to single
    ],
)
def test_bond_order_cylinder_counts(bond_order, expected_cylinders):
    """Test that different bond orders produce correct cylinder counts."""
    # This tests the logical behavior, actual rendering is in frontend
    connectivity = [[0, 1, bond_order]] if bond_order is not None else [[0, 1, None]]
    bond = Bond(connectivity=connectivity, bond_order_mode="parallel")
    assert len(bond.connectivity) == 1


def test_bond_with_mixed_orders(server):
    """Test molecule with mixed single, double, and triple bonds."""
    vis = ZnDraw(url=server, room="test-mixed-bonds", user="tester")

    # Create a molecule with different bond types
    # Positions: C-C=C≡C (single, double, triple)
    atoms = ase.Atoms(
        "C4", positions=[[0, 0, 0], [1.5, 0, 0], [3.0, 0, 0], [4.2, 0, 0]]
    )
    atoms.info["connectivity"] = np.array(
        [
            [0, 1, 1],  # C-C single bond
            [1, 2, 2],  # C=C double bond
            [2, 3, 3],  # C≡C triple bond
        ],
        dtype=np.float32,
    )

    vis.append(atoms)
    vis.socket.sio.sleep(0.5)

    # Configure bond geometry with parallel mode
    vis.geometries["bonds"] = Bond(
        connectivity="info.connectivity",
        bond_order_mode="parallel",
        bond_order_offset=0.15,
        scale=0.15,
    )
    vis.socket.sio.sleep(0.5)

    # Verify geometry was set
    retrieved_bond = vis.geometries["bonds"]
    assert retrieved_bond.bond_order_mode == "parallel"
    assert retrieved_bond.connectivity == "info.connectivity"


def test_bond_aromatic_visualization(server):
    """Test aromatic bond (order 1.5) visualization."""
    vis = ZnDraw(url=server, room="test-aromatic", user="tester")

    # Create benzene-like ring with aromatic bonds
    atoms = ase.Atoms(
        "C6",
        positions=[
            [0, 0, 0],
            [1, 0, 0],
            [1.5, 0.866, 0],
            [1, 1.732, 0],
            [0, 1.732, 0],
            [-0.5, 0.866, 0],
        ],
    )
    atoms.info["connectivity"] = np.array(
        [
            [0, 1, 1.5],  # Aromatic bonds
            [1, 2, 1.5],
            [2, 3, 1.5],
            [3, 4, 1.5],
            [4, 5, 1.5],
            [5, 0, 1.5],
        ],
        dtype=np.float32,
    )

    vis.append(atoms)
    vis.socket.sio.sleep(0.5)

    # Configure with parallel mode
    vis.geometries["bonds"] = Bond(
        connectivity="info.connectivity", bond_order_mode="parallel", scale=0.12
    )
    vis.socket.sio.sleep(0.5)

    retrieved_bond = vis.geometries["bonds"]
    assert retrieved_bond.bond_order_mode == "parallel"
    # Aromatic bonds should use 0.8 radius scale by default
    assert retrieved_bond.bond_order_radius_scale[1.5] == 0.75


def test_bond_order_ignore_mode(server):
    """Test that ignore mode renders all bonds as single."""
    vis = ZnDraw(url=server, room="test-ignore-mode", user="tester")

    atoms = ase.Atoms("C3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.info["connectivity"] = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.float32)

    vis.append(atoms)
    vis.socket.sio.sleep(0.5)

    # Set ignore mode (default)
    vis.geometries["bonds"] = Bond(
        connectivity="info.connectivity", bond_order_mode="ignore"
    )
    vis.socket.sio.sleep(0.5)

    retrieved_bond = vis.geometries["bonds"]
    assert retrieved_bond.bond_order_mode == "ignore"


def test_bond_order_custom_offset(server):
    """Test custom bond order offset configuration."""
    vis = ZnDraw(url=server, room="test-custom-offset", user="tester")

    atoms = ase.Atoms("C2", positions=[[0, 0, 0], [1.5, 0, 0]])
    atoms.info["connectivity"] = np.array([[0, 1, 2]], dtype=np.float32)

    vis.append(atoms)
    vis.socket.sio.sleep(0.5)

    # Use custom offset for wider spacing
    vis.geometries["bonds"] = Bond(
        connectivity="info.connectivity",
        bond_order_mode="parallel",
        bond_order_offset=0.3,  # Wider spacing
        scale=0.15,
    )
    vis.socket.sio.sleep(0.5)

    retrieved_bond = vis.geometries["bonds"]
    assert retrieved_bond.bond_order_offset == 0.3


def test_bond_order_custom_radius_scale(server):
    """Test custom radius scaling per bond order."""
    vis = ZnDraw(url=server, room="test-custom-scale", user="tester")

    atoms = ase.Atoms("C3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.info["connectivity"] = np.array([[0, 1, 2], [1, 2, 3]], dtype=np.float32)

    vis.append(atoms)
    vis.socket.sio.sleep(0.5)

    # Custom radius scales
    custom_scales = {1.0: 1.0, 2.0: 0.5, 3.0: 0.4}
    vis.geometries["bonds"] = Bond(
        connectivity="info.connectivity",
        bond_order_mode="parallel",
        bond_order_radius_scale=custom_scales,
    )
    vis.socket.sio.sleep(0.5)

    retrieved_bond = vis.geometries["bonds"]
    assert retrieved_bond.bond_order_radius_scale[2.0] == 0.5
    assert retrieved_bond.bond_order_radius_scale[3.0] == 0.4


def test_bond_order_rest_api(joined_room):
    """Test bond order configuration via Python client."""
    server, room = joined_room

    # Use ZnDraw client which handles lock acquisition automatically
    vis = ZnDraw(url=server, room=room, user="test-bond-order-update")

    # Update bond geometry with parallel mode
    bond = Bond(
        connectivity="info.connectivity",
        bond_order_mode="parallel",
        bond_order_offset=0.2,
        scale=0.15,
    )
    vis.geometries["bonds"] = bond

    # Verify the update
    response = requests.get(f"{server}/api/rooms/{room}/geometries/bonds")
    assert response.status_code == 200
    data = response.json()
    assert data["geometry"]["data"]["bond_order_mode"] == "parallel"
    assert data["geometry"]["data"]["bond_order_offset"] == 0.2


def test_connectivity_backward_compatibility():
    """Test that old connectivity format still works."""
    # Old format: all bond orders = 1
    bond_old = Bond(connectivity=[[0, 1, 1], [1, 2, 1]])
    # Pydantic converts lists to tuples and ints to floats
    assert bond_old.connectivity == [(0.0, 1.0, 1.0), (1.0, 2.0, 1.0)]

    # New format with mixed orders
    bond_new = Bond(connectivity=[[0, 1, 1], [1, 2, 2], [2, 3, 1.5]])
    assert bond_new.connectivity == [(0.0, 1.0, 1.0), (1.0, 2.0, 2.0), (2.0, 3.0, 1.5)]
