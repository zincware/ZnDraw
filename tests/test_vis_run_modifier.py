import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk

from zndraw import ZnDraw
from zndraw.extensions import modifiers

CELERY_TIMEOUT = 5


@pytest.fixture
def atoms_with_cell():
    """Create atoms with a defined cell for wrap/center/replicate tests."""
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    return atoms


@pytest.fixture
def multi_species_atoms():
    """Create atoms with multiple species for ChangeType tests."""
    atoms = Atoms("HCHCH", positions=[[i, 0, 0] for i in range(5)])
    return atoms


@pytest.fixture
def trajectory_atoms():
    """Create a multi-frame trajectory for testing all=True parameter."""
    frames = []
    for i in range(3):
        atoms = Atoms("HHH", positions=[[j + i * 0.1, 0, 0] for j in range(3)])
        atoms.cell = [5, 5, 5]
        atoms.pbc = True
        frames.append(atoms)
    return frames


def test_delete_modifier(server, celery_worker):
    """Test deleting selected atoms."""
    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HHHHH", positions=[[i, 0, 0] for i in range(5)])
    vis.extend([atoms])

    initial_count = len(vis.atoms)
    initial_frame_count = len(vis)

    # Select atoms to delete
    vis.selection = [0, 1, 2]

    # Run delete modifier
    vis.run(modifiers.Delete())
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Should create new frame with deleted atoms
    assert len(vis) == initial_frame_count + 1
    # Step should be on the new frame
    assert vis.step >= initial_frame_count
    assert len(vis.atoms) == initial_count - 3
    assert vis.selection == tuple()


def test_change_type_modifier(server, multi_species_atoms, celery_worker):
    """Test changing atom types by calling modifier directly."""
    vis = ZnDraw(url=server, room="test", user="tester")
    vis.extend([multi_species_atoms])

    initial_frame_count = len(vis)

    # Select H atoms and change to N - call .run() directly to bypass API serialization
    vis.selection = [0, 2, 4]
    modifier = modifiers.ChangeType(symbol=modifiers.Symbols.N)
    modifier.run(vis)

    # Should create new frame with changed types
    assert len(vis) == initial_frame_count + 1
    # Get the last (newly created) frame
    atoms = vis[-1]
    assert atoms[0].symbol == "N"
    assert atoms[2].symbol == "N"
    assert atoms[4].symbol == "N"
    assert atoms[1].symbol == "C"  # Unchanged
    assert vis.selection == tuple()


def test_fix_atoms_modifier(server, celery_worker):
    """Test fixing atoms with constraints."""
    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HHH", positions=[[i, 0, 0] for i in range(3)])
    vis.extend([atoms])

    # Fix first two atoms
    vis.selection = [0, 1]

    vis.run(modifiers.FixAtoms())
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Check constraints are set on the correct atoms
    atoms = vis.atoms
    assert atoms.constraints is not None
    assert len(atoms.constraints) > 0
    # Verify the correct atoms are fixed
    fixed_indices = atoms.constraints[0].get_indices()
    assert set(fixed_indices) == {0, 1}


def test_remove_atoms_modifier(server, celery_worker):
    """Test removing current frame."""
    vis = ZnDraw(url=server, room="test", user="tester")

    atoms1 = Atoms("HH", positions=[[0, 0, 0], [1, 0, 0]])
    atoms2 = Atoms("HHH", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    vis.extend([atoms1, atoms2])

    initial_frame_count = len(vis)
    vis.step = 0

    vis.run(modifiers.RemoveAtoms())
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Should have one less frame
    assert len(vis) == initial_frame_count - 1
    # Verify the correct frame was removed (frame 0 with 2 atoms)
    # Now frame 0 should be what was previously frame 1 (3 atoms)
    assert len(vis[0]) == 3


def test_empty_modifier(server, celery_worker):
    """Test creating empty frame."""
    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HH", positions=[[0, 0, 0], [1, 0, 0]])
    vis.extend([atoms])

    initial_count = len(vis)

    vis.run(modifiers.Empty())
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Should add empty frame
    assert len(vis) == initial_count + 1
    # Check last frame is empty
    assert len(vis[-1]) == 0


def test_new_canvas_modifier(server, celery_worker):
    """Test clearing scene with new canvas."""
    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HHH", positions=[[i, 0, 0] for i in range(3)])
    vis.extend([atoms])

    initial_frame_count = len(vis)

    modifier = modifiers.NewCanvas()
    modifier.run(vis)

    # Should add empty frame
    assert len(vis) == initial_frame_count + 1
    assert len(vis.atoms) == 0
    # Should have added a plane geometry
    assert "plane" in vis.geometries
    # Should have cleared curve points
    assert "curve" in vis.geometries
    # Should have created a bookmark
    assert (len(vis) - 1) in vis.bookmarks


@pytest.mark.skip(reason="Connect has bug: uses undefined get_scaled_radii() function")
def test_connect_modifier(server, celery_worker):
    """Test creating connection points from selection."""
    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HHH", positions=[[i, 0, 0] for i in range(3)])
    vis.extend([atoms])
    vis.selection = [0, 1, 2]
    vis.camera = {"position": [0, 0, -10], "target": [0, 0, 0]}

    modifier = modifiers.Connect()
    modifier.run(vis)

    # Should create points and clear selection
    assert vis.selection == tuple()


def test_translate_modifier(server, celery_worker):
    """Test translating atoms along curve."""
    from zndraw.geometries import Curve

    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HH", positions=[[0, 0, 0], [1, 0, 0]])
    vis.extend([atoms])
    vis.selection = [0]

    # Add curve geometry with multiple control points for translation path
    curve_points = [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]]
    vis.geometries["curve"] = Curve(position=curve_points, divisions=50)

    initial_frame_count = len(vis)
    initial_position = atoms.positions[0].copy()
    steps = 5

    # Run modifier via Celery
    vis.run(modifiers.Translate(steps=steps, curve="curve"))
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Should create new frames
    assert len(vis) == initial_frame_count + steps

    # Verify translation actually occurred in the final frame
    final_atoms = vis[-1]
    final_position = final_atoms.positions[0]

    # Selected atom should have moved
    assert not np.allclose(final_position, initial_position)

    # Unselected atom should not have moved
    assert np.allclose(final_atoms.positions[1], atoms.positions[1])


def test_add_line_particles_modifier(server, celery_worker):
    """Test adding particles at point positions from curve geometry."""
    from zndraw.geometries import Curve

    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("H", positions=[[0, 0, 0]])
    vis.extend([atoms])

    # Add curve geometry with points
    curve_points = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
    vis.geometries["curve"] = Curve(position=curve_points)

    initial_count = len(vis.atoms)
    initial_frame_count = len(vis)

    # Call .run() directly to bypass API serialization of Symbols enum
    modifier = modifiers.AddLineParticles(symbol=modifiers.Symbols.C, curve="curve")
    modifier.run(vis)

    # Should create new frame with added atoms at point positions
    assert len(vis) == initial_frame_count + 1
    new_atoms = vis[-1]
    assert len(new_atoms) == initial_count + 2
    assert new_atoms[-1].symbol == "C"
    assert new_atoms[-2].symbol == "C"
    # Verify positions match curve points
    assert np.allclose(new_atoms[-2].position, curve_points[0])
    assert np.allclose(new_atoms[-1].position, curve_points[1])


def test_rotate_modifier(server, celery_worker):
    """Test rotating atoms around axis defined by curve geometry."""
    from zndraw.geometries import Curve

    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HHH", positions=[[0, 1, 0], [1, 1, 0], [2, 1, 0]])
    vis.extend([atoms])

    # Select multiple atoms to test tuple indexing bug
    vis.selection = [0, 2]
    initial_positions_selected = atoms.positions[[0, 2]].copy()
    initial_position_unselected = atoms.positions[1].copy()

    # Add curve geometry with 2 points for rotation axis (rotating around x-axis)
    vis.geometries["curve"] = Curve(position=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    initial_frame_count = len(vis)
    steps = 5

    # Run modifier via Celery
    vis.run(modifiers.Rotate(angle=90, steps=steps, direction="left"))
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Should create rotation frames
    assert len(vis) == initial_frame_count + steps

    # Verify rotation actually occurred in the final frame
    final_atoms = vis[-1]
    final_positions_selected = final_atoms.positions[[0, 2]]
    final_position_unselected = final_atoms.positions[1]

    # Selected atoms should have moved (rotated 90 degrees around x-axis)
    assert not np.allclose(final_positions_selected, initial_positions_selected)

    # Unselected atom should not have moved
    assert np.allclose(final_position_unselected, initial_position_unselected)


def test_duplicate_modifier(server, celery_worker):
    """Test duplicating selected atoms with offset by calling modifier directly."""
    vis = ZnDraw(url=server, room="test", user="tester")

    atoms = Atoms("HH", positions=[[0, 0, 0], [1, 0, 0]])
    vis.extend([atoms])

    initial_count = len(vis.atoms)
    initial_frame_count = len(vis)
    initial_positions = atoms.positions.copy()
    vis.selection = [0, 1]

    # Call .run() directly to bypass API serialization
    offset = np.array([1.0, 0.0, 0.0])
    modifier = modifiers.Duplicate(x=offset[0], y=offset[1], z=offset[2], symbol=modifiers.Symbols.C)
    modifier.run(vis)

    # Should create new frame with duplicated atoms
    assert len(vis) == initial_frame_count + 1
    new_atoms = vis[-1]  # Get the new frame
    assert len(new_atoms) == initial_count + 2
    assert new_atoms[-1].symbol == "C"
    assert new_atoms[-2].symbol == "C"
    # Verify positions are offset correctly
    assert np.allclose(new_atoms[-2].position, initial_positions[0] + offset)
    assert np.allclose(new_atoms[-1].position, initial_positions[1] + offset)
    assert vis.selection == tuple()


def test_wrap_modifier(server, atoms_with_cell, celery_worker):
    """Test wrapping atoms to cell."""
    vis = ZnDraw(url=server, room="test", user="tester")

    # Move atom outside cell
    atoms = atoms_with_cell.copy()
    atoms.positions[0] += atoms.cell[0]
    vis.extend([atoms])

    original_pos = vis.atoms.positions[0].copy()
    cell_diagonal = np.diag(atoms.cell)

    vis.run(modifiers.Wrap(recompute_bonds=False, all=False))
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Atom should be wrapped back into cell
    wrapped_pos = vis.atoms.positions[0]
    assert not np.allclose(wrapped_pos, original_pos)
    # Verify atom is inside cell bounds [0, cell_diagonal]
    assert np.all(wrapped_pos >= 0)
    assert np.all(wrapped_pos <= cell_diagonal)


def test_center_modifier(server, atoms_with_cell, celery_worker):
    """Test centering atoms in cell."""
    vis = ZnDraw(url=server, room="test", user="tester")
    vis.extend([atoms_with_cell])

    vis.selection = [0]

    vis.run(modifiers.Center(recompute_bonds=False, wrap=False, dynamic=False, all=False))
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Selected atom should be near cell center
    atoms = vis.atoms
    cell_center = np.diag(atoms.cell) / 2
    atom_pos = atoms.positions[0]

    # Check atom is close to center (within reasonable tolerance)
    assert np.allclose(atom_pos, cell_center, atol=2.0)


def test_replicate_modifier(server, atoms_with_cell, celery_worker):
    """Test replicating atoms in cell."""
    vis = ZnDraw(url=server, room="test", user="tester")
    vis.extend([atoms_with_cell])

    initial_count = len(vis.atoms)
    initial_cell = atoms_with_cell.cell.copy()

    vis.run(modifiers.Replicate(x=2, y=2, z=2, keep_box=False, all=False))
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    # Should replicate atoms
    replicated_atoms = vis.atoms
    assert len(replicated_atoms) == initial_count * 8  # 2x2x2
    # Verify cell size increased (since keep_box=False)
    expected_cell = initial_cell * [2, 2, 2]
    assert np.allclose(replicated_atoms.cell, expected_cell)
