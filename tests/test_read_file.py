import ase
import pytest

from zndraw import ZnDraw
from zndraw.app import tasks
from zndraw.app.tasks import calculate_adaptive_resolution

# It is assumed that you have fixtures named 'server', 's22_xyz', 's22_h5',
# and 's22' defined in a conftest.py file.


@pytest.mark.parametrize(
    "file_fixture_name", ["s22_xyz", "s22_h5", "s22_db", "s22_json_db"]
)
def test_worker_read_file(server, s22, file_fixture_name, request):
    """
    Test that reading a complete file (.xyz, .h5, .db, or .json) results in the correct
    number of structures in the ZnDraw visualization.

    This test combines:
    - test_worker_read_file_xyz
    - test_worker_read_file_h5
    - test_worker_read_file_db
    - test_worker_read_file_json_db
    """
    # 'request.getfixturevalue' allows us to dynamically get a fixture by its name.
    file = request.getfixturevalue(file_fixture_name)

    tasks.read_file(file=file, room="test_room", server_url=server)
    vis = ZnDraw(room="test_room", url=server, user="tester")

    assert len(vis) == len(s22)


@pytest.mark.parametrize(
    "file_fixture_name, kwargs, expected_slice",
    [
        # --- Positive Slicing ---
        ("s22_xyz", {"start": 5}, slice(5, None, None)),
        ("s22_h5", {"start": 5}, slice(5, None, None)),
        ("s22_db", {"start": 5}, slice(5, None, None)),
        ("s22_xyz", {"stop": 5}, slice(None, 5, None)),
        ("s22_h5", {"stop": 5}, slice(None, 5, None)),
        ("s22_db", {"stop": 5}, slice(None, 5, None)),
        ("s22_xyz", {"step": 2}, slice(None, None, 2)),
        ("s22_h5", {"step": 2}, slice(None, None, 2)),
        ("s22_db", {"step": 2}, slice(None, None, 2)),
        ("s22_xyz", {"start": 2, "stop": 8, "step": 2}, slice(2, 8, 2)),
        ("s22_h5", {"start": 2, "stop": 8, "step": 2}, slice(2, 8, 2)),
        ("s22_db", {"start": 2, "stop": 8, "step": 2}, slice(2, 8, 2)),
        # --- Negative Slicing (Not supported for databases) ---
        ("s22_xyz", {"start": -5}, slice(-5, None, None)),  # Last 5 items
        ("s22_h5", {"start": -5}, slice(-5, None, None)),
        ("s22_xyz", {"stop": -5}, slice(None, -5, None)),  # All but the last 5
        ("s22_h5", {"stop": -5}, slice(None, -5, None)),
        ("s22_xyz", {"step": -1}, slice(None, None, -1)),  # Reversed
        # ("s22_h5", {"step": -1}, slice(None, None, -1)),
        ("s22_xyz", {"step": -2}, slice(None, None, -2)),  # Reversed with a step
        # ("s22_h5", {"step": -2}, slice(None, None, -2)),
    ],
)
def test_worker_read_file_slicing(
    server, s22, file_fixture_name, kwargs, expected_slice, request
):
    """
    Test that reading a file with slicing parameters (start, stop, step)
    works correctly for both positive and negative values.

    This test combines and extends:
    - test_worker_read_file_xyz_start
    - test_worker_read_file_xyz_stop
    - test_worker_read_file_xyz_step
    - test_worker_read_file_h5_start
    - test_worker_read_file_h5_stop
    - test_worker_read_file_h5_step
    """
    file = request.getfixturevalue(file_fixture_name)

    # Unpack the keyword arguments for start, stop, step
    tasks.read_file(file=file, room="test_room", server_url=server, **kwargs)
    vis = ZnDraw(room="test_room", url=server, user="tester")

    expected_data = s22[expected_slice]
    assert len(vis) == len(expected_data), f"Failed with parameters: {kwargs}"

    # Check that the content of each structure matches
    for atoms_vis, atoms_expected in zip(vis, expected_data):
        assert atoms_vis == atoms_expected


def test_worker_read_file_db_empty(server, tmp_path):
    """Test reading an empty ASE database."""
    import ase.db

    # Create empty database
    db_path = tmp_path / "empty.db"
    db = ase.db.connect(str(db_path))
    # Don't write any structures

    # Should not raise an error, just return early
    tasks.read_file(file=str(db_path), room="test_room_empty", server_url=server)
    vis = ZnDraw(room="test_room_empty", url=server, user="tester")

    assert len(vis) == 0


def test_worker_read_file_db_out_of_bounds_start(server, s22_db):
    """Test reading a database with start index exceeding database size."""
    import ase.db

    # Get actual database size
    db = ase.db.connect(s22_db)
    n_rows = db.count()

    # Try to read with start index beyond database size
    with pytest.raises(ValueError, match="Start row .* exceeds database size"):
        tasks.read_file(
            file=s22_db,
            room="test_room_oob",
            server_url=server,
            start=n_rows + 5,  # Beyond database
        )


def test_worker_read_file_db_stop_exceeds_size(server, s22, s22_db):
    """Test reading a database with stop index exceeding database size.

    This should warn but still work, using the end of the database.
    """
    import ase.db

    # Get actual database size
    db = ase.db.connect(s22_db)
    n_rows = db.count()

    # Read with stop beyond database size - should work but warn
    tasks.read_file(
        file=s22_db,
        room="test_room_stop_exceed",
        server_url=server,
        stop=n_rows + 10,  # Beyond database
    )
    vis = ZnDraw(room="test_room_stop_exceed", url=server, user="tester")

    # Should have loaded all structures (stop clamped to n_rows)
    assert len(vis) == len(s22)


@pytest.mark.parametrize("invalid_step", [0, -1, -2])
def test_worker_read_file_db_invalid_step(server, s22_db, invalid_step):
    """Test reading a database with invalid step values (0 or negative)."""
    with pytest.raises(ValueError, match="Step must be a positive integer"):
        tasks.read_file(
            file=s22_db,
            room="test_room_invalid_step",
            server_url=server,
            step=invalid_step,
        )


def test_worker_read_file_db_partial_range(server, s22, s22_db):
    """Test reading a specific range from the middle of a database."""
    # Read structures 3-7 (0-indexed: rows with IDs 4-7 in database)
    tasks.read_file(
        file=s22_db,
        room="test_room_partial",
        server_url=server,
        start=3,
        stop=7,
    )
    vis = ZnDraw(room="test_room_partial", url=server, user="tester")

    expected_data = s22[3:7]
    assert len(vis) == len(expected_data)

    # Check content matches
    for atoms_vis, atoms_expected in zip(vis, expected_data):
        assert atoms_vis == atoms_expected


@pytest.mark.parametrize(
    "num_particles, expected_resolution",
    [
        (100, 16),  # Small system
        (500, 16),  # Still small
        (1000, 14),  # Boundary case
        (2500, 14),  # Medium-small
        (5000, 12),  # Medium
        (10000, 10),  # Large
        (25000, 8),  # Very large
        (100000, 8),  # Huge
        (500000, 6),  # Massive
    ],
)
def test_calculate_adaptive_resolution(num_particles, expected_resolution):
    """Test that adaptive resolution calculation produces expected values."""
    resolution = calculate_adaptive_resolution(num_particles)
    assert resolution == expected_resolution
    assert 6 <= resolution <= 16, "Resolution should be between 6 and 16"


def test_adaptive_resolution_applied_large_system(server, tmp_path):
    """Test that adaptive resolution is applied when loading a large system."""
    # Create a file with a large number of atoms
    atoms = ase.Atoms("H" * 15000, positions=[[i, 0, 0] for i in range(15000)])
    xyz_file = tmp_path / "large_system.xyz"
    ase.io.write(xyz_file, atoms)

    # Load the file
    tasks.read_file(file=str(xyz_file), room="test_large_system", server_url=server)
    vis = ZnDraw(room="test_large_system", url=server, user="tester")

    # Check that resolution was reduced
    particles_geometry = vis.geometries["particles"]
    assert particles_geometry.resolution < 16, (
        "Resolution should be reduced for large systems"
    )
    assert particles_geometry.resolution == 10, (
        "Expected resolution of 10 for 15000 particles"
    )


def test_adaptive_resolution_not_applied_small_system(server, s22_xyz):
    """Test that adaptive resolution is NOT applied for small systems."""
    # Load a small system (s22 has small molecules)
    tasks.read_file(file=s22_xyz, room="test_small_system", server_url=server)
    vis = ZnDraw(room="test_small_system", url=server, user="tester")

    # Check that resolution remains at default
    particles_geometry = vis.geometries["particles"]
    assert particles_geometry.resolution == 16, (
        "Resolution should remain 16 for small systems"
    )
