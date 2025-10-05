import pytest

from zndraw import ZnDraw
from zndraw.app import tasks

# It is assumed that you have fixtures named 'server', 's22_xyz', 's22_h5',
# and 's22' defined in a conftest.py file.


@pytest.mark.parametrize("file_fixture_name", ["s22_xyz", "s22_h5"])
def test_worker_read_file(server, s22, file_fixture_name, request):
    """
    Test that reading a complete file (.xyz or .h5) results in the correct
    number of structures in the ZnDraw visualization.

    This test combines:
    - test_worker_read_file_xyz
    - test_worker_read_file_h5
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
        ("s22_xyz", {"stop": 5}, slice(None, 5, None)),
        ("s22_h5", {"stop": 5}, slice(None, 5, None)),
        ("s22_xyz", {"step": 2}, slice(None, None, 2)),
        ("s22_h5", {"step": 2}, slice(None, None, 2)),
        ("s22_xyz", {"start": 2, "stop": 8, "step": 2}, slice(2, 8, 2)),
        ("s22_h5", {"start": 2, "stop": 8, "step": 2}, slice(2, 8, 2)),
        # --- Negative Slicing (New Tests) ---
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
