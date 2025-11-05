import pytest
from ase import Atoms

from zndraw import ZnDraw
from zndraw.extensions import selections

CELERY_TIMEOUT = 5

TEST_CASES = [
    (selections.NoneSelection, lambda vis: tuple(), {}, "NoneSelection"),
    (selections.All, lambda vis: tuple(range(len(vis[0]))), {}, "All"),
    (selections.Invert, lambda vis: tuple(range(3, len(vis[0]))), {}, "Invert"),
    (
        selections.Range,
        lambda vis: tuple(range(0, 5, 1)),
        {"start": 0, "end": 5, "step": 1},
        "Range 0-5",
    ),
    (selections.UpdateSelection, lambda vis: (0, 1, 2), {}, "UpdateSelection"),
]


@pytest.mark.parametrize(
    "selection_class, expected_result_func, kwargs, test_id",
    TEST_CASES,
    ids=[case[3] for case in TEST_CASES],  # Use the test_id for clearer reporting
)
def test_selections(
    server, s22, celery_worker, selection_class, expected_result_func, kwargs, test_id
):
    """A single, parameterized test for all selection types."""
    vis = ZnDraw(url=server, room="test", user="tester")
    vis.extend(s22)

    expected_selection = expected_result_func(vis)
    vis.selection = [0, 1, 2]
    assert vis.selection == (0, 1, 2)

    selection_instance = selection_class(**kwargs)
    selection_instance.run(vis)
    assert vis.selection == expected_selection

    # --- Test run via Celery ---
    vis.selection = [0, 1, 2]
    assert vis.selection == (0, 1, 2)

    vis.run(selection_instance)
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    assert vis.selection == expected_selection


def test_random_selection(server, s22, celery_worker):
    """Test random selection with seeded randomness for reproducibility."""
    vis = ZnDraw(url=server, room="test", user="tester")
    vis.extend(s22)

    count = 5
    expected_selection = (0, 1, 2, 5, 7)  # Deterministic with seed=42

    # Test direct run
    selection_instance = selections.Random(count=count, seed=42)
    selection_instance.run(vis)
    assert set(vis.selection) == set(expected_selection)

    # Test run via Celery
    selection_instance = selections.Random(count=count, seed=42)
    vis.run(selection_instance)
    vis.socket.sio.sleep(CELERY_TIMEOUT)
    assert set(vis.selection) == set(expected_selection)


def test_identical_species_selection(server, celery_worker):
    """Test selecting all atoms of the same species as selected atoms."""
    vis = ZnDraw(url=server, room="test", user="tester")

    # Create structure: H-C-H-C-H (indices 0,1,2,3,4)
    atoms = Atoms("HCHCH", positions=[[i, 0, 0] for i in range(5)])
    vis.extend([atoms])

    # Select one H atom (index 0), should select all H atoms (0, 2, 4)
    vis.selection = [0]
    expected_selection = (0, 2, 4)

    # Test direct run
    selection_instance = selections.IdenticalSpecies()
    selection_instance.run(vis)
    assert set(vis.selection) == set(expected_selection)

    # Test run via Celery - select one C atom, should get all C atoms (1, 3)
    vis.selection = [1]
    expected_selection_c = (1, 3)

    vis.run(selections.IdenticalSpecies())
    vis.socket.sio.sleep(CELERY_TIMEOUT)
    assert set(vis.selection) == set(expected_selection_c)


def test_connected_particles_selection(server, celery_worker):
    """Test selecting connected particles in a graph with disconnected components."""
    vis = ZnDraw(url=server, room="test", user="tester")

    # Create two disconnected components:
    # Component 1: atoms 0-1-2 (chain)
    # Component 2: atoms 3-4 (chain)
    atoms = Atoms("HHHHH", positions=[[i, 0, 0] for i in range(5)])
    connectivity = [
        [0, 1, 1.0],  # 0 connected to 1
        [1, 2, 1.0],  # 1 connected to 2
        [3, 4, 1.0],  # 3 connected to 4 (separate component)
    ]
    atoms.info["connectivity"] = connectivity
    vis.extend([atoms])

    # Select atom 0, should get entire component 1 (0, 1, 2)
    vis.selection = [0]
    expected_selection = (0, 1, 2)

    # Test direct run
    selection_instance = selections.ConnectedParticles()
    selection_instance.run(vis)
    assert set(vis.selection) == set(expected_selection)

    # Test run via Celery - select atom 3, should get component 2 (3, 4)
    vis.selection = [3]
    expected_selection_2 = (3, 4)

    vis.run(selections.ConnectedParticles())
    vis.socket.sio.sleep(CELERY_TIMEOUT)
    assert set(vis.selection) == set(expected_selection_2)


def test_neighbour_selection(server, celery_worker):
    """Test selecting nth order neighbours in a linear chain."""
    vis = ZnDraw(url=server, room="test", user="tester")

    # Create linear chain: 0-1-2-3-4
    atoms = Atoms("HHHHH", positions=[[i, 0, 0] for i in range(5)])
    connectivity = [
        [0, 1, 1.0],
        [1, 2, 1.0],
        [2, 3, 1.0],
        [3, 4, 1.0],
    ]
    atoms.info["connectivity"] = connectivity
    vis.extend([atoms])

    # Test 1st order neighbour: select atom 2, should get 1, 2, 3
    vis.selection = [2]
    expected_selection_1st = (1, 2, 3)

    selection_instance = selections.Neighbour(order=1)
    selection_instance.run(vis)
    assert set(vis.selection) == set(expected_selection_1st)

    # Test 2nd order neighbour via Celery: select atom 2, should get 0, 1, 2, 3, 4
    vis.selection = [2]
    expected_selection_2nd = (0, 1, 2, 3, 4)

    vis.run(selections.Neighbour(order=2))
    vis.socket.sio.sleep(CELERY_TIMEOUT)
    assert set(vis.selection) == set(expected_selection_2nd)
