import pytest
from zndraw import ZnDraw
from zndraw.extensions import selections

CELERY_TIMEOUT = 2

TEST_CASES = [
    (selections.NoneSelection, lambda vis: frozenset(), {}, "NoneSelection"),
    (selections.All, lambda vis: frozenset(range(len(vis[0]))), {}, "All"),
    (selections.Invert, lambda vis: frozenset(range(3, len(vis[0]))), {}, "Invert"),
    (selections.Range, lambda vis: frozenset(range(0, 5, 1)), {"start": 0, "end": 5, "step": 1}, "Range 0-5"),
]

@pytest.mark.parametrize(
    "selection_class, expected_result_func, kwargs, test_id",
    TEST_CASES,
    ids=[case[3] for case in TEST_CASES] # Use the test_id for clearer reporting
)
def test_selections(server, s22, celery_worker, selection_class, expected_result_func, kwargs, test_id):
    """A single, parameterized test for all selection types."""
    vis = ZnDraw(url=server, room="test", user="tester")
    vis.extend(s22)
    
    expected_selection = expected_result_func(vis)
    vis.selection = [0, 1, 2]
    assert vis.selection == frozenset([0, 1, 2])

    selection_instance = selection_class(**kwargs)
    selection_instance.run(vis)
    assert vis.selection == expected_selection

    # --- Test run via Celery ---
    vis.selection = [0, 1, 2]
    assert vis.selection == frozenset([0, 1, 2])

    vis.run(selection_instance)
    vis.socket.sio.sleep(CELERY_TIMEOUT)

    assert vis.selection == expected_selection
