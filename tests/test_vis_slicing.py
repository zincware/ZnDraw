import numpy as np
import pytest
from ase.build import molecule

from zndraw import ZnDraw


def test_getitem_slice(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    # Basic slicing
    slice1 = vis[5:15]
    assert len(slice1) == 10
    for i in range(10):
        assert slice1[i] == s22[5 + i]
    slice2 = vis[:10]
    assert len(slice2) == 10
    for i in range(10):
        assert slice2[i] == s22[i]
    slice3 = vis[10:]
    assert len(slice3) == len(s22) - 10
    for i in range(len(s22) - 10):
        assert slice3[i] == s22[10 + i]
    slice4 = vis[-10:]
    assert len(slice4) == 10
    for i in range(10):
        assert slice4[i] == s22[len(s22) - 10 + i]
    slice5 = vis[:-10]
    assert len(slice5) == len(s22) - 10
    for i in range(len(s22) - 10):
        assert slice5[i] == s22[i]

    # Slicing with step
    slice6 = vis[0:20:2]
    assert len(slice6) == 10
    for i in range(10):
        assert slice6[i] == s22[0 + i * 2]
    slice7 = vis[1:20:3]
    assert len(slice7) == 7
    for i in range(7):
        assert slice7[i] == s22[1 + i * 3]
    slice8 = vis[::5]
    assert len(slice8) == (len(s22) + 4) // 5
    for i in range(len(slice8)):
        assert slice8[i] == s22[i * 5]
    slice9 = vis[::-1]
    assert len(slice9) == len(s22)
    for i in range(len(s22)):
        assert slice9[i] == s22[len(s22) - 1 - i]
    slice10 = vis[-1:-11:-1]
    assert len(slice10) == 10
    for i in range(10):
        assert slice10[i] == s22[len(s22) - 1 - i]


def test_getitem_slice_invalid(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    # Invalid slice types
    with pytest.raises(TypeError):
        _ = vis["invalid"]
    with pytest.raises(TypeError):
        _ = vis[5.5]
    with pytest.raises(TypeError):
        _ = vis[5:15:1.5]

    # Out of range slices should not raise errors
    slice1 = vis[1000:1010]
    assert len(slice1) == 0
    slice2 = vis[-1000:-990]
    assert len(slice2) == 0

    # Slices with step 0 should raise ValueError
    with pytest.raises(ValueError):
        _ = vis[::0]
    with pytest.raises(ValueError):
        _ = vis[5:15:0]


def test_getitem_index(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    # Valid indices
    for i in range(len(s22)):
        assert vis[i] == s22[i]
    for i in range(-len(s22), 0):
        assert vis[i] == s22[len(s22) + i]
    for i in range(len(s22)):
        assert vis[i] == s22[np.array(i)]

    # Invalid indices
    with pytest.raises(IndexError):
        _ = vis[len(s22)]
    with pytest.raises(IndexError):
        _ = vis[-len(s22) - 1]


def test_getitem_list(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    slice = vis[[1, 2, 3]]
    assert len(slice) == 3
    for i in range(3):
        assert slice[i] == s22[i + 1]

    slice = vis[[0, -1, len(s22) - 2]]
    assert len(slice) == 3
    assert slice[0] == s22[0]
    assert slice[1] == s22[len(s22) - 1]
    assert slice[2] == s22[len(s22) - 2]

    slice = vis[np.array([0, 2, -1])]
    assert len(slice) == 3
    assert slice[0] == s22[0]
    assert slice[1] == s22[2]
    assert slice[2] == s22[len(s22) - 1]

    with pytest.raises(IndexError):
        _ = vis[[0, len(s22)]]
    with pytest.raises(IndexError):
        _ = vis[[-len(s22) - 1, 0]]
    with pytest.raises(TypeError):
        _ = vis[[1, 2.5, 3]]
    with pytest.raises(TypeError):
        _ = vis[[1, "2", 3]]
    with pytest.raises(TypeError):
        _ = vis[[1, [2], 3]]


def test_delitem_slice(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    assert len(vis) == len(s22)
    del vis[:]
    assert len(vis) == 0

    vis.extend(s22)
    assert len(vis) == len(s22)

    del vis[5:15]
    assert len(vis) == len(s22) - 10
    for i in range(5):
        assert vis[i] == s22[i]
    for i in range(5, len(vis)):
        assert vis[i] == s22[i + 10]

    # reset
    del vis[:]
    vis.extend(s22)

    del vis[:-10]
    assert len(vis) == 10
    for i in range(10):
        assert vis[i] == s22[len(s22) - 10 + i]

    # reset
    del vis[:]
    vis.extend(s22)

    del vis[::2]
    assert len(vis) == (len(s22) + 1) // 2
    for i in range(len(vis)):
        assert vis[i] == s22[i * 2 + 1]


def test_delitem_index(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)
    assert len(vis) == len(s22)

    del vis[0]
    assert len(vis) == len(s22) - 1
    for i in range(len(vis)):
        assert vis[i] == s22[i + 1]

    # reset
    del vis[:]
    vis.extend(s22)

    del vis[-1]
    assert len(vis) == len(s22) - 1
    for i in range(len(vis)):
        assert vis[i] == s22[i]

    # reset
    del vis[:]
    vis.extend(s22)


def test_invalid_delitem_index(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)
    assert len(vis) == len(s22)

    with pytest.raises(IndexError):
        del vis[len(s22)]
    with pytest.raises(IndexError):
        del vis[-len(s22) - 1]
    with pytest.raises(TypeError):
        del vis[5.5]
    with pytest.raises(TypeError):
        del vis["invalid"]

    assert len(vis) == len(s22)


def test_delitem_list(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    del vis[[1, 2, 3]]
    assert len(vis) == len(s22) - 3
    assert vis[0] == s22[0]
    for i in range(1, len(vis)):
        assert vis[i] == s22[i + 3]

    # reset
    del vis[:]
    vis.extend(s22)

    del vis[[0, -1]]
    assert len(vis) == len(s22) - 2
    assert vis[0] == s22[1]
    for i in range(1, len(vis)):
        assert vis[i] == s22[i + 1]

    # reset
    del vis[:]
    vis.extend(s22)

    del vis[np.array([0, -1])]
    assert len(vis) == len(s22) - 2
    assert vis[0] == s22[1]
    for i in range(1, len(vis)):
        assert vis[i] == s22[i + 1]

    # reset
    del vis[:]
    vis.extend(s22)

    with pytest.raises(IndexError):
        del vis[[0, len(vis)]]
    with pytest.raises(IndexError):
        del vis[[-len(vis) - 1, 0]]
    with pytest.raises(TypeError):
        del vis[[1, 2.5, 3]]
    with pytest.raises(TypeError):
        del vis[[1, "2", 3]]
    with pytest.raises(TypeError):
        del vis[[1, [2], 3]]


def test_setitem_slice(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    new_atoms = molecule("H2O")

    vis[5:15] = [new_atoms.copy() for _ in range(10)]
    assert len(vis) == len(s22)
    for i in range(5):
        assert vis[i] == s22[i]
    for i in range(5, 15):
        assert vis[i] == new_atoms
    for i in range(15, len(vis)):
        assert vis[i] == s22[i]

    # reset
    del vis[:]
    vis.extend(s22)

    vis[:10] = [new_atoms.copy() for _ in range(10)]
    assert len(vis) == len(s22)
    for i in range(10):
        assert vis[i] == new_atoms
    for i in range(10, len(vis)):
        assert vis[i] == s22[i]

    # reset
    del vis[:]
    vis.extend(s22)

    vis[-10:] = [new_atoms.copy() for _ in range(10)]
    assert len(vis) == len(s22)
    for i in range(len(vis) - 10):
        assert vis[i] == s22[i]
    for i in range(len(vis) - 10, len(vis)):
        assert vis[i] == new_atoms

    # reset
    del vis[:]
    vis.extend(s22)

    vis[::2] = [new_atoms.copy() for _ in range((len(s22) + 1) // 2)]
    assert len(vis) == len(s22)
    for i in range((len(s22) + 1) // 2):
        assert vis[i * 2] == new_atoms
    for i in range(len(s22) // 2):
        assert vis[i * 2 + 1] == s22[i * 2 + 1]


def test_setitem_index(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    new_atoms = molecule("H2O")

    vis[0] = new_atoms
    assert len(vis) == len(s22)
    assert vis[0] == new_atoms
    for i in range(1, len(vis)):
        assert vis[i] == s22[i]

    # reset
    del vis[:]
    vis.extend(s22)

    vis[-1] = new_atoms
    assert len(vis) == len(s22)
    for i in range(len(vis) - 1):
        assert vis[i] == s22[i]
    assert vis[-1] == new_atoms

    # reset
    del vis[:]
    vis.extend(s22)

    vis[np.array(0)] = new_atoms
    assert len(vis) == len(s22)
    assert vis[0] == new_atoms
    for i in range(1, len(vis)):
        assert vis[i] == s22[i]


def test_invalid_setitem_index(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)
    assert len(vis) == len(s22)

    new_atoms = molecule("H2O")

    with pytest.raises(IndexError):
        vis[len(s22)] = new_atoms
    with pytest.raises(IndexError):
        vis[-len(s22) - 1] = new_atoms
    with pytest.raises(TypeError):
        vis[5.5] = new_atoms
    with pytest.raises(TypeError):
        vis["invalid"] = new_atoms
    with pytest.raises(TypeError):
        vis[0] = "not an atoms"

    assert len(vis) == len(s22)
    for i in range(len(vis)):
        assert vis[i] == s22[i]


def test_setitem_list(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    new_atoms = molecule("H2O")

    vis[[1, 2, 3]] = [new_atoms.copy() for _ in range(3)]
    assert len(vis) == len(s22)
    assert vis[0] == s22[0]
    for i in range(1, 4):
        assert vis[i] == new_atoms
    for i in range(4, len(vis)):
        assert vis[i] == s22[i]

    # reset
    del vis[:]
    vis.extend(s22)

    vis[[0, -1]] = [new_atoms.copy() for _ in range(2)]
    assert len(vis) == len(s22)
    assert vis[0] == new_atoms
    for i in range(1, len(vis) - 1):
        assert vis[i] == s22[i]
    assert vis[-1] == new_atoms

    # reset
    del vis[:]
    vis.extend(s22)

    vis[np.array([0, -1])] = [new_atoms.copy() for _ in range(2)]
    assert len(vis) == len(s22)
    assert vis[0] == new_atoms
    for i in range(1, len(vis) - 1):
        assert vis[i] == s22[i]
    assert vis[-1] == new_atoms

    # reset
    del vis[:]
    vis.extend(s22)

    with pytest.raises(IndexError):
        vis[[0, len(vis)]] = [new_atoms.copy() for _ in range(2)]
    with pytest.raises(IndexError):
        vis[[-len(vis) - 1, 0]] = [new_atoms.copy() for _ in range(2)]
    with pytest.raises(TypeError):
        vis[[1, 2.5, 3]] = [new_atoms.copy() for _ in range(3)]
    with pytest.raises(TypeError):
        vis[[1, "2", 3]] = [new_atoms.copy() for _ in range(3)]
    with pytest.raises(TypeError):
        vis[[1, [2], 3]] = [new_atoms.copy() for _ in range(3)]
    with pytest.raises(TypeError):
        vis[[1, 2, 3]] = "not a list"
    with pytest.raises(ValueError):
        vis[[1, 2, 3]] = [new_atoms.copy() for _ in range(2)]


def test_setitem_slice_unequal_length(server, s22):
    """Test assigning a list of a different size to a simple slice."""
    new_atoms = molecule("H2O")

    # --- Test Case 1: Grow the list ---
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)
    original_len = len(vis)

    # Replace 2 elements (at index 5 and 6) with 5 new elements
    vis[5:7] = [new_atoms.copy() for _ in range(5)]

    assert len(vis) == original_len - 2 + 5
    # Check elements before the slice
    for i in range(5):
        assert vis[i] == s22[i]
    # Check the newly inserted elements
    for i in range(5, 10):
        assert vis[i] == new_atoms
    # Check elements after the slice
    for i in range(10, len(vis)):
        assert vis[i] == s22[i - 3]  # original index was i-5+2 = i-3

    # --- Test Case 2: Shrink the list ---
    # reset
    del vis[:]
    vis.extend(s22)
    original_len = len(vis)

    # Replace 10 elements (from index 5 to 14) with 2 new elements
    vis[5:15] = [new_atoms.copy() for _ in range(2)]

    assert len(vis) == original_len - 10 + 2
    # Check elements before the slice
    for i in range(5):
        assert vis[i] == s22[i]
    # Check the newly inserted elements
    for i in range(5, 7):
        assert vis[i] == new_atoms
    # Check elements after the slice
    for i in range(7, len(vis)):
        assert vis[i] == s22[i + 8]  # original index was i-2+10 = i+8


def test_setitem_extended_slice_invalid_length(server, s22):
    """Test ValueError when assigning a list of incorrect size to an extended slice."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)

    new_atoms = molecule("H2O")

    # The slice vis[::2] has length (len(s22) + 1) // 2
    slice_len = (len(s22) + 1) // 2

    # This should fail because the list has length 2, but the slice is longer
    with pytest.raises(ValueError):
        vis[::2] = [new_atoms.copy() for _ in range(2)]

    # This should fail because the list is too long
    with pytest.raises(ValueError):
        vis[::2] = [new_atoms.copy() for _ in range(slice_len + 5)]

    # This should fail for a different step
    with pytest.raises(ValueError):
        vis[5:15:3] = [new_atoms.copy()]  # slice has len 4, assigned list has len 1


def test_setitem_slice_insertion(server, s22):
    """Test inserting items using an empty slice assignment, e.g., vis[5:5] = ..."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)
    original_len = len(vis)
    new_atoms = molecule("H2O")

    # --- Insert in the middle ---
    insert_data = [new_atoms.copy() for _ in range(3)]
    vis[5:5] = insert_data

    assert len(vis) == original_len + 3
    # Check elements before the insertion point
    for i in range(5):
        assert vis[i] == s22[i]
    # Check the newly inserted elements
    for i in range(3):
        assert vis[5 + i] == new_atoms
    # Check elements after the insertion point
    for i in range(5, original_len):
        assert vis[i + 3] == s22[i]

    # --- Insert at the beginning ---
    # reset
    del vis[:]
    vis.extend(s22)
    vis[0:0] = insert_data
    assert len(vis) == original_len + 3
    for i in range(3):
        assert vis[i] == new_atoms
    for i in range(original_len):
        assert vis[i + 3] == s22[i]

    # --- Insert at the end ---
    # reset
    del vis[:]
    vis.extend(s22)
    vis[len(vis) : len(vis)] = insert_data
    assert len(vis) == original_len + 3
    for i in range(original_len):
        assert vis[i] == s22[i]
    for i in range(3):
        assert vis[original_len + i] == new_atoms


def test_delitem_slice_edge_cases(server, s22):
    """Test edge cases for __delitem__ like negative steps and empty slices."""
    # --- Test deletion with a negative step ---
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)
    original_len = len(s22)

    del vis[::-2]  # Delete every other element from the end

    expected_len = original_len // 2
    assert len(vis) == expected_len
    # The remaining items should be s22[0], s22[2], s22[4], ...
    for i in range(expected_len):
        assert vis[i] == s22[i * 2]

    # --- Test that deleting an empty slice is a no-op ---
    # reset
    del vis[:]
    vis.extend(s22)

    del vis[5:2]  # This slice is empty
    assert len(vis) == original_len

    del vis[100:200]  # This slice is also empty (out of bounds)
    assert len(vis) == original_len

    # Verify content is unchanged
    for i in range(original_len):
        assert vis[i] == s22[i]


def test_setitem_slice_negative_step(server, s22):
    """Test assigning values to a slice with a negative step."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.cache = None
    vis.extend(s22)
    new_atoms = molecule("H2O")

    # The slice vis[10:0:-2] corresponds to indices [10, 8, 6, 4, 2]
    # It has 5 elements.
    slice_len = 5

    new_values = [new_atoms.copy() for _ in range(slice_len)]
    vis[10:0:-2] = new_values

    assert len(vis) == len(s22)
    # Check that the items at the affected indices have been replaced
    assert vis[10] == new_atoms
    assert vis[8] == new_atoms
    assert vis[6] == new_atoms
    assert vis[4] == new_atoms
    assert vis[2] == new_atoms

    # Check that other items are untouched
    assert vis[0] == s22[0]
    assert vis[1] == s22[1]
    assert vis[3] == s22[3]
    assert vis[11] == s22[11]

    # Test that a length mismatch still raises ValueError
    with pytest.raises(ValueError):
        vis[::-1] = [new_atoms.copy()]  # Mismatched length


# def test_delitem_list_with_duplicates(server, s22):
#     """Test deleting with a list of indices containing duplicates."""
#     vis = ZnDraw(url=server, room="testroom", user="testuser")
#     vis.client._cache = LRUCache(maxsize=1)
#     vis.extend(s22)
#     original_len = len(s22)

#     # A robust implementation should probably forbid duplicates to avoid ambiguity
#     # in a transactional client-server setting.
#     with pytest.raises(ValueError, match="Duplicate indices are not allowed in bulk operations"):
#          del vis[[5, 10, 5]]

#     # Ensure no change was made
#     assert len(vis) == original_len


# def test_setitem_list_with_duplicates(server, s22):
#     """Test setting with a list of indices containing duplicates."""
#     vis = ZnDraw(url=server, room="testroom", user="testuser")
#     vis.client._cache = LRUCache(maxsize=1)
#     vis.extend(s22)
#     original_len = len(s22)
#     new_atoms = molecule("H2O")

#     with pytest.raises(ValueError, match="Duplicate indices are not allowed in bulk operations"):
#          vis[[5, 10, 5]] = [new_atoms.copy() for _ in range(3)]

#     # Ensure no change was made
#     assert len(vis) == original_len

# TODO: slice beyond available range
