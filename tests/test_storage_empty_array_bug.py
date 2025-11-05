"""Tests for the division by zero bug when storing empty arrays in Zarr.

This test file addresses the error:
ZeroDivisionError: division by zero
That occurs in zarr.core.common.ceildiv when chunk size is 0.

The fix ensures that chunks never have 0 dimensions, even when array shapes do.
"""

import numpy as np
import numpy.testing as npt
import pytest
import zarr
from zarr.storage import MemoryStore

from zndraw.storage import ZarrStorageSequence, extend_zarr


def test_empty_array_initial_storage():
    """Test storing an empty array as the first entry works correctly."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # Should work without division by zero error
    store.append({"empty": np.array([])})

    # Verify it was stored
    assert len(store) == 1
    result = store[0]
    assert "empty" in result
    npt.assert_array_equal(result["empty"], np.array([]))
    assert result["empty"].shape == (0,)


def test_empty_2d_array_initial_storage():
    """Test storing an empty 2D array works correctly."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # Empty 2D array with shape (0, 3)
    empty_2d = np.array([]).reshape(0, 3)
    store.append({"empty_2d": empty_2d})

    assert len(store) == 1
    result = store[0]
    npt.assert_array_equal(result["empty_2d"], empty_2d)
    assert result["empty_2d"].shape == (0, 3)


def test_empty_array_after_non_empty():
    """Test storing an empty array after a non-empty one works correctly."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # First store a non-empty array
    store.append({"data": np.array([1, 2, 3])})

    # Then store an empty array
    store.append({"data": np.array([])})

    assert len(store) == 2
    npt.assert_array_equal(store[0]["data"], np.array([1, 2, 3]))
    npt.assert_array_equal(store[1]["data"], np.array([]))


def test_mixed_empty_and_non_empty_arrays():
    """Test storing a frame with both empty and non-empty arrays."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    store.append({
        "positions": np.array([[1.0, 2.0, 3.0]]),
        "empty_field": np.array([]),
    })

    assert len(store) == 1
    result = store[0]
    npt.assert_array_equal(result["positions"], np.array([[1.0, 2.0, 3.0]]))
    npt.assert_array_equal(result["empty_field"], np.array([]))


def test_empty_array_in_extend():
    """Test that extend_zarr handles empty arrays correctly."""
    root = zarr.group(store=MemoryStore())

    extend_zarr(root, [{"empty": np.array([])}])

    # Verify the array was created with valid chunks
    arr = root["empty"]
    assert arr.shape == (1, 0)
    # Chunks should have been corrected to avoid 0 dimensions
    assert all(c > 0 for c in arr.chunks), f"Chunks should not contain 0: {arr.chunks}"


@pytest.mark.parametrize("shape", [
    (0,),           # 1D empty
    (0, 3),         # 2D empty with fixed second dimension
    (0, 0),         # 2D completely empty
    (0, 2, 3),      # 3D empty with fixed dimensions
])
def test_various_empty_shapes(shape):
    """Test various empty array shapes are all handled correctly."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    empty_array = np.zeros(shape)
    store.append({"empty": empty_array})

    assert len(store) == 1
    result = store[0]
    npt.assert_array_equal(result["empty"], empty_array)

    # Verify chunks don't have 0 dimensions
    arr = root["empty"]
    assert all(c > 0 for c in arr.chunks), f"Chunks should not contain 0: {arr.chunks}"


def test_chunk_validation_on_creation():
    """Test that chunks are validated to not contain 0 when creating arrays."""
    root = zarr.group(store=MemoryStore())

    # Store data with 0 dimension
    extend_zarr(root, [{"data": np.array([]).reshape(0, 3)}])

    # Check that the array was created with valid chunks
    arr = root["data"]
    assert arr.shape == (1, 0, 3)
    # All chunk dimensions should be > 0
    assert all(c > 0 for c in arr.chunks), f"Chunks must be > 0, got: {arr.chunks}"


def test_chunk_validation_on_resize():
    """Test that resizing arrays with 0 dimensions doesn't cause errors."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # Create initial entry with empty array
    store.append({"data": np.array([]).reshape(0, 3)})

    # Append another empty array (may trigger resize)
    store.append({"data": np.array([]).reshape(0, 3)})

    assert len(store) == 2
    npt.assert_array_equal(store[0]["data"], np.array([]).reshape(0, 3))
    npt.assert_array_equal(store[1]["data"], np.array([]).reshape(0, 3))
