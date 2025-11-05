"""Test for race condition bug in batch write optimization.

This test exposes the bug where frames are skipped when all their keys
are batch-written, because the empty dict evaluates to False.
"""
import numpy as np
import pytest
import zarr
from zarr.storage import MemoryStore

from zndraw.storage import (
    ZarrStorageSequence,
    extend_zarr,
    read_zarr,
)


def test_batch_write_all_keys_batchable():
    """Test that frames are not skipped when all keys are batchable.

    When all keys in a frame are batchable (same shape, same dtype across
    all frames), they get batch-written and removed from the entry dict.
    This leaves an empty dict {}, which is falsy in Python.

    The bug: `if entry:` check causes these frames to be skipped entirely.

    This test creates 20 frames (> 10 to trigger optimization) where all
    keys are batchable, and verifies that all frames can be read back.
    """
    store = MemoryStore()
    root = zarr.open_group(store=store, mode="w")

    # Create data where ALL keys are batchable:
    # - Same dtype (float64)
    # - Same shape (3,) for positions, () for energy
    # - Appears in ALL frames
    # - Non-object dtype
    data = []
    for i in range(20):  # > 10 to trigger batch optimization
        frame = {
            "positions": np.array([1.0 * i, 2.0 * i, 3.0 * i], dtype=np.float64),
            "energy": np.array(i * 10.0, dtype=np.float64),
            "numbers": np.array([6, 8, 7], dtype=np.int64),
        }
        data.append(frame)

    # Extend zarr with this data
    extend_zarr(root, data)

    # Verify all frames were written
    assert "__valid_keys__" in root
    assert root["__valid_keys__"].shape[0] == 20

    # Verify we can read back all frames
    for i in range(20):
        frame = read_zarr(root, i)

        # Check positions
        assert "positions" in frame
        expected_positions = np.array([1.0 * i, 2.0 * i, 3.0 * i], dtype=np.float64)
        np.testing.assert_array_equal(frame["positions"], expected_positions)

        # Check energy
        assert "energy" in frame
        expected_energy = i * 10.0
        assert frame["energy"] == expected_energy

        # Check numbers
        assert "numbers" in frame
        expected_numbers = np.array([6, 8, 7], dtype=np.int64)
        np.testing.assert_array_equal(frame["numbers"], expected_numbers)


def test_batch_write_partial_batchable_keys():
    """Test frames where only some keys are batchable.

    This ensures that after batch-writing some keys, the remaining
    keys are still processed correctly via the recursive path.
    """
    store = MemoryStore()
    root = zarr.open_group(store=store, mode="w")

    # Create data where:
    # - "positions" is batchable (same shape across all frames)
    # - "metadata" is NOT batchable (different structure per frame)
    data = []
    for i in range(15):  # > 10 to trigger optimization
        frame = {
            "positions": np.array([1.0 * i, 2.0 * i, 3.0 * i], dtype=np.float64),
            "metadata": {
                "id": i,
                "name": f"frame_{i}",
                # Different nested structure to prevent batching
                "custom": {"value": i * 2} if i % 2 == 0 else None,
            }
        }
        data.append(frame)

    # Extend zarr
    extend_zarr(root, data)

    # Verify all frames
    assert root["__valid_keys__"].shape[0] == 15

    for i in range(15):
        frame = read_zarr(root, i)

        # Check positions (was batch-written)
        assert "positions" in frame
        expected_positions = np.array([1.0 * i, 2.0 * i, 3.0 * i], dtype=np.float64)
        np.testing.assert_array_equal(frame["positions"], expected_positions)

        # Check metadata (was recursively processed)
        assert "metadata" in frame
        assert frame["metadata"]["id"] == i
        assert frame["metadata"]["name"] == f"frame_{i}"
        if i % 2 == 0:
            assert frame["metadata"]["custom"]["value"] == i * 2
        else:
            assert frame["metadata"]["custom"] is None


def test_batch_write_with_zarr_storage_sequence():
    """Test batch write through ZarrStorageSequence interface.

    This is the actual user-facing API, so we need to ensure it works
    correctly with the batch optimization.
    """
    store = MemoryStore()
    root = zarr.open_group(store=store, mode="w")

    # Create initial data
    initial_data = [
        {
            "positions": np.array([0.0, 0.0, 0.0], dtype=np.float64),
            "energy": np.array(0.0, dtype=np.float64),
        }
    ]
    extend_zarr(root, initial_data)

    storage = ZarrStorageSequence(root)

    # Extend with many frames (> 10 to trigger optimization)
    new_frames = []
    for i in range(1, 21):
        frame = {
            "positions": np.array([1.0 * i, 2.0 * i, 3.0 * i], dtype=np.float64),
            "energy": np.array(i * 10.0, dtype=np.float64),
        }
        new_frames.append(frame)

    storage.extend(new_frames)

    # Verify all frames (initial + extended)
    assert len(storage) == 21

    # Check first frame (initial)
    frame0 = storage[0]
    np.testing.assert_array_equal(frame0["positions"], [0.0, 0.0, 0.0])
    assert frame0["energy"] == 0.0

    # Check extended frames
    for i in range(1, 21):
        frame = storage[i]
        expected_positions = np.array([1.0 * i, 2.0 * i, 3.0 * i], dtype=np.float64)
        np.testing.assert_array_equal(frame["positions"], expected_positions)
        assert frame["energy"] == i * 10.0
