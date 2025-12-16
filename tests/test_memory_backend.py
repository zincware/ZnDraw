"""Tests for InMemoryStorageBackend."""

import numpy as np
import pytest

from zndraw.storage import InMemoryStorageBackend


@pytest.fixture
def backend():
    """Create a fresh in-memory storage backend."""
    return InMemoryStorageBackend()


@pytest.fixture
def populated_backend(backend):
    """Create a backend with some test data."""
    frames = [
        {b"arrays.positions": b"pos1", b"info.energy": b"e1"},
        {b"arrays.positions": b"pos2", b"info.energy": b"e2"},
        {b"arrays.positions": b"pos3", b"info.energy": b"e3"},
    ]
    backend.extend(frames)
    return backend


def test_empty_initially(backend):
    """Test that backend starts empty."""
    assert len(backend) == 0


def test_extend_single_frame(backend):
    """Test extending with a single frame."""
    backend.extend([{b"key": b"value"}])
    assert len(backend) == 1


def test_extend_multiple_frames(backend):
    """Test extending with multiple frames."""
    frames = [{b"key": b"value1"}, {b"key": b"value2"}, {b"key": b"value3"}]
    backend.extend(frames)
    assert len(backend) == 3


def test_extend_empty_list(backend):
    """Test extending with empty list does nothing."""
    backend.extend([])
    assert len(backend) == 0


def test_get_single_index(populated_backend):
    """Test getting a single frame by index."""
    result = populated_backend.get(0)
    assert result == {b"arrays.positions": b"pos1", b"info.energy": b"e1"}


def test_get_negative_index(populated_backend):
    """Test getting a frame with negative index."""
    result = populated_backend.get(-1)
    assert result == {b"arrays.positions": b"pos3", b"info.energy": b"e3"}


def test_get_multiple_indices(populated_backend):
    """Test getting multiple frames by list of indices."""
    result = populated_backend.get([0, 2])
    assert len(result) == 2
    assert result[0] == {b"arrays.positions": b"pos1", b"info.energy": b"e1"}
    assert result[1] == {b"arrays.positions": b"pos3", b"info.energy": b"e3"}


def test_get_slice(populated_backend):
    """Test getting frames with slice."""
    result = populated_backend.get(slice(0, 2))
    assert len(result) == 2
    assert result[0] == {b"arrays.positions": b"pos1", b"info.energy": b"e1"}
    assert result[1] == {b"arrays.positions": b"pos2", b"info.energy": b"e2"}


def test_get_numpy_array_index(populated_backend):
    """Test getting frames with numpy array indices."""
    indices = np.array([0, 2])
    result = populated_backend.get(indices)
    assert len(result) == 2


def test_get_numpy_scalar_index(populated_backend):
    """Test getting frame with numpy scalar index."""
    index = np.int64(1)
    result = populated_backend.get(index)
    assert result == {b"arrays.positions": b"pos2", b"info.energy": b"e2"}


def test_get_with_key_filter(populated_backend):
    """Test getting frames with key filtering."""
    result = populated_backend.get(0, keys=["arrays.positions"])
    assert result == {b"arrays.positions": b"pos1"}
    assert b"info.energy" not in result


def test_get_index_out_of_bounds(populated_backend):
    """Test that out of bounds index raises IndexError."""
    with pytest.raises(IndexError):
        populated_backend.get(10)


def test_get_negative_index_out_of_bounds(populated_backend):
    """Test that negative out of bounds index raises IndexError."""
    with pytest.raises(IndexError):
        populated_backend.get(-10)


def test_get_available_keys(populated_backend):
    """Test getting available keys for a frame."""
    keys = populated_backend.get_available_keys(0)
    assert set(keys) == {"arrays.positions", "info.energy"}


def test_get_available_keys_negative_index(populated_backend):
    """Test getting available keys with negative index."""
    keys = populated_backend.get_available_keys(-1)
    assert set(keys) == {"arrays.positions", "info.energy"}


def test_get_available_keys_out_of_bounds(populated_backend):
    """Test that out of bounds index raises IndexError."""
    with pytest.raises(IndexError):
        populated_backend.get_available_keys(10)


def test_append(backend):
    """Test appending a single frame."""
    backend.append({b"key": b"value"})
    assert len(backend) == 1
    assert backend.get(0) == {b"key": b"value"}


def test_getitem(populated_backend):
    """Test __getitem__ works like get."""
    result = populated_backend[0]
    assert result == {b"arrays.positions": b"pos1", b"info.energy": b"e1"}


def test_setitem_raises(populated_backend):
    """Test that __setitem__ raises NotImplementedError."""
    with pytest.raises(NotImplementedError):
        populated_backend[0] = {b"key": b"value"}


def test_delitem_raises(populated_backend):
    """Test that __delitem__ raises NotImplementedError."""
    with pytest.raises(NotImplementedError):
        del populated_backend[0]


def test_insert_raises(populated_backend):
    """Test that insert raises NotImplementedError."""
    with pytest.raises(NotImplementedError):
        populated_backend.insert(0, {b"key": b"value"})
