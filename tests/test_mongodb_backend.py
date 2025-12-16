"""Tests for MongoDB storage backend."""

import uuid

import numpy as np
import pytest
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

from zndraw.storage import MongoDBStorageBackend

# MongoDB connection URL for tests
MONGODB_TEST_URL = "mongodb://root:example@localhost:27017/"
MONGODB_TEST_DATABASE = "zndraw_test"


def is_mongodb_available() -> bool:
    """Check if MongoDB is available for testing."""
    client = None
    try:
        client = MongoClient(MONGODB_TEST_URL, serverSelectionTimeoutMS=1000)
        client.admin.command("ping")
        return True
    except ServerSelectionTimeoutError:
        return False
    except Exception:
        # Other connection errors (auth, network, etc.)
        return False
    finally:
        if client:
            client.close()


# Skip all tests in this module if MongoDB is unavailable
pytestmark = pytest.mark.skipif(
    not is_mongodb_available(),
    reason="MongoDB not available at mongodb://root:example@localhost:27017/",
)


@pytest.fixture
def mongodb_client():
    """Create MongoDB client for test setup/cleanup."""
    client = MongoClient(MONGODB_TEST_URL)
    yield client
    client.close()


@pytest.fixture
def clean_collection(mongodb_client):
    """Fixture that provides a clean collection name and cleans up after test."""
    collection_name = f"test_room_{uuid.uuid4().hex[:8]}"
    yield collection_name
    # Cleanup
    mongodb_client[MONGODB_TEST_DATABASE].drop_collection(collection_name)


@pytest.fixture
def storage(clean_collection):
    """Create a MongoDBStorageBackend instance for testing."""
    return MongoDBStorageBackend(
        uri=MONGODB_TEST_URL,
        database=MONGODB_TEST_DATABASE,
        room_id=clean_collection,
    )


@pytest.fixture
def sample_frame():
    """Create a sample frame for testing."""
    return {
        b"arrays.positions": np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]).tobytes(),
        b"arrays.numbers": np.array([1, 6]).tobytes(),
        b"info.energy": b"-10.5",
    }


@pytest.fixture
def sample_frames():
    """Create multiple sample frames for testing."""
    frames = []
    for i in range(5):
        frame = {
            b"arrays.positions": np.array(
                [[float(i), 0.0, 0.0], [float(i) + 1, 1.0, 1.0]]
            ).tobytes(),
            b"arrays.numbers": np.array([1, 6]).tobytes(),
            b"info.energy": f"-{i * 10}.5".encode(),
        }
        frames.append(frame)
    return frames


def test_storage_empty_initially(storage):
    """Test that storage is empty when first created."""
    assert len(storage) == 0


def test_extend_single_frame(storage, sample_frame):
    """Test extending storage with a single frame."""
    storage.extend([sample_frame])
    assert len(storage) == 1


def test_extend_multiple_frames(storage, sample_frames):
    """Test extending storage with multiple frames."""
    storage.extend(sample_frames)
    assert len(storage) == 5


def test_append_single_frame(storage, sample_frame):
    """Test appending a single frame."""
    storage.append(sample_frame)
    assert len(storage) == 1


def test_get_single_frame(storage, sample_frame):
    """Test retrieving a single frame by index."""
    storage.extend([sample_frame])
    retrieved = storage.get(0)

    assert b"arrays.positions" in retrieved
    assert b"arrays.numbers" in retrieved
    assert b"info.energy" in retrieved


def test_get_negative_index(storage, sample_frames):
    """Test retrieving a frame with negative index."""
    storage.extend(sample_frames)
    last_frame = storage.get(-1)
    first_frame = storage.get(0)

    # Last frame should have different positions than first
    assert last_frame[b"arrays.positions"] != first_frame[b"arrays.positions"]


def test_get_with_list_of_indices(storage, sample_frames):
    """Test retrieving multiple frames by list of indices."""
    storage.extend(sample_frames)
    retrieved = storage.get([0, 2, 4])

    assert len(retrieved) == 3


def test_get_with_slice(storage, sample_frames):
    """Test retrieving frames with slice."""
    storage.extend(sample_frames)
    retrieved = storage.get(slice(1, 4))

    assert len(retrieved) == 3


def test_get_with_numpy_array(storage, sample_frames):
    """Test retrieving frames with numpy array of indices."""
    storage.extend(sample_frames)
    indices = np.array([0, 2, 4])
    retrieved = storage.get(indices)

    assert len(retrieved) == 3


def test_get_with_numpy_scalar(storage, sample_frame):
    """Test retrieving frame with numpy scalar index."""
    storage.extend([sample_frame])
    index = np.int64(0)
    retrieved = storage.get(index)

    assert isinstance(retrieved, dict)
    assert b"arrays.positions" in retrieved


def test_get_with_key_filter(storage, sample_frame):
    """Test retrieving frame with key filtering."""
    storage.extend([sample_frame])
    retrieved = storage.get(0, keys=["arrays.positions"])

    assert b"arrays.positions" in retrieved
    assert b"arrays.numbers" not in retrieved
    assert b"info.energy" not in retrieved


def test_get_out_of_bounds_raises(storage, sample_frame):
    """Test that out of bounds index raises IndexError."""
    storage.extend([sample_frame])

    with pytest.raises(IndexError):
        storage.get(10)


def test_get_negative_out_of_bounds_raises(storage, sample_frame):
    """Test that negative out of bounds index raises IndexError."""
    storage.extend([sample_frame])

    with pytest.raises(IndexError):
        storage.get(-10)


def test_getitem_single(storage, sample_frame):
    """Test __getitem__ with single index."""
    storage.extend([sample_frame])
    retrieved = storage[0]

    assert isinstance(retrieved, dict)
    assert b"arrays.positions" in retrieved


def test_getitem_slice(storage, sample_frames):
    """Test __getitem__ with slice."""
    storage.extend(sample_frames)
    retrieved = storage[1:4]

    assert len(retrieved) == 3


def test_get_available_keys(storage, sample_frame):
    """Test listing available keys for a frame."""
    storage.extend([sample_frame])
    keys = storage.get_available_keys(0)

    assert "arrays.positions" in keys
    assert "arrays.numbers" in keys
    assert "info.energy" in keys


def test_get_available_keys_out_of_bounds(storage, sample_frame):
    """Test get_available_keys with out of bounds index."""
    storage.extend([sample_frame])

    with pytest.raises(IndexError):
        storage.get_available_keys(10)


def test_sequential_extend(storage, sample_frame):
    """Test multiple extend calls maintain correct indices."""
    storage.extend([sample_frame])
    storage.extend([sample_frame])
    storage.extend([sample_frame])

    assert len(storage) == 3

    # Verify all frames are retrievable
    for i in range(3):
        frame = storage.get(i)
        assert b"arrays.positions" in frame


def test_setitem_raises_not_implemented(storage, sample_frame):
    """Test that __setitem__ raises NotImplementedError."""
    storage.extend([sample_frame])

    with pytest.raises(NotImplementedError):
        storage[0] = sample_frame


def test_delitem_raises_not_implemented(storage, sample_frame):
    """Test that __delitem__ raises NotImplementedError."""
    storage.extend([sample_frame])

    with pytest.raises(NotImplementedError):
        del storage[0]


def test_insert_raises_not_implemented(storage, sample_frame):
    """Test that insert raises NotImplementedError."""
    with pytest.raises(NotImplementedError):
        storage.insert(0, sample_frame)


def test_extend_empty_list(storage):
    """Test extending with empty list does nothing."""
    storage.extend([])
    assert len(storage) == 0


@pytest.mark.parametrize("num_frames", [10, 100, 500])
def test_large_batch_extend(storage, num_frames):
    """Test extending with large batch of frames."""
    frames = []
    for i in range(num_frames):
        frame = {
            b"arrays.positions": np.random.rand(10, 3).tobytes(),
            b"info.index": str(i).encode(),
        }
        frames.append(frame)

    storage.extend(frames)
    assert len(storage) == num_frames

    # Verify random access
    for idx in [0, num_frames // 2, num_frames - 1]:
        frame = storage.get(idx)
        assert frame[b"info.index"] == str(idx).encode()
