"""Unit tests for FrameIndexManager."""

import pytest
import redis

from zndraw.app.frame_index_manager import FrameIndexManager


@pytest.fixture
def redis_client():
    """Create a Redis client and clean up after test."""
    client = redis.Redis(host="localhost", port=6379, decode_responses=True)
    test_key = "test:frame:indices"
    yield client, test_key
    # Cleanup
    client.delete(test_key)


def test_append_first_frame(redis_client):
    """Test appending the first frame."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    score = manager.append("room:0")
    assert score == 1.0
    assert manager.get_count() == 1


def test_append_multiple_frames(redis_client):
    """Test appending multiple frames."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    scores = []
    for i in range(5):
        score = manager.append(f"room:{i}")
        scores.append(score)

    assert manager.get_count() == 5
    assert scores == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_insert_at_beginning(redis_client):
    """Test inserting at the beginning."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Setup: add frames 0, 1, 2
    for i in range(3):
        manager.append(f"room:{i}")

    # Insert at position 0
    score = manager.insert(0, "room:new")
    assert score == 0.0  # 1.0 - 1.0

    # Verify order
    frames = manager.get_all()
    assert len(frames) == 4
    assert frames[0] == "room:new"
    assert frames[1] == "room:0"


def test_insert_at_middle(redis_client):
    """Test inserting in the middle."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Setup: add frames 0, 1, 2
    for i in range(3):
        manager.append(f"room:{i}")

    # Insert at position 1 (between 0 and 1)
    score = manager.insert(1, "room:new")
    assert score == 1.5  # (1.0 + 2.0) / 2

    # Verify order
    frames = manager.get_all()
    assert len(frames) == 4
    assert frames[0] == "room:0"
    assert frames[1] == "room:new"
    assert frames[2] == "room:1"


def test_insert_at_end(redis_client):
    """Test inserting at the end."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Setup
    for i in range(3):
        manager.append(f"room:{i}")

    # Insert at the end (position 3)
    score = manager.insert(3, "room:new")
    assert score == 4.0  # 3.0 + 1.0

    # Verify order
    frames = manager.get_all()
    assert len(frames) == 4
    assert frames[3] == "room:new"


def test_delete_member(redis_client):
    """Test deleting a member."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Setup
    for i in range(3):
        manager.append(f"room:{i}")

    # Delete middle frame
    result = manager.delete("room:1")
    assert result == 1

    # Verify
    frames = manager.get_all()
    assert len(frames) == 2
    assert frames[0] == "room:0"
    assert frames[1] == "room:2"


def test_delete_at_position(redis_client):
    """Test deleting by position."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Setup
    for i in range(3):
        manager.append(f"room:{i}")

    # Delete position 1
    result = manager.delete_at_position(1)
    assert result == 1

    # Verify
    frames = manager.get_all()
    assert len(frames) == 2
    assert frames[0] == "room:0"
    assert frames[1] == "room:2"


def test_renormalize(redis_client):
    """Test renormalization to contiguous integers."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Create non-contiguous scores by inserting between frames
    manager.append("room:0")  # score: 1.0
    manager.append("room:1")  # score: 2.0
    manager.insert(1, "room:new")  # score: 1.5

    # Renormalize
    count = manager.renormalize()
    assert count == 3

    # Verify scores are now 0, 1, 2
    frames = manager.get_all(withscores=True)
    assert frames[0][1] == 0.0
    assert frames[1][1] == 1.0
    assert frames[2][1] == 2.0

    # Verify order is preserved
    assert frames[0][0] == "room:0"
    assert frames[1][0] == "room:new"
    assert frames[2][0] == "room:1"


def test_get_range(redis_client):
    """Test getting a range of frames."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Setup
    for i in range(5):
        manager.append(f"room:{i}")

    # Get range
    frames = manager.get_range(1, 3)
    assert len(frames) == 3
    assert frames[0] == "room:1"
    assert frames[1] == "room:2"
    assert frames[2] == "room:3"


def test_insert_into_empty(redis_client):
    """Test inserting into an empty set."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    score = manager.insert(0, "room:0")
    assert score == 1.0
    assert manager.get_count() == 1


def test_negative_position_insert(redis_client):
    """Test inserting with negative position."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    # Setup
    for i in range(3):
        manager.append(f"room:{i}")

    # Insert at position -1 (before last)
    score = manager.insert(-1, "room:new")

    # Verify order
    frames = manager.get_all()
    assert len(frames) == 4
    # Negative position -1 means insert before the last element
    # With 3 elements (0, 1, 2), position -1 is position 3 (after conversion)
    # So it inserts at the end
    assert frames[3] == "room:new"


def test_multiple_inserts_same_position(redis_client):
    """Test multiple insertions at the same position."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    manager.append("room:0")  # score: 1.0
    manager.append("room:1")  # score: 2.0

    # Insert at position 1 multiple times
    manager.insert(1, "room:a")  # score: 1.5
    manager.insert(1, "room:b")  # score: 1.25

    frames = manager.get_all()
    assert len(frames) == 4
    assert frames[0] == "room:0"
    assert frames[1] == "room:b"
    assert frames[2] == "room:a"
    assert frames[3] == "room:1"


def test_delete_nonexistent(redis_client):
    """Test deleting a non-existent member."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    result = manager.delete("room:999")
    assert result == 0


def test_delete_position_out_of_range(redis_client):
    """Test deleting at an out-of-range position."""
    r, key = redis_client
    manager = FrameIndexManager(r, key)

    manager.append("room:0")

    result = manager.delete_at_position(10)
    assert result == 0
