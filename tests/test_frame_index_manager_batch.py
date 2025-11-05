"""Test batch append performance improvements for FrameIndexManager."""
import time

import pytest
from znsocket import MemoryStorage

from zndraw.app.frame_index_manager import FrameIndexManager


def test_append_batch_correctness():
    """Test that append_batch produces the same result as multiple append calls."""
    redis = MemoryStorage()
    key = "test:frames"

    # Test with append
    manager1 = FrameIndexManager(redis, key + ":1")
    for i in range(10):
        manager1.append(f"frame:{i}")
    result1 = manager1.get_all(withscores=True)

    # Test with append_batch
    manager2 = FrameIndexManager(redis, key + ":2")
    members = [f"frame:{i}" for i in range(10)]
    manager2.append_batch(members)
    result2 = manager2.get_all(withscores=True)

    # Compare results (members should be the same, scores should be sequential)
    assert len(result1) == len(result2)
    for (member1, score1), (member2, score2) in zip(result1, result2):
        assert member1 == member2
        # Scores should be sequential integers starting from 1.0
        assert score1 == score2


def test_append_batch_empty():
    """Test that append_batch handles empty list correctly."""
    redis = MemoryStorage()
    manager = FrameIndexManager(redis, "test:empty")

    scores = manager.append_batch([])
    assert scores == []
    assert len(manager) == 0


def test_append_batch_after_existing():
    """Test that append_batch works correctly after existing frames."""
    redis = MemoryStorage()
    manager = FrameIndexManager(redis, "test:after")

    # Add some frames normally
    manager.append("frame:0")
    manager.append("frame:1")

    # Now batch append more
    manager.append_batch(["frame:2", "frame:3", "frame:4"])

    # Verify all frames are in order
    all_frames = manager.get_all()
    assert len(all_frames) == 5
    for i, frame in enumerate(all_frames):
        expected = f"frame:{i}"
        assert frame == expected or frame == expected.encode()


def test_append_batch_performance():
    """Test that append_batch is significantly faster than multiple appends."""
    redis = MemoryStorage()

    # Measure time for individual appends
    manager1 = FrameIndexManager(redis, "test:perf:individual")
    start = time.time()
    for i in range(1000):
        manager1.append(f"frame:{i}")
    individual_time = time.time() - start

    # Measure time for batch append
    manager2 = FrameIndexManager(redis, "test:perf:batch")
    members = [f"frame:{i}" for i in range(1000)]
    start = time.time()
    manager2.append_batch(members)
    batch_time = time.time() - start

    # Batch should be much faster (at least 2x, typically much more)
    print(f"\nIndividual appends (1000): {individual_time:.4f}s")
    print(f"Batch append (1000): {batch_time:.4f}s")
    print(f"Speedup: {individual_time / batch_time:.1f}x")

    assert batch_time < individual_time / 2, (
        f"Batch append should be at least 2x faster, "
        f"but was only {individual_time / batch_time:.1f}x faster"
    )


def test_append_batch_scores_sequential():
    """Test that batch append assigns sequential scores."""
    redis = MemoryStorage()
    manager = FrameIndexManager(redis, "test:sequential")

    members = [f"frame:{i}" for i in range(100)]
    scores = manager.append_batch(members)

    # Scores should be sequential floats starting from 1.0
    for i, score in enumerate(scores):
        assert score == 1.0 + i

    # Verify in Redis as well
    all_frames = manager.get_all(withscores=True)
    for i, (member, score) in enumerate(all_frames):
        assert score == 1.0 + i
