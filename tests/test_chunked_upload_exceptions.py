"""Tests for exception handling in chunked upload functionality."""
import ase
import numpy as np
import pytest
from unittest.mock import Mock, patch

from zndraw.zndraw import LocalSettings, ZnDraw


def test_keyboard_interrupt_not_caught():
    """Test that KeyboardInterrupt is not caught by retry logic.

    The chunked upload should catch specific exceptions (IOError, ConnectionError,
    TimeoutError) but NOT KeyboardInterrupt, which should propagate immediately.
    """
    # Create a mock ZnDraw instance
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(max_retries=3, retry_delay=0.1)

    # Create test data
    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from asebytes import encode
    from zndraw.utils import update_colors_and_radii
    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(5)]

    # Mock extend_frames to raise KeyboardInterrupt
    vis.extend_frames = Mock(side_effect=KeyboardInterrupt("User cancelled"))

    # Calculate chunks (this should work)
    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Try to upload - KeyboardInterrupt should propagate immediately
    # without retries
    with pytest.raises(KeyboardInterrupt):
        # Simulate the upload loop from extend()
        for chunk_idx, (chunk, chunk_size) in enumerate(zip(chunks, chunk_sizes)):
            for attempt in range(vis.local.max_retries + 1):
                try:
                    vis.extend_frames(chunk)
                    break
                except KeyboardInterrupt:
                    # Should re-raise immediately, not retry
                    raise
                except Exception as e:
                    if attempt == vis.local.max_retries:
                        raise

    # Verify extend_frames was called only ONCE (no retries)
    assert vis.extend_frames.call_count == 1


def test_connection_error_triggers_retry():
    """Test that ConnectionError triggers retry logic with exponential backoff."""
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(max_retries=3, retry_delay=0.01)  # Fast for testing

    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from asebytes import encode
    from zndraw.utils import update_colors_and_radii
    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(2)]

    # Mock extend_frames to fail twice, then succeed
    call_count = 0
    def mock_extend_frames(chunk):
        nonlocal call_count
        call_count += 1
        if call_count <= 2:
            raise ConnectionError("Connection failed")
        # Success on third attempt
        return None

    vis.extend_frames = Mock(side_effect=mock_extend_frames)

    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Upload should succeed after retries
    import time
    with patch('time.sleep') as mock_sleep:  # Mock sleep to speed up test
        for chunk_idx, (chunk, chunk_size) in enumerate(zip(chunks, chunk_sizes)):
            for attempt in range(vis.local.max_retries + 1):
                try:
                    vis.extend_frames(chunk)
                    break
                except (ConnectionError, IOError, TimeoutError) as e:
                    if attempt == vis.local.max_retries:
                        raise
                    else:
                        delay = vis.local.retry_delay * (2**attempt)
                        time.sleep(delay)

    # Should have been called 3 times (2 failures + 1 success)
    assert vis.extend_frames.call_count == 3

    # Should have slept twice (after first two failures)
    assert mock_sleep.call_count == 2


def test_max_retries_exceeded_raises_error():
    """Test that after max_retries, the exception is re-raised."""
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(max_retries=2, retry_delay=0.01)

    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from asebytes import encode
    from zndraw.utils import update_colors_and_radii
    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(2)]

    # Mock extend_frames to always fail
    vis.extend_frames = Mock(side_effect=ConnectionError("Always fails"))

    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Should raise ConnectionError after exhausting retries
    with patch('time.sleep'):  # Mock sleep to speed up test
        with pytest.raises(ConnectionError, match="Always fails"):
            for chunk_idx, (chunk, chunk_size) in enumerate(zip(chunks, chunk_sizes)):
                for attempt in range(vis.local.max_retries + 1):
                    try:
                        vis.extend_frames(chunk)
                        break
                    except (ConnectionError, IOError, TimeoutError) as e:
                        if attempt == vis.local.max_retries:
                            raise
                        else:
                            delay = vis.local.retry_delay * (2**attempt)
                            import time
                            time.sleep(delay)

    # Should have been called max_retries + 1 times
    assert vis.extend_frames.call_count == vis.local.max_retries + 1


def test_specific_exceptions_only():
    """Test that only specific exceptions (ConnectionError, IOError, TimeoutError) trigger retries."""
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(max_retries=3, retry_delay=0.01)

    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from asebytes import encode
    from zndraw.utils import update_colors_and_radii
    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(2)]

    # Mock extend_frames to raise ValueError (not a retryable exception)
    vis.extend_frames = Mock(side_effect=ValueError("Bad value"))

    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # ValueError should propagate immediately without retries
    with pytest.raises(ValueError, match="Bad value"):
        for chunk_idx, (chunk, chunk_size) in enumerate(zip(chunks, chunk_sizes)):
            for attempt in range(vis.local.max_retries + 1):
                try:
                    vis.extend_frames(chunk)
                    break
                except (ConnectionError, IOError, TimeoutError) as e:
                    if attempt == vis.local.max_retries:
                        raise
                    else:
                        delay = vis.local.retry_delay * (2**attempt)
                        import time
                        time.sleep(delay)

    # Should have been called only ONCE (no retries for ValueError)
    assert vis.extend_frames.call_count == 1


def test_exponential_backoff_timing():
    """Test that exponential backoff increases delay correctly."""
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(max_retries=3, retry_delay=1.0)

    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from asebytes import encode
    from zndraw.utils import update_colors_and_radii
    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(2)]

    # Mock extend_frames to fail 3 times, then succeed
    call_count = 0
    def mock_extend_frames(chunk):
        nonlocal call_count
        call_count += 1
        if call_count <= 3:
            raise ConnectionError("Connection failed")
        return None

    vis.extend_frames = Mock(side_effect=mock_extend_frames)

    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Track sleep delays
    sleep_delays = []
    def track_sleep(delay):
        sleep_delays.append(delay)

    with patch('time.sleep', side_effect=track_sleep):
        for chunk_idx, (chunk, chunk_size) in enumerate(zip(chunks, chunk_sizes)):
            for attempt in range(vis.local.max_retries + 1):
                try:
                    vis.extend_frames(chunk)
                    break
                except (ConnectionError, IOError, TimeoutError) as e:
                    if attempt == vis.local.max_retries:
                        raise
                    else:
                        delay = vis.local.retry_delay * (2**attempt)
                        import time
                        time.sleep(delay)

    # Verify exponential backoff: 1.0, 2.0, 4.0
    assert len(sleep_delays) == 3
    assert sleep_delays[0] == 1.0  # 1.0 * 2^0
    assert sleep_delays[1] == 2.0  # 1.0 * 2^1
    assert sleep_delays[2] == 4.0  # 1.0 * 2^2
