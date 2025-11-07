"""Tests for chunked upload functionality."""
import ase
import numpy as np
import pytest

from zndraw.zndraw import LocalSettings


def test_local_settings_defaults():
    """Test that LocalSettings has correct default values."""
    settings = LocalSettings()

    assert settings.target_chunk_size_bytes == 500_000
    assert settings.show_progress is True
    assert settings.max_retries == 3
    assert settings.retry_delay == 1.0


def test_local_settings_validation():
    """Test that LocalSettings validates input correctly."""
    # Valid settings
    settings = LocalSettings(
        target_chunk_size_bytes=1_000_000,
        show_progress=False,
        max_retries=5,
        retry_delay=2.0,
    )
    assert settings.target_chunk_size_bytes == 1_000_000

    # Invalid: target_chunk_size_bytes must be > 0
    with pytest.raises(ValueError):
        LocalSettings(target_chunk_size_bytes=0)

    with pytest.raises(ValueError):
        LocalSettings(target_chunk_size_bytes=-1)

    # Invalid: max_retries must be >= 0
    with pytest.raises(ValueError):
        LocalSettings(max_retries=-1)

    # Invalid: retry_delay must be >= 0
    with pytest.raises(ValueError):
        LocalSettings(retry_delay=-1.0)


def test_calculate_chunk_boundaries_single_chunk():
    """Test chunk calculation when all frames fit in one chunk."""
    from zndraw.zndraw import ZnDraw
    from unittest.mock import Mock

    # Create a mock ZnDraw instance
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(target_chunk_size_bytes=10_000_000)  # 10MB

    # Create small test data (will fit in one chunk)
    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from asebytes import encode
    from zndraw.utils import update_colors_and_radii

    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(10)]

    # Call the method
    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Should be single chunk
    assert len(chunks) == 1
    assert len(chunks[0]) == 10
    assert chunk_sizes[0] > 0


def test_calculate_chunk_boundaries_multiple_chunks():
    """Test chunk calculation when frames need to be split."""
    from zndraw.zndraw import ZnDraw
    from unittest.mock import Mock

    # Create a mock ZnDraw instance
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(target_chunk_size_bytes=2_000)  # 2KB - very small

    # Create test data
    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from asebytes import encode
    from zndraw.utils import update_colors_and_radii

    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(10)]

    # Call the method
    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Should have multiple chunks
    assert len(chunks) > 1

    # All chunks should have sizes
    assert len(chunks) == len(chunk_sizes)

    # Total frames should match
    total_frames = sum(len(chunk) for chunk in chunks)
    assert total_frames == 10

    # Each chunk size should be reasonable
    for size in chunk_sizes:
        assert size > 0
        assert size <= vis.local.target_chunk_size_bytes * 2  # Allow some overflow


def test_calculate_chunk_boundaries_exact_sizes():
    """Test that chunk sizes are calculated from exact encoded sizes."""
    from zndraw.zndraw import ZnDraw
    from asebytes import encode
    import msgpack
    from unittest.mock import Mock

    # Create a mock ZnDraw instance
    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(target_chunk_size_bytes=500_000)

    # Create test data
    atoms = ase.Atoms("H2O", positions=np.random.rand(3, 3))
    from zndraw.utils import update_colors_and_radii

    update_colors_and_radii(atoms)
    dicts = [encode(atoms) for _ in range(100)]

    # Calculate expected size manually (dicts are already msgpack format)
    packed = [msgpack.packb(d) for d in dicts]
    expected_total_size = sum(len(p) for p in packed)

    # Call the method
    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Total size should match
    actual_total_size = sum(chunk_sizes)
    assert actual_total_size == expected_total_size


def test_calculate_chunk_boundaries_empty_list():
    """Test chunk calculation with empty list."""
    from zndraw.zndraw import ZnDraw
    from unittest.mock import Mock

    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings()

    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, [])

    assert len(chunks) == 0
    assert len(chunk_sizes) == 0


def test_local_settings_modification():
    """Test that LocalSettings can be modified after creation."""
    settings = LocalSettings()

    # Modify settings
    settings.target_chunk_size_bytes = 1_000_000
    settings.show_progress = False
    settings.max_retries = 5
    settings.retry_delay = 2.5

    assert settings.target_chunk_size_bytes == 1_000_000
    assert settings.show_progress is False
    assert settings.max_retries == 5
    assert settings.retry_delay == 2.5


def test_single_frame_exceeds_chunk_size():
    """Test that a single frame larger than target_chunk_size_bytes is handled correctly."""
    from zndraw.zndraw import ZnDraw
    from unittest.mock import Mock

    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(target_chunk_size_bytes=100)  # Very small, 100 bytes

    # Create a large atoms object that will exceed the chunk size
    # 100 atoms with positions will be > 100 bytes when encoded
    positions = np.random.rand(100, 3) * 10
    atoms = ase.Atoms(numbers=[6] * 100, positions=positions)

    from asebytes import encode
    from zndraw.utils import update_colors_and_radii
    update_colors_and_radii(atoms)

    # Create just one frame
    dicts = [encode(atoms)]

    # Calculate chunks - single frame should be in a chunk by itself
    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Should have exactly 1 chunk even though it exceeds target size
    assert len(chunks) == 1
    assert len(chunks[0]) == 1

    # The chunk size will exceed the target (that's expected and OK)
    assert chunk_sizes[0] > vis.local.target_chunk_size_bytes


def test_chunk_size_realistic():
    """Test chunk calculation with realistic data (28 atoms like in the user's case)."""
    from zndraw.zndraw import ZnDraw
    from unittest.mock import Mock

    vis = Mock(spec=ZnDraw)
    vis.local = LocalSettings(target_chunk_size_bytes=500_000)  # 500KB

    # Create atoms similar to user's case (28 atoms)
    positions = np.random.rand(28, 3) * 10
    numbers = [6] * 20 + [8] * 4 + [7] * 4
    atoms = ase.Atoms(numbers=numbers, positions=positions)

    from asebytes import encode
    from zndraw.utils import update_colors_and_radii

    update_colors_and_radii(atoms)

    # Create 1000 frames like in the user's scenario
    dicts = [encode(atoms) for _ in range(1000)]

    # Calculate chunks
    chunks, chunk_sizes = ZnDraw._calculate_chunk_boundaries(vis, dicts)

    # Should have multiple chunks (1000 frames at ~1.6KB each = ~1.6MB total)
    # With 500KB target, should get ~3-4 chunks
    assert 3 <= len(chunks) <= 5

    # Total frames should be 1000
    total_frames = sum(len(chunk) for chunk in chunks)
    assert total_frames == 1000

    # Each chunk should be roughly 500KB or less (allow some overhead)
    for size in chunk_sizes:
        # Allow first frame to exceed target (no prior chunk to put it in)
        assert size <= vis.local.target_chunk_size_bytes * 2

    print(f"\nRealistic test results:")
    print(f"  Total frames: 1000")
    print(f"  Number of chunks: {len(chunks)}")
    print(f"  Chunk sizes: {[f'{s/1024:.1f}KB' for s in chunk_sizes]}")
    print(f"  Frames per chunk: {[len(c) for c in chunks]}")
