"""Tests for room name generation utility."""

import pytest

from zndraw.utils import generate_room_name


def test_generate_room_name_short_filename():
    """Test that short filenames are used as-is."""
    result = generate_room_name("data.xyz", None)
    assert result == "data.xyz"


def test_generate_room_name_truncates_long_filename():
    """Test that long filenames are truncated to max_length."""
    long_name = "my_very_long_file_name_that_exceeds_twenty_chars.xyz"
    result = generate_room_name(long_name, None, max_length=20)
    assert len(result) <= 20
    assert result == "my_very_long_file_na"


def test_generate_room_name_collision_adds_hash(redis_client):
    """Test that collision detection adds hash suffix."""
    # Create existing room
    redis_client.zadd("room:test_file.xyz:trajectory:indices", {"frame_0": 0})

    # Generate name for file that would create same room name
    result = generate_room_name("test_file.xyz", redis_client, max_length=20)
    assert result.startswith("test_file.xyz_")
    # Hash should be appended
    assert len(result) > len("test_file.xyz")
    assert "_" in result


@pytest.mark.parametrize(
    "filename,max_len,expected",
    [
        ("a.xyz", 20, "a.xyz"),
        ("my_file.pdb", 20, "my_file.pdb"),
        ("data/structure.h5", 20, "data/structure.h5"),
        ("very_long_filename_that_exceeds_limit.xyz", 20, "very_long_filename_t"),
    ],
)
def test_generate_room_name_various_inputs(filename, max_len, expected):
    """Test room name generation with various inputs."""
    result = generate_room_name(filename, None, max_length=max_len)
    assert result == expected


def test_generate_room_name_no_extension():
    """Test filename without extension."""
    result = generate_room_name("myfile", None)
    assert result == "myfile"


def test_generate_room_name_with_extension():
    """Test filename with extension is kept."""
    result = generate_room_name(".xyz", None)
    assert isinstance(result, str)
    assert result == ".xyz"


def test_generate_room_name_custom_max_length():
    """Test with custom max length."""
    result = generate_room_name("very_long_filename.xyz", None, max_length=10)
    assert len(result) <= 10
    assert result == "very_long_"


def test_generate_room_name_unicode():
    """Test with unicode characters."""
    result = generate_room_name("café_données.xyz", None)
    # Should handle gracefully
    assert isinstance(result, str)


def test_generate_room_name_with_path():
    """Test that path is included if given."""
    result = generate_room_name("path/to/file.xyz", None)
    assert result == "path/to/file.xyz"


def test_generate_room_name_no_redis_client():
    """Test without Redis client (no collision check)."""
    result = generate_room_name("test.xyz", redis_client=None)
    assert result == "test.xyz"


def test_generate_room_name_multiple_collisions(redis_client):
    """Test handling multiple files with same truncated name."""
    # Create first room
    truncated = "my_very_long_file_na"  # 20 chars
    redis_client.zadd(f"room:{truncated}:trajectory:indices", {"frame_0": 0})

    # Generate for file that truncates to same name
    result1 = generate_room_name(
        "my_very_long_file_name_A.xyz", redis_client, max_length=20
    )
    result2 = generate_room_name(
        "my_very_long_file_name_B.xyz", redis_client, max_length=20
    )

    # Both should get hash suffixes but different hashes
    assert "_" in result1
    assert "_" in result2
    assert result1 != result2
    assert len(result1) <= 20
    assert len(result2) <= 20


def test_generate_room_name_preserves_content():
    """Test that filename content is preserved."""
    result = generate_room_name("my-file-name.xyz", None)
    assert result == "my-file-name.xyz"


def test_generate_room_name_special_chars():
    """Test filename with special characters."""
    result = generate_room_name("@#$%.xyz", None)
    assert isinstance(result, str)
    assert result == "@#$%.xyz"


def test_generate_room_name_dots():
    """Test filename with multiple dots."""
    result = generate_room_name("file.v2.0.xyz", None)
    assert result == "file.v2.0.xyz"


def test_generate_room_name_numbers():
    """Test filename with numbers."""
    result = generate_room_name("file_123_test.xyz", None)
    assert result == "file_123_test.xyz"


def test_generate_room_name_uppercase():
    """Test that uppercase is preserved."""
    result = generate_room_name("MyFile.xyz", None)
    assert result == "MyFile.xyz"


def test_generate_room_name_existing_underscore():
    """Test filename that already has underscores."""
    result = generate_room_name("my_existing_name.xyz", None)
    assert result == "my_existing_name.xyz"


def test_generate_room_name_hash_suffix_length():
    """Test that hash suffix is correct length."""
    redis_client = pytest.importorskip("redis").Redis()
    redis_client.zadd("room:test.xyz:trajectory:indices", {"frame_0": 0})

    result = generate_room_name("test.xyz", redis_client)
    # Should have _XXXX suffix where XXXX is 4 hex chars
    assert result.startswith("test.xyz_")
    suffix = result.split("_")[-1]
    assert len(suffix) == 4
    assert all(c in "0123456789abcdef" for c in suffix)
