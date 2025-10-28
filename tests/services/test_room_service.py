"""Tests for RoomService.

Tests follow CLAUDE.md principles:
- Each test is a function (not a method)
- Tests are specific and test only one thing
- Uses pytest.mark.parametrize to avoid duplication
"""

import json

import pytest

from zndraw.services.room_service import RoomService


def test_room_exists_returns_false_for_new_room(redis_client):
    """Test room_exists returns False for nonexistent room."""
    service = RoomService(redis_client)
    assert service.room_exists("nonexistent") is False


def test_room_exists_returns_true_for_existing_room(redis_client):
    """Test room_exists returns True for existing room."""
    service = RoomService(redis_client)
    redis_client.set("room:test:current_frame", 0)
    assert service.room_exists("test") is True


def test_create_empty_room(redis_client):
    """Test creating an empty room with defaults."""
    service = RoomService(redis_client)
    result = service.create_room("test", "user1")

    assert result["created"] is True
    assert result["frameCount"] == 0
    assert redis_client.exists("room:test:current_frame") > 0
    # Verify default geometries created
    assert redis_client.exists("room:test:geometries") > 0


def test_create_empty_room_with_description(redis_client):
    """Test creating room with description."""
    service = RoomService(redis_client)
    result = service.create_room("test", "user1", description="Test Room")

    assert result["created"] is True
    assert redis_client.get("room:test:description") == "Test Room"


def test_create_empty_room_initializes_metadata(redis_client):
    """Test room creation initializes all required metadata."""
    service = RoomService(redis_client)
    service.create_room("test", "user1")

    assert redis_client.get("room:test:current_frame") == "0"
    assert redis_client.get("room:test:locked") == "0"
    assert redis_client.get("room:test:hidden") == "0"


def test_create_empty_room_creates_default_geometries(redis_client):
    """Test room creation creates all default geometries."""
    service = RoomService(redis_client)
    service.create_room("test", "user1")

    geometries = redis_client.hgetall("room:test:geometries")
    assert "particles" in geometries
    assert "bonds" in geometries
    assert "curve" in geometries
    assert "cell" in geometries
    assert "floor" in geometries

    # Verify geometry structure
    particles = json.loads(geometries["particles"])
    assert "type" in particles
    assert "data" in particles
    assert particles["type"] == "Sphere"


@pytest.mark.parametrize(
    "invalid_room_id,expected_error",
    [
        ("room:with:colons", "cannot contain ':'"),
        ("room with spaces", "invalid characters"),
        ("room@special", "invalid characters"),
        ("room#hash", "invalid characters"),
    ],
)
def test_create_room_validates_room_id(redis_client, invalid_room_id, expected_error):
    """Test room ID validation rejects invalid characters."""
    service = RoomService(redis_client)

    with pytest.raises(ValueError, match=expected_error):
        service.create_room(invalid_room_id, "user1")


@pytest.mark.parametrize(
    "valid_room_id",
    [
        "simple",
        "room-with-dashes",
        "room_with_underscores",
        "Room123",
        "UPPERCASE",
        "MixedCase123",
        "room.with.dots",
        "file.xyz",
        "data.structure.h5",
    ],
)
def test_create_room_accepts_valid_room_ids(redis_client, valid_room_id):
    """Test room ID validation accepts valid characters."""
    service = RoomService(redis_client)
    result = service.create_room(valid_room_id, "user1")

    assert result["created"] is True


def test_create_room_from_copy(redis_client):
    """Test creating room by copying from existing room."""
    service = RoomService(redis_client)

    # Create source room with data
    redis_client.zadd("room:source:trajectory:indices", {"frame0": 0, "frame1": 1})
    redis_client.hset("room:source:geometries", "test", "data")
    redis_client.hset("room:source:bookmarks", "0", "bookmark")

    # Copy room
    result = service.create_room("copy", "user1", copy_from="source")

    assert result["created"] is True
    assert result["frameCount"] == 2
    assert redis_client.exists("room:copy:trajectory:indices") > 0
    assert redis_client.exists("room:copy:geometries") > 0
    assert redis_client.exists("room:copy:bookmarks") > 0


def test_create_room_from_copy_shares_frame_data(redis_client):
    """Test copying room shares frame data (doesn't duplicate)."""
    service = RoomService(redis_client)

    # Create source room
    redis_client.zadd("room:source:trajectory:indices", {"frame0": 0, "frame1": 1})

    # Copy room
    service.create_room("copy", "user1", copy_from="source")

    # Verify indices are copied
    source_indices = redis_client.zrange("room:source:trajectory:indices", 0, -1)
    copy_indices = redis_client.zrange("room:copy:trajectory:indices", 0, -1)
    assert source_indices == copy_indices


def test_create_room_from_copy_initializes_metadata(redis_client):
    """Test copying room initializes metadata for new room."""
    service = RoomService(redis_client)

    # Create source room
    redis_client.zadd("room:source:trajectory:indices", {"frame0": 0})

    # Copy room
    service.create_room("copy", "user1", copy_from="source")

    # Verify new room has its own metadata
    assert redis_client.get("room:copy:current_frame") == "0"
    assert redis_client.get("room:copy:locked") == "0"
    assert redis_client.get("room:copy:hidden") == "0"


def test_create_room_from_nonexistent_source_raises_error(redis_client):
    """Test copying from nonexistent room raises ValueError."""
    service = RoomService(redis_client)

    with pytest.raises(ValueError, match="Source room 'nonexistent' not found"):
        service.create_room("test", "user1", copy_from="nonexistent")


def test_get_frame_count(redis_client):
    """Test getting frame count from room."""
    service = RoomService(redis_client)
    redis_client.zadd(
        "room:test:trajectory:indices", {"frame0": 0, "frame1": 1, "frame2": 2}
    )

    assert service.get_frame_count("test") == 3


def test_get_frame_count_empty_room(redis_client):
    """Test getting frame count from room with no frames."""
    service = RoomService(redis_client)
    assert service.get_frame_count("test") == 0


def test_get_frame_count_nonexistent_room(redis_client):
    """Test getting frame count from nonexistent room returns 0."""
    service = RoomService(redis_client)
    assert service.get_frame_count("nonexistent") == 0


@pytest.mark.parametrize(
    "redis_value,expected",
    [
        ("0", 0),
        ("5", 5),
        ("100", 100),
        (None, 0),
        ("-1", 0),
        ("-100", 0),
        ("invalid", 0),
        ("", 0),
    ],
)
def test_get_current_frame_handles_values(redis_client, redis_value, expected):
    """Test get_current_frame handles various valid and invalid values."""
    service = RoomService(redis_client)

    if redis_value is not None:
        redis_client.set("room:test:current_frame", redis_value)

    assert service.get_current_frame("test") == expected


def test_get_current_frame_positive_values(redis_client):
    """Test get_current_frame returns positive frame numbers correctly."""
    service = RoomService(redis_client)

    for value in [0, 1, 10, 100, 1000]:
        redis_client.set("room:test:current_frame", value)
        assert service.get_current_frame("test") == value


def test_create_room_uses_pipeline(redis_client, monkeypatch):
    """Test room creation uses Redis pipeline for atomic operations."""
    service = RoomService(redis_client)

    # Track pipeline usage
    pipeline_created = []
    original_pipeline = redis_client.pipeline

    def mock_pipeline(*args, **kwargs):
        pipe = original_pipeline(*args, **kwargs)
        pipeline_created.append(pipe)
        return pipe

    monkeypatch.setattr(redis_client, "pipeline", mock_pipeline)

    service.create_room("test", "user1")

    # Verify pipeline was used
    assert len(pipeline_created) > 0


def test_create_room_from_copy_uses_pipeline(redis_client, monkeypatch):
    """Test room copy uses Redis pipeline for efficient copying."""
    service = RoomService(redis_client)

    # Create source room
    redis_client.zadd("room:source:trajectory:indices", {"frame0": 0})

    # Track pipeline usage
    pipeline_created = []
    original_pipeline = redis_client.pipeline

    def mock_pipeline(*args, **kwargs):
        pipe = original_pipeline(*args, **kwargs)
        pipeline_created.append(pipe)
        return pipe

    monkeypatch.setattr(redis_client, "pipeline", mock_pipeline)

    service.create_room("copy", "user1", copy_from="source")

    # Verify pipeline was used
    assert len(pipeline_created) > 0
