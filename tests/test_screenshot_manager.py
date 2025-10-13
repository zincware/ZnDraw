"""Tests for screenshot management functionality."""

import tempfile

import pytest

from zndraw.screenshot_manager import Screenshot, ScreenshotManager


@pytest.fixture
def temp_storage():
    """Create a temporary storage directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


def test_screenshot_manager_initialization(temp_storage):
    """Test ScreenshotManager initializes correctly."""
    manager = ScreenshotManager("test_room", temp_storage)

    assert manager.room_id == "test_room"
    assert manager.screenshot_dir.exists()
    assert manager.screenshot_dir.is_dir()


def test_save_screenshot(temp_storage):
    """Test saving a screenshot."""
    manager = ScreenshotManager("test_room", temp_storage)

    image_data = b"fake image data"

    screenshot = manager.save(image_data, "png", width=1920, height=1080)

    assert screenshot.format == "png"
    assert screenshot.size == len(image_data)
    assert screenshot.width == 1920
    assert screenshot.height == 1080
    assert screenshot.id == 1

    screenshot_files = list(manager.screenshot_dir.glob("*.png"))
    assert len(screenshot_files) == 1

    saved_file = screenshot_files[0]
    assert saved_file.read_bytes() == image_data


def test_sequential_ids(temp_storage):
    """Test IDs increment sequentially."""
    manager = ScreenshotManager("test_room", temp_storage)

    s1 = manager.save(b"image 1", "png")
    s2 = manager.save(b"image 2", "png")
    s3 = manager.save(b"image 3", "png")

    assert s1.id == 1
    assert s2.id == 2
    assert s3.id == 3


def test_invalid_format(temp_storage):
    """Test invalid format raises ValueError."""
    manager = ScreenshotManager("test_room", temp_storage)

    with pytest.raises(ValueError, match="Invalid format"):
        manager.save(b"test", "gif")


def test_get_screenshot(temp_storage):
    """Test retrieving a screenshot by ID."""
    manager = ScreenshotManager("test_room", temp_storage)

    image_data = b"test image"
    screenshot = manager.save(image_data, "png")

    result = manager.get(screenshot.id)

    assert result is not None
    filepath, metadata = result
    assert filepath.exists()
    assert metadata.id == screenshot.id
    assert metadata.format == "png"


def test_get_nonexistent_screenshot(temp_storage):
    """Test getting a nonexistent screenshot returns None."""
    manager = ScreenshotManager("test_room", temp_storage)

    result = manager.get(999)
    assert result is None


def test_list_screenshots(temp_storage):
    """Test listing screenshots."""
    manager = ScreenshotManager("test_room", temp_storage)

    for i in range(5):
        manager.save(f"image_{i}".encode(), "png")

    screenshots = manager.list(limit=10)
    assert len(screenshots) == 5

    screenshots = manager.list(limit=2)
    assert len(screenshots) == 2

    screenshots = manager.list(limit=2, offset=2)
    assert len(screenshots) == 2


def test_count_screenshots(temp_storage):
    """Test counting screenshots."""
    manager = ScreenshotManager("test_room", temp_storage)

    assert manager.count() == 0

    for i in range(3):
        manager.save(f"image_{i}".encode(), "png")

    assert manager.count() == 3


def test_delete_screenshot(temp_storage):
    """Test deleting a screenshot."""
    manager = ScreenshotManager("test_room", temp_storage)

    screenshot = manager.save(b"test image", "png")

    assert manager.count() == 1

    result = manager.delete(screenshot.id)
    assert result is True

    assert manager.count() == 0
    assert manager.get(screenshot.id) is None


def test_delete_nonexistent_screenshot(temp_storage):
    """Test deleting a nonexistent screenshot returns False."""
    manager = ScreenshotManager("test_room", temp_storage)

    result = manager.delete(999)
    assert result is False


def test_multiple_formats(temp_storage):
    """Test saving screenshots in different formats."""
    manager = ScreenshotManager("test_room", temp_storage)

    formats = ["png", "jpeg", "webp"]

    for fmt in formats:
        screenshot = manager.save(f"image_{fmt}".encode(), fmt)
        assert screenshot.format == fmt

        result = manager.get(screenshot.id)
        assert result is not None
        filepath, _ = result
        assert filepath.suffix == f".{fmt}"


def test_screenshot_ordering(temp_storage):
    """Test screenshots are listed by ID descending (newest first)."""
    manager = ScreenshotManager("test_room", temp_storage)

    screenshot_ids = []
    for i in range(3):
        screenshot = manager.save(f"image_{i}".encode(), "png")
        screenshot_ids.append(screenshot.id)

    screenshots = manager.list()

    assert screenshots[0].id == screenshot_ids[-1]
    assert screenshots[-1].id == screenshot_ids[0]
