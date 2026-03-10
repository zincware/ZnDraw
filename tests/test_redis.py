"""Unit tests for Redis key utilities."""

from zndraw.redis import RedisKey


class TestEditLockKey:
    """Tests for edit lock Redis key."""

    def test_edit_lock_key(self) -> None:
        """edit_lock returns correct key format."""
        assert RedisKey.edit_lock("test-room-42") == "room:test-room-42:edit-lock"

    def test_edit_lock_key_different_rooms(self) -> None:
        """Different rooms produce different edit lock keys."""
        key1 = RedisKey.edit_lock("room-1")
        key2 = RedisKey.edit_lock("room-2")
        assert key1 != key2
        assert "room-1" in key1
        assert "room-2" in key2
