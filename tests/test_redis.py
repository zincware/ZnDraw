"""Unit tests for Redis key utilities."""

import pytest

from zndraw.redis import RedisKey


class TestRedisKey:
    """Tests for RedisKey patterns and parsing."""

    def test_presence_sid_roundtrip(self) -> None:
        """Builder and parser must stay in sync."""
        room_id = "test-room-42"
        sid = "abc123XYZ"

        # Build key
        key = RedisKey.presence_sid(room_id, sid)

        # Parse should extract the same SID
        parsed_sid = RedisKey.parse_presence_sid(key)
        assert parsed_sid == sid

    def test_presence_sid_pattern_matches_key(self) -> None:
        """Pattern should match keys created by presence_sid."""
        room_id = "test-room-42"
        sid = "test-sid"

        key = RedisKey.presence_sid(room_id, sid)
        pattern = RedisKey.presence_sid_pattern(room_id)

        # Pattern is presence:room:42:sid:* - key should start with prefix
        prefix = pattern.rstrip("*")
        assert key.startswith(prefix)

    @pytest.mark.parametrize(
        "invalid_key",
        [
            "",
            "presence",
            "presence:room:1:sid",  # Missing SID
            "presence:room:1:user:123",  # Wrong format (old style)
            "presence:WRONG:1:sid:abc",  # Invalid segment
            "other:room:1:sid:abc",  # Wrong prefix
            "presence:room:1:sid:abc:extra",  # Too many parts
        ],
    )
    def test_parse_presence_sid_rejects_invalid(self, invalid_key: str) -> None:
        """Parser should return None for invalid keys."""
        assert RedisKey.parse_presence_sid(invalid_key) is None


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
