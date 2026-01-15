"""Tests for ClientService.

Tests follow CLAUDE.md principles:
- Each test is a function (not a method)
- Tests are specific and test only one thing
- Uses pytest.mark.parametrize to avoid duplication
"""

import pytest

from zndraw.app.redis_keys import RoomKeys, UserKeys
from zndraw.services.client_service import ClientService


def test_remove_user_from_room(redis_client):
    """Test removing user from room's user set."""
    service = ClientService(redis_client)
    redis_client.sadd(RoomKeys("room1").users(), "alice")
    service.remove_user_from_room("room1", "alice")

    assert not redis_client.sismember(RoomKeys("room1").users(), "alice")


def test_remove_user_from_room_idempotent(redis_client):
    """Test removing nonexistent user is idempotent."""
    service = ClientService(redis_client)
    # Should not raise even if user not in set
    service.remove_user_from_room("room1", "alice")
    service.remove_user_from_room("room1", "alice")

    assert not redis_client.sismember(RoomKeys("room1").users(), "alice")


def test_remove_user_from_room_leaves_others(redis_client):
    """Test removing one user doesn't affect other users."""
    service = ClientService(redis_client)
    room_keys = RoomKeys("room1")
    redis_client.sadd(room_keys.users(), "alice")
    redis_client.sadd(room_keys.users(), "bob")
    redis_client.sadd(room_keys.users(), "charlie")

    service.remove_user_from_room("room1", "bob")

    assert redis_client.scard(room_keys.users()) == 2
    assert redis_client.sismember(room_keys.users(), "alice")
    assert redis_client.sismember(room_keys.users(), "charlie")
    assert not redis_client.sismember(room_keys.users(), "bob")


def test_remove_user_from_nonexistent_room(redis_client):
    """Test removing user from room that doesn't exist."""
    service = ClientService(redis_client)
    # Should not raise error
    service.remove_user_from_room("nonexistent-room", "alice")


def test_get_room_users(redis_client):
    """Test getting all users in a room."""
    service = ClientService(redis_client)
    room_keys = RoomKeys("room1")
    redis_client.sadd(room_keys.users(), "alice")
    redis_client.sadd(room_keys.users(), "bob")

    users = service.get_room_users("room1")
    assert "alice" in users
    assert "bob" in users
    assert len(users) == 2


def test_get_room_users_empty_room(redis_client):
    """Test getting users from room with no users."""
    service = ClientService(redis_client)

    users = service.get_room_users("empty-room")
    assert len(users) == 0
    assert isinstance(users, set)


def test_get_room_users_nonexistent_room(redis_client):
    """Test getting users from nonexistent room returns empty set."""
    service = ClientService(redis_client)

    users = service.get_room_users("nonexistent")
    assert len(users) == 0
    assert isinstance(users, set)


def test_get_room_users_returns_set(redis_client):
    """Test get_room_users returns a set."""
    service = ClientService(redis_client)
    redis_client.sadd(RoomKeys("room1").users(), "alice")

    users = service.get_room_users("room1")
    assert isinstance(users, set)


@pytest.mark.parametrize(
    "user_names",
    [
        ["alice"],
        ["alice", "bob"],
        ["alice", "bob", "charlie", "david", "eve"],
        ["user-1", "user-2", "user-3", "user-4", "user-5", "user-6", "user-7"],
    ],
)
def test_get_room_users_multiple_users(redis_client, user_names):
    """Test getting users works with various numbers of users."""
    service = ClientService(redis_client)
    room_keys = RoomKeys("room1")

    for user_name in user_names:
        redis_client.sadd(room_keys.users(), user_name)

    users = service.get_room_users("room1")
    assert len(users) == len(user_names)
    for user_name in user_names:
        assert user_name in users


def test_update_user_and_room_membership(redis_client):
    """Test atomic update of room membership and visit tracking."""
    service = ClientService(redis_client)
    keys = UserKeys("alice")
    service.update_user_and_room_membership("alice", "physics-lab")

    # Verify room membership and visit tracking
    assert redis_client.sismember(RoomKeys("physics-lab").users(), "alice")
    assert redis_client.sismember(keys.visited_rooms(), "physics-lab")


def test_update_user_and_room_membership_uses_pipeline(redis_client, monkeypatch):
    """Test update_user_and_room_membership uses pipeline for atomicity."""
    service = ClientService(redis_client)

    # Track pipeline usage
    pipeline_created = []
    original_pipeline = redis_client.pipeline

    def mock_pipeline(*args, **kwargs):
        pipe = original_pipeline(*args, **kwargs)
        pipeline_created.append(pipe)
        return pipe

    monkeypatch.setattr(redis_client, "pipeline", mock_pipeline)

    service.update_user_and_room_membership("alice", "room1")

    # Verify pipeline was used for atomicity
    assert len(pipeline_created) > 0


def test_update_user_and_room_membership_multiple_times(redis_client):
    """Test user can join multiple rooms sequentially."""
    service = ClientService(redis_client)
    keys = UserKeys("alice")

    service.update_user_and_room_membership("alice", "room1")
    service.update_user_and_room_membership("alice", "room2")
    service.update_user_and_room_membership("alice", "room3")

    # User should be in all room sets
    assert redis_client.sismember(RoomKeys("room1").users(), "alice")
    assert redis_client.sismember(RoomKeys("room2").users(), "alice")
    assert redis_client.sismember(RoomKeys("room3").users(), "alice")

    # User should have visited all rooms
    assert redis_client.sismember(keys.visited_rooms(), "room1")
    assert redis_client.sismember(keys.visited_rooms(), "room2")
    assert redis_client.sismember(keys.visited_rooms(), "room3")


def test_update_user_and_room_membership_same_room_twice(redis_client):
    """Test updating user to same room twice is idempotent."""
    service = ClientService(redis_client)

    service.update_user_and_room_membership("alice", "room1")
    service.update_user_and_room_membership("alice", "room1")

    # User should only appear once in room set (set semantics)
    assert redis_client.scard(RoomKeys("room1").users()) == 1


def test_user_in_multiple_rooms(redis_client):
    """Test user can be tracked in multiple rooms simultaneously."""
    service = ClientService(redis_client)

    service.update_user_and_room_membership("alice", "room1")
    service.update_user_and_room_membership("alice", "room2")
    service.update_user_and_room_membership("alice", "room3")

    assert redis_client.sismember(RoomKeys("room1").users(), "alice")
    assert redis_client.sismember(RoomKeys("room2").users(), "alice")
    assert redis_client.sismember(RoomKeys("room3").users(), "alice")


def test_multiple_users_multiple_rooms(redis_client):
    """Test complex scenario with multiple users in multiple rooms."""
    service = ClientService(redis_client)

    # Alice in room1 and room2
    redis_client.sadd(RoomKeys("room1").users(), "alice")
    redis_client.sadd(RoomKeys("room2").users(), "alice")

    # Bob in room2 and room3
    redis_client.sadd(RoomKeys("room2").users(), "bob")
    redis_client.sadd(RoomKeys("room3").users(), "bob")

    # Charlie in all three rooms
    redis_client.sadd(RoomKeys("room1").users(), "charlie")
    redis_client.sadd(RoomKeys("room2").users(), "charlie")
    redis_client.sadd(RoomKeys("room3").users(), "charlie")

    # Verify room1 has alice and charlie
    assert service.get_room_users("room1") == {"alice", "charlie"}

    # Verify room2 has all three users
    assert service.get_room_users("room2") == {"alice", "bob", "charlie"}

    # Verify room3 has bob and charlie
    assert service.get_room_users("room3") == {"bob", "charlie"}


def test_room_isolation(redis_client):
    """Test that different rooms maintain separate user sets."""
    service = ClientService(redis_client)

    redis_client.sadd(RoomKeys("room1").users(), "alice")
    redis_client.sadd(RoomKeys("room2").users(), "bob")

    # Rooms should be completely isolated
    assert service.get_room_users("room1") == {"alice"}
    assert service.get_room_users("room2") == {"bob"}


def test_redis_key_format_room(redis_client):
    """Test that room user sets use correct Redis key format."""
    service = ClientService(redis_client)
    service.update_user_and_room_membership("alice", "room1")

    # Should use room:<roomName>:users format
    assert redis_client.exists(RoomKeys("room1").users())


def test_concurrent_room_updates(redis_client):
    """Test that concurrent updates to different rooms don't interfere."""
    service = ClientService(redis_client)

    # Simulate concurrent operations
    redis_client.sadd(RoomKeys("room1").users(), "alice")
    redis_client.sadd(RoomKeys("room2").users(), "bob")
    redis_client.sadd(RoomKeys("room1").users(), "charlie")
    redis_client.sadd(RoomKeys("room2").users(), "david")

    # All operations should succeed independently
    assert service.get_room_users("room1") == {"alice", "charlie"}
    assert service.get_room_users("room2") == {"bob", "david"}


def test_remove_all_users_from_room(redis_client):
    """Test removing all users from a room."""
    service = ClientService(redis_client)
    room_keys = RoomKeys("room1")

    # Add multiple users
    redis_client.sadd(room_keys.users(), "alice")
    redis_client.sadd(room_keys.users(), "bob")
    redis_client.sadd(room_keys.users(), "charlie")

    # Remove all users
    service.remove_user_from_room("room1", "alice")
    service.remove_user_from_room("room1", "bob")
    service.remove_user_from_room("room1", "charlie")

    # Room should be empty
    assert len(service.get_room_users("room1")) == 0


@pytest.mark.parametrize(
    "user_name",
    [
        "user-with-dashes",
        "user_with_underscores",
        "user.with.dots",
        "user@email.com",
        "user123",
        "User",  # Capital letter
        "user-123_test.name@example",  # Complex combination
    ],
)
def test_special_characters_in_username(redis_client, user_name):
    """Test usernames with various special characters."""
    service = ClientService(redis_client)
    room_keys = RoomKeys("room1")

    # Add user to room via service
    service.update_user_and_room_membership(user_name, "room1")
    assert redis_client.sismember(room_keys.users(), user_name)

    # Verify get_room_users works with special characters
    users = service.get_room_users("room1")
    assert user_name in users

    # Verify remove works with special characters
    service.remove_user_from_room("room1", user_name)
    assert not redis_client.sismember(room_keys.users(), user_name)
