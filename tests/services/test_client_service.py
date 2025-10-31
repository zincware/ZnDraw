"""Tests for ClientService.

Tests follow CLAUDE.md principles:
- Each test is a function (not a method)
- Tests are specific and test only one thing
- Uses pytest.mark.parametrize to avoid duplication
"""

import pytest

from zndraw.services.client_service import ClientService


def test_update_user_room(redis_client):
    """Test updating user's current room."""
    service = ClientService(redis_client)
    service.update_user_room("alice", "physics-lab")

    assert redis_client.hget("user:alice", "currentRoom") == "physics-lab"


def test_update_user_room_creates_user_hash(redis_client):
    """Test updating user room creates user hash if it doesn't exist."""
    service = ClientService(redis_client)
    service.update_user_room("new-user", "room1")

    assert redis_client.exists("user:new-user") > 0
    assert redis_client.hget("user:new-user", "currentRoom") == "room1"


def test_update_user_room_overwrites_previous(redis_client):
    """Test updating user room overwrites previous value."""
    service = ClientService(redis_client)
    service.update_user_room("bob", "room1")
    service.update_user_room("bob", "room2")

    assert redis_client.hget("user:bob", "currentRoom") == "room2"


def test_update_user_room_multiple_users(redis_client):
    """Test updating rooms for multiple users independently."""
    service = ClientService(redis_client)
    service.update_user_room("alice", "room1")
    service.update_user_room("bob", "room2")
    service.update_user_room("charlie", "room1")

    assert redis_client.hget("user:alice", "currentRoom") == "room1"
    assert redis_client.hget("user:bob", "currentRoom") == "room2"
    assert redis_client.hget("user:charlie", "currentRoom") == "room1"


def test_update_user_room_empty_room_name(redis_client):
    """Test updating with empty room name."""
    service = ClientService(redis_client)
    service.update_user_room("alice", "")

    assert redis_client.hget("user:alice", "currentRoom") == ""


def test_update_user_room_special_characters(redis_client):
    """Test room names with special characters."""
    service = ClientService(redis_client)
    special_rooms = [
        "room-with-dashes",
        "room_with_underscores",
        "room.with.dots",
        "room:with:colons",
        "room/with/slashes",
    ]

    for i, room_name in enumerate(special_rooms):
        user_name = f"user{i}"
        service.update_user_room(user_name, room_name)
        assert redis_client.hget(f"user:{user_name}", "currentRoom") == room_name


def test_add_user_to_room(redis_client):
    """Test adding user to room's user set."""
    service = ClientService(redis_client)
    service.add_user_to_room("chemistry-lab", "alice")

    assert redis_client.sismember("room:chemistry-lab:users", "alice")


def test_add_user_to_room_creates_set(redis_client):
    """Test adding user creates room user set if it doesn't exist."""
    service = ClientService(redis_client)
    service.add_user_to_room("new-room", "alice")

    assert redis_client.exists("room:new-room:users") > 0
    assert redis_client.scard("room:new-room:users") == 1


def test_add_user_to_room_idempotent(redis_client):
    """Test adding same user twice is idempotent."""
    service = ClientService(redis_client)
    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room1", "alice")

    assert redis_client.scard("room:room1:users") == 1
    assert redis_client.sismember("room:room1:users", "alice")


def test_add_multiple_users_to_room(redis_client):
    """Test adding multiple users to same room."""
    service = ClientService(redis_client)
    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room1", "bob")
    service.add_user_to_room("room1", "charlie")

    assert redis_client.scard("room:room1:users") == 3
    assert redis_client.sismember("room:room1:users", "alice")
    assert redis_client.sismember("room:room1:users", "bob")
    assert redis_client.sismember("room:room1:users", "charlie")


def test_add_user_to_room_case_sensitive(redis_client):
    """Test that usernames are case-sensitive."""
    service = ClientService(redis_client)
    service.add_user_to_room("room1", "Alice")
    service.add_user_to_room("room1", "alice")

    # Should have 2 distinct users
    assert redis_client.scard("room:room1:users") == 2


def test_remove_user_from_room(redis_client):
    """Test removing user from room's user set."""
    service = ClientService(redis_client)
    service.add_user_to_room("room1", "alice")
    service.remove_user_from_room("room1", "alice")

    assert not redis_client.sismember("room:room1:users", "alice")


def test_remove_user_from_room_idempotent(redis_client):
    """Test removing nonexistent user is idempotent."""
    service = ClientService(redis_client)
    # Should not raise even if user not in set
    service.remove_user_from_room("room1", "alice")
    service.remove_user_from_room("room1", "alice")

    assert not redis_client.sismember("room:room1:users", "alice")


def test_remove_user_from_room_leaves_others(redis_client):
    """Test removing one user doesn't affect other users."""
    service = ClientService(redis_client)
    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room1", "bob")
    service.add_user_to_room("room1", "charlie")

    service.remove_user_from_room("room1", "bob")

    assert redis_client.scard("room:room1:users") == 2
    assert redis_client.sismember("room:room1:users", "alice")
    assert redis_client.sismember("room:room1:users", "charlie")
    assert not redis_client.sismember("room:room1:users", "bob")


def test_remove_user_from_nonexistent_room(redis_client):
    """Test removing user from room that doesn't exist."""
    service = ClientService(redis_client)
    # Should not raise error
    service.remove_user_from_room("nonexistent-room", "alice")


def test_get_room_users(redis_client):
    """Test getting all users in a room."""
    service = ClientService(redis_client)
    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room1", "bob")

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
    service.add_user_to_room("room1", "alice")

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

    for user_name in user_names:
        service.add_user_to_room("room1", user_name)

    users = service.get_room_users("room1")
    assert len(users) == len(user_names)
    for user_name in user_names:
        assert user_name in users


def test_update_user_and_room_membership(redis_client):
    """Test atomic update of user room and room membership."""
    service = ClientService(redis_client)
    service.update_user_and_room_membership("alice", "physics-lab")

    # Verify both operations completed
    assert redis_client.hget("user:alice", "currentRoom") == "physics-lab"
    assert redis_client.sismember("room:physics-lab:users", "alice")


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

    service.update_user_and_room_membership("alice", "room1")
    service.update_user_and_room_membership("alice", "room2")
    service.update_user_and_room_membership("alice", "room3")

    # User's current room should be room3 (last one)
    assert redis_client.hget("user:alice", "currentRoom") == "room3"

    # User should be in all room sets
    assert redis_client.sismember("room:room1:users", "alice")
    assert redis_client.sismember("room:room2:users", "alice")
    assert redis_client.sismember("room:room3:users", "alice")


def test_update_user_and_room_membership_same_room_twice(redis_client):
    """Test updating user to same room twice is idempotent."""
    service = ClientService(redis_client)

    service.update_user_and_room_membership("alice", "room1")
    service.update_user_and_room_membership("alice", "room1")

    assert redis_client.hget("user:alice", "currentRoom") == "room1"
    assert redis_client.scard("room:room1:users") == 1


def test_user_in_multiple_rooms(redis_client):
    """Test user can be tracked in multiple rooms simultaneously."""
    service = ClientService(redis_client)

    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room2", "alice")
    service.add_user_to_room("room3", "alice")

    assert redis_client.sismember("room:room1:users", "alice")
    assert redis_client.sismember("room:room2:users", "alice")
    assert redis_client.sismember("room:room3:users", "alice")


def test_multiple_users_multiple_rooms(redis_client):
    """Test complex scenario with multiple users in multiple rooms."""
    service = ClientService(redis_client)

    # Alice in room1 and room2
    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room2", "alice")

    # Bob in room2 and room3
    service.add_user_to_room("room2", "bob")
    service.add_user_to_room("room3", "bob")

    # Charlie in all three rooms
    service.add_user_to_room("room1", "charlie")
    service.add_user_to_room("room2", "charlie")
    service.add_user_to_room("room3", "charlie")

    # Verify room1 has alice and charlie
    assert service.get_room_users("room1") == {"alice", "charlie"}

    # Verify room2 has all three users
    assert service.get_room_users("room2") == {"alice", "bob", "charlie"}

    # Verify room3 has bob and charlie
    assert service.get_room_users("room3") == {"bob", "charlie"}


def test_room_isolation(redis_client):
    """Test that different rooms maintain separate user sets."""
    service = ClientService(redis_client)

    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room2", "bob")

    # Rooms should be completely isolated
    assert service.get_room_users("room1") == {"alice"}
    assert service.get_room_users("room2") == {"bob"}


def test_redis_key_format_user(redis_client):
    """Test that user data uses correct Redis key format."""
    service = ClientService(redis_client)
    service.update_user_room("alice", "room1")

    # Should use user:<userName> format
    assert redis_client.exists("user:alice")
    assert not redis_client.exists("client:alice")


def test_redis_key_format_room(redis_client):
    """Test that room user sets use correct Redis key format."""
    service = ClientService(redis_client)
    service.add_user_to_room("room1", "alice")

    # Should use room:<roomName>:users format
    assert redis_client.exists("room:room1:users")


def test_concurrent_room_updates(redis_client):
    """Test that concurrent updates to different rooms don't interfere."""
    service = ClientService(redis_client)

    # Simulate concurrent operations
    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room2", "bob")
    service.add_user_to_room("room1", "charlie")
    service.add_user_to_room("room2", "david")

    # All operations should succeed independently
    assert service.get_room_users("room1") == {"alice", "charlie"}
    assert service.get_room_users("room2") == {"bob", "david"}


def test_user_room_consistency(redis_client):
    """Test consistency between user's currentRoom and room membership."""
    service = ClientService(redis_client)

    # Use atomic update
    service.update_user_and_room_membership("alice", "room1")

    # Both should be consistent
    current_room = redis_client.hget("user:alice", "currentRoom")
    is_member = redis_client.sismember(f"room:{current_room}:users", "alice")

    assert current_room == "room1"
    assert is_member


def test_remove_all_users_from_room(redis_client):
    """Test removing all users from a room."""
    service = ClientService(redis_client)

    # Add multiple users
    service.add_user_to_room("room1", "alice")
    service.add_user_to_room("room1", "bob")
    service.add_user_to_room("room1", "charlie")

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

    service.add_user_to_room("room1", user_name)
    assert redis_client.sismember("room:room1:users", user_name)

    service.update_user_room(user_name, "room1")
    assert redis_client.hget(f"user:{user_name}", "currentRoom") == "room1"

    users = service.get_room_users("room1")
    assert user_name in users

    service.remove_user_from_room("room1", user_name)
    assert not redis_client.sismember("room:room1:users", user_name)
