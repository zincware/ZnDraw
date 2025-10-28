"""Tests for ClientService.

Tests follow CLAUDE.md principles:
- Each test is a function (not a method)
- Tests are specific and test only one thing
- Uses pytest.mark.parametrize to avoid duplication
"""

import pytest

from zndraw.services.client_service import ClientService


def test_update_client_room(redis_client):
    """Test updating client's current room."""
    service = ClientService(redis_client)
    service.update_client_room("client1", "room1")

    assert redis_client.hget("client:client1", "currentRoom") == "room1"


def test_update_client_room_creates_client_hash(redis_client):
    """Test updating client room creates client hash if it doesn't exist."""
    service = ClientService(redis_client)
    service.update_client_room("newclient", "room1")

    assert redis_client.exists("client:newclient") > 0


def test_update_client_room_overwrites_previous(redis_client):
    """Test updating client room overwrites previous value."""
    service = ClientService(redis_client)
    service.update_client_room("client1", "room1")
    service.update_client_room("client1", "room2")

    assert redis_client.hget("client:client1", "currentRoom") == "room2"


def test_update_client_room_multiple_clients(redis_client):
    """Test updating rooms for multiple clients independently."""
    service = ClientService(redis_client)
    service.update_client_room("client1", "room1")
    service.update_client_room("client2", "room2")

    assert redis_client.hget("client:client1", "currentRoom") == "room1"
    assert redis_client.hget("client:client2", "currentRoom") == "room2"


def test_add_client_to_room(redis_client):
    """Test adding client to room's client set."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")

    assert redis_client.sismember("room:room1:clients", "client1")


def test_add_client_to_room_creates_set(redis_client):
    """Test adding client creates room client set if it doesn't exist."""
    service = ClientService(redis_client)
    service.add_client_to_room("newroom", "client1")

    assert redis_client.exists("room:newroom:clients") > 0


def test_add_client_to_room_idempotent(redis_client):
    """Test adding same client twice is idempotent."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room1", "client1")

    assert redis_client.scard("room:room1:clients") == 1


def test_add_multiple_clients_to_room(redis_client):
    """Test adding multiple clients to same room."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room1", "client2")
    service.add_client_to_room("room1", "client3")

    assert redis_client.scard("room:room1:clients") == 3


def test_remove_client_from_room(redis_client):
    """Test removing client from room's client set."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.remove_client_from_room("room1", "client1")

    assert not redis_client.sismember("room:room1:clients", "client1")


def test_remove_client_from_room_idempotent(redis_client):
    """Test removing nonexistent client is idempotent."""
    service = ClientService(redis_client)
    # Should not raise even if client not in set
    service.remove_client_from_room("room1", "client1")

    assert not redis_client.sismember("room:room1:clients", "client1")


def test_remove_client_from_room_leaves_others(redis_client):
    """Test removing one client doesn't affect other clients."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room1", "client2")
    service.add_client_to_room("room1", "client3")

    service.remove_client_from_room("room1", "client2")

    assert redis_client.scard("room:room1:clients") == 2
    assert redis_client.sismember("room:room1:clients", "client1")
    assert redis_client.sismember("room:room1:clients", "client3")


def test_get_room_clients(redis_client):
    """Test getting all clients in a room."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room1", "client2")

    clients = service.get_room_clients("room1")
    assert "client1" in clients
    assert "client2" in clients
    assert len(clients) == 2


def test_get_room_clients_empty_room(redis_client):
    """Test getting clients from room with no clients."""
    service = ClientService(redis_client)

    clients = service.get_room_clients("room1")
    assert len(clients) == 0
    assert isinstance(clients, set)


def test_get_room_clients_nonexistent_room(redis_client):
    """Test getting clients from nonexistent room returns empty set."""
    service = ClientService(redis_client)

    clients = service.get_room_clients("nonexistent")
    assert len(clients) == 0


def test_get_room_clients_returns_set(redis_client):
    """Test get_room_clients returns a set."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")

    clients = service.get_room_clients("room1")
    assert isinstance(clients, set)


@pytest.mark.parametrize(
    "client_ids",
    [
        ["client1"],
        ["client1", "client2"],
        ["client1", "client2", "client3", "client4", "client5"],
    ],
)
def test_get_room_clients_multiple_clients(redis_client, client_ids):
    """Test getting clients works with various numbers of clients."""
    service = ClientService(redis_client)

    for client_id in client_ids:
        service.add_client_to_room("room1", client_id)

    clients = service.get_room_clients("room1")
    assert len(clients) == len(client_ids)
    for client_id in client_ids:
        assert client_id in clients


def test_update_client_and_room_membership(redis_client):
    """Test atomic update of client room and room membership."""
    service = ClientService(redis_client)
    service.update_client_and_room_membership("client1", "room1")

    # Verify both operations completed
    assert redis_client.hget("client:client1", "currentRoom") == "room1"
    assert redis_client.sismember("room:room1:clients", "client1")


def test_update_client_and_room_membership_uses_pipeline(redis_client, monkeypatch):
    """Test update_client_and_room_membership uses pipeline for atomicity."""
    service = ClientService(redis_client)

    # Track pipeline usage
    pipeline_created = []
    original_pipeline = redis_client.pipeline

    def mock_pipeline(*args, **kwargs):
        pipe = original_pipeline(*args, **kwargs)
        pipeline_created.append(pipe)
        return pipe

    monkeypatch.setattr(redis_client, "pipeline", mock_pipeline)

    service.update_client_and_room_membership("client1", "room1")

    # Verify pipeline was used
    assert len(pipeline_created) > 0


def test_update_client_and_room_membership_multiple_times(redis_client):
    """Test client can join multiple rooms sequentially."""
    service = ClientService(redis_client)

    service.update_client_and_room_membership("client1", "room1")
    service.update_client_and_room_membership("client1", "room2")

    # Client's current room should be room2
    assert redis_client.hget("client:client1", "currentRoom") == "room2"
    # Client should be in both room sets
    assert redis_client.sismember("room:room1:clients", "client1")
    assert redis_client.sismember("room:room2:clients", "client1")


def test_client_in_multiple_rooms(redis_client):
    """Test client can be tracked in multiple rooms simultaneously."""
    service = ClientService(redis_client)

    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room2", "client1")
    service.add_client_to_room("room3", "client1")

    assert redis_client.sismember("room:room1:clients", "client1")
    assert redis_client.sismember("room:room2:clients", "client1")
    assert redis_client.sismember("room:room3:clients", "client1")


def test_multiple_clients_multiple_rooms(redis_client):
    """Test complex scenario with multiple clients in multiple rooms."""
    service = ClientService(redis_client)

    # Client1 in room1 and room2
    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room2", "client1")

    # Client2 in room2 and room3
    service.add_client_to_room("room2", "client2")
    service.add_client_to_room("room3", "client2")

    # Verify room1 has only client1
    assert service.get_room_clients("room1") == {"client1"}

    # Verify room2 has both clients
    assert service.get_room_clients("room2") == {"client1", "client2"}

    # Verify room3 has only client2
    assert service.get_room_clients("room3") == {"client2"}
