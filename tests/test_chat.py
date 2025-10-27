import time

import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw


def test_rest_get_chat_messages_empty(server):
    """Test fetching messages from an empty room"""
    room = "test-chat-room"
    # Create room first
    response = requests.post(
        f"{server}/api/rooms/{room}/join", json={}, headers=get_jwt_auth_headers(server)
    )
    assert response.status_code == 200

    # Fetch chat messages
    response = requests.get(f"{server}/api/rooms/{room}/chat/messages")
    assert response.status_code == 200
    data = response.json()
    assert data["messages"] == []
    assert data["metadata"]["totalCount"] == 0
    assert data["metadata"]["hasMore"] is False
    assert data["metadata"]["oldestTimestamp"] is None
    assert data["metadata"]["newestTimestamp"] is None


def test_chat_message_create(server):
    """Test creating a chat message via socket"""
    room = "test-chat-room"
    user = "test-user"

    vis = ZnDraw(url=server, room=room, user=user)

    # Send a message
    response = vis.log("Hello **world**!")
    assert response["success"] is True
    assert "message" in response
    message = response["message"]
    assert message["content"] == "Hello **world**!"
    assert message["author"]["id"] == user
    assert message["roomId"] == room
    assert message["isEdited"] is False


def test_chat_message_edit(server):
    """Test editing a chat message"""
    room = "test-chat-room"
    user = "test-user"

    vis = ZnDraw(url=server, room=room, user=user)

    # Send a message
    response = vis.log("Original message")
    assert response["success"] is True
    message_id = response["message"]["id"]

    # Edit the message
    response = vis.edit_message(message_id, "Edited message")
    assert response["success"] is True
    assert response["message"]["content"] == "Edited message"
    assert response["message"]["isEdited"] is True
    assert response["message"]["id"] == message_id


def test_chat_message_edit_unauthorized(server):
    """Test that users cannot edit other users' messages"""
    room = "test-chat-room"
    user1 = "user1"
    user2 = "user2"

    # User 1 sends a message
    vis1 = ZnDraw(url=server, room=room, user=user1)
    response = vis1.log("User 1 message")
    assert response["success"] is True
    message_id = response["message"]["id"]

    # User 2 tries to edit user 1's message
    vis2 = ZnDraw(url=server, room=room, user=user2)
    response = vis2.edit_message(message_id, "Hacked message")
    assert response["success"] is False
    assert "error" in response


def test_get_chat_messages_rest(server):
    """Test REST API for fetching messages"""
    room = "test-chat-room"
    user = "test-user"

    vis = ZnDraw(url=server, room=room, user=user)

    # Send multiple messages
    for i in range(5):
        vis.log(f"Message {i + 1}")
        time.sleep(0.01)  # Small delay to ensure different timestamps

    # Fetch messages
    response = vis.get_messages(limit=10)
    assert len(response["messages"]) == 5
    assert response["metadata"]["totalCount"] == 5
    assert response["metadata"]["hasMore"] is False

    # Check order (newest first)
    assert response["messages"][0]["content"] == "Message 5"
    assert response["messages"][4]["content"] == "Message 1"


def test_chat_pagination(server):
    """Test pagination with before/after parameters"""
    room = "test-chat-room"
    user = "test-user"

    vis = ZnDraw(url=server, room=room, user=user)

    # Send 10 messages
    for i in range(10):
        vis.log(f"Message {i + 1}")
        time.sleep(0.01)

    # Fetch first page (limit 5)
    response = vis.get_messages(limit=5)
    assert len(response["messages"]) == 5
    assert response["metadata"]["hasMore"] is True
    oldest_timestamp = response["metadata"]["oldestTimestamp"]

    # Fetch next page using before parameter
    response = vis.get_messages(limit=5, before=oldest_timestamp)
    assert len(response["messages"]) == 5
    assert response["metadata"]["hasMore"] is False


def test_chat_room_isolation(server):
    """Test messages are isolated per room"""
    room1 = "room1"
    room2 = "room2"
    user = "test-user"

    vis1 = ZnDraw(url=server, room=room1, user=user)
    vis1.log("Message in room 1")

    vis2 = ZnDraw(url=server, room=room2, user=user)
    vis2.log("Message in room 2")

    # Check room 1 messages
    response1 = vis1.get_messages()
    assert len(response1["messages"]) == 1
    assert response1["messages"][0]["content"] == "Message in room 1"

    # Check room 2 messages
    response2 = vis2.get_messages()
    assert len(response2["messages"]) == 1
    assert response2["messages"][0]["content"] == "Message in room 2"


def test_chat_python_client(server):
    """Test vis.log() method"""
    room = "test-chat-room"
    user = "test-user"

    vis = ZnDraw(url=server, room=room, user=user)

    # Test sending a message
    response = vis.log("Test message with **markdown**")
    assert response["success"] is True
    assert "message" in response

    # Fetch and verify
    messages = vis.get_messages()
    assert len(messages["messages"]) == 1
    assert messages["messages"][0]["content"] == "Test message with **markdown**"
