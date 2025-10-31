import json
import time
import typing as t


def create_message(
    redis_client,
    room_id: str,
    author_user_name: str,
    content: str,
) -> dict:
    """Create and store a new chat message.

    Args:
        redis_client: Redis client instance
        room_id: Room identifier
        author_user_name: Username of the message author
        content: Message content (supports markdown)

    Returns:
        dict: The created message object
    """
    # Generate message ID using atomic counter
    counter_key = f"room:{room_id}:chat:counter"
    message_number = redis_client.incr(counter_key)
    message_id = f"msg_{room_id}_{message_number}"

    # Create message object
    timestamp = int(time.time() * 1000)  # Unix timestamp in milliseconds
    message = {
        "id": message_id,
        "roomId": room_id,
        "author": {"id": author_user_name},
        "content": content,
        "createdAt": timestamp,
        "updatedAt": timestamp,
        "isEdited": False,
    }

    # Store message data and index
    data_key = f"room:{room_id}:chat:data"
    index_key = f"room:{room_id}:chat:index"

    with redis_client.pipeline() as pipe:
        # Store message data as JSON
        pipe.hset(data_key, message_id, json.dumps(message))
        # Add to sorted set with timestamp as score
        pipe.zadd(index_key, {message_id: timestamp})
        pipe.execute()

    return message


def get_message(redis_client, room_id: str, message_id: str) -> t.Optional[dict]:
    """Fetch a message by ID.

    Args:
        redis_client: Redis client instance
        room_id: Room identifier
        message_id: Message identifier

    Returns:
        dict | None: The message object, or None if not found
    """
    data_key = f"room:{room_id}:chat:data"
    message_json = redis_client.hget(data_key, message_id)

    if message_json is None:
        return None

    return json.loads(message_json)


def update_message(
    redis_client, room_id: str, message_id: str, content: str
) -> t.Optional[dict]:
    """Update message content.

    Args:
        redis_client: Redis client instance
        room_id: Room identifier
        message_id: Message identifier
        content: New message content

    Returns:
        dict | None: The updated message object, or None if message not found
    """
    # Fetch existing message
    message = get_message(redis_client, room_id, message_id)

    if message is None:
        return None

    # Update message fields
    message["content"] = content
    message["updatedAt"] = int(time.time() * 1000)
    message["isEdited"] = True

    # Store updated message
    data_key = f"room:{room_id}:chat:data"
    redis_client.hset(data_key, message_id, json.dumps(message))

    return message
