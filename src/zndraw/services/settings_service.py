"""User settings management per room.

This service handles user-specific settings in rooms.
Settings are stored per user per room to allow personalized configurations.
"""

import json

from redis import Redis

from zndraw.settings import RoomConfig


class SettingsService:
    """Handles user settings in rooms.

    Settings are stored per user per room, allowing each user
    to have their own personalized configuration in each room.

    Parameters
    ----------
    redis_client : Redis
        Redis client instance for data operations
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def _get_settings_key(self, room_id: str, user_name: str) -> str:
        return f"room:{room_id}:user:{user_name}:settings"

    def get_all(self, room_id: str, user_name: str) -> dict[str, dict]:
        """Get all settings categories for a user.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            Username

        Returns
        -------
        dict[str, dict]
            All settings data with defaults applied for missing categories
        """
        settings_key = self._get_settings_key(room_id, user_name)
        stored: dict[bytes, bytes] = self.r.hgetall(settings_key)  # type: ignore

        # Build result with defaults for each category (skip inherited callback field)
        result = {}
        for category, field_info in RoomConfig.model_fields.items():
            if field_info.default_factory is None:
                continue  # Skip non-settings fields like 'callback'
            # Redis may return string or bytes keys depending on decode_responses setting
            stored_value = stored.get(category) or stored.get(category.encode())
            if stored_value:
                # Handle both string and bytes values
                if isinstance(stored_value, bytes):
                    stored_value = stored_value.decode()
                result[category] = json.loads(stored_value)
            else:
                # Use default_factory to create default instance
                # by_alias=True ensures schema field names match (e.g., "camera" not "camera_type")
                result[category] = field_info.default_factory().model_dump(
                    by_alias=True
                )
        return result

    def update_all(self, room_id: str, user_name: str, data: dict[str, dict]) -> None:
        """Update multiple settings categories for a user.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            Username
        data : dict[str, dict]
            Settings data to store, keyed by category name
        """
        settings_key = self._get_settings_key(room_id, user_name)
        mapping = {category: json.dumps(values) for category, values in data.items()}
        if mapping:
            self.r.hset(settings_key, mapping=mapping)
