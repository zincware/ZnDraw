"""Session settings management per room.

This service handles session-specific settings in rooms.
Settings are stored per session per room to allow different configurations
per browser window/tab.
"""

import json

from redis import Redis

from zndraw.settings import RoomConfig


class SettingsService:
    """Handles session settings in rooms.

    Settings are stored per session per room, allowing each browser
    window/tab to have its own configuration.

    Parameters
    ----------
    redis_client : Redis
        Redis client instance for data operations
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def _get_settings_key(self, room_id: str, session_id: str) -> str:
        return f"room:{room_id}:session_settings:{session_id}"

    def get_all(self, room_id: str, session_id: str) -> dict[str, dict]:
        """Get all settings categories for a session.

        Parameters
        ----------
        room_id : str
            Room identifier
        session_id : str
            Session identifier

        Returns
        -------
        dict[str, dict]
            All settings data with defaults applied for missing categories
        """
        settings_key = self._get_settings_key(room_id, session_id)
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

    def update_all(self, room_id: str, session_id: str, data: dict[str, dict]) -> None:
        """Update multiple settings categories for a session.

        Parameters
        ----------
        room_id : str
            Room identifier
        session_id : str
            Session identifier
        data : dict[str, dict]
            Settings data to store, keyed by category name
        """
        settings_key = self._get_settings_key(room_id, session_id)
        mapping = {category: json.dumps(values) for category, values in data.items()}
        if mapping:
            self.r.hset(settings_key, mapping=mapping)
