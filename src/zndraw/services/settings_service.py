"""User settings management per room.

This service handles user-specific settings in rooms.
Settings are stored per user per room to allow personalized configurations.
"""

import json
import logging

from redis import Redis

log = logging.getLogger(__name__)


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

    def initialize_defaults(self, room_id: str, user_name: str) -> None:
        """Initialize default settings for a user in a room.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            Username to initialize settings for
        """
        from zndraw.settings import RoomConfig

        default_config = RoomConfig()
        config_dict = default_config.model_dump()

        # Store each settings category
        settings_key = f"room:{room_id}:user:{user_name}:settings"
        for category_name, category_data in config_dict.items():
            self.r.hset(settings_key, category_name, json.dumps(category_data))

        log.info(
            f"Initialized default settings for user '{user_name}' in room '{room_id}'"
        )

    def get_all(self, room_id: str, user_name: str) -> dict:
        """Get all settings for a user in a room.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            Username

        Returns
        -------
        dict
            Dictionary of settings categories
        """
        settings_key = f"room:{room_id}:user:{user_name}:settings"
        raw_settings: dict = self.r.hgetall(settings_key)  # type: ignore
        return {key: json.loads(value) for key, value in raw_settings.items()}

    def get_category(self, room_id: str, user_name: str, category: str) -> dict | None:
        """Get specific settings category for a user.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            Username
        category : str
            Settings category name

        Returns
        -------
        dict | None
            Settings data for the category, or None if not found
        """
        settings_key = f"room:{room_id}:user:{user_name}:settings"
        value: str | None = self.r.hget(settings_key, category)  # type: ignore
        return json.loads(value) if value else None

    def update_category(
        self, room_id: str, user_name: str, category: str, data: dict
    ) -> None:
        """Update specific settings category for a user.

        Parameters
        ----------
        room_id : str
            Room identifier
        user_name : str
            Username
        category : str
            Settings category name
        data : dict
            Settings data to store
        """
        settings_key = f"room:{room_id}:user:{user_name}:settings"
        self.r.hset(settings_key, category, json.dumps(data))
        log.debug(
            f"Updated settings category '{category}' for user '{user_name}' in room '{room_id}'"
        )
