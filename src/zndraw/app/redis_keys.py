"""Redis key management for ZnDraw."""

from dataclasses import dataclass


@dataclass
class ExtensionKeys:
    """Redis keys for an extension.

    Centralizes key construction to avoid duplication and errors.
    """

    schema: str
    idle_workers: str
    progressing_workers: str
    queue: str

    @classmethod
    def for_extension(
        cls, room_id: str, category: str, extension: str
    ) -> "ExtensionKeys":
        """Create ExtensionKeys for a specific extension.

        Args:
            room_id: The room identifier
            category: The extension category (e.g., 'modifiers', 'selections')
            extension: The extension name

        Returns:
            ExtensionKeys instance with all relevant Redis keys
        """
        base = f"room:{room_id}:extensions:{category}"
        return cls(
            schema=base,
            idle_workers=f"{base}:{extension}:idle_workers",
            progressing_workers=f"{base}:{extension}:progressing_workers",
            queue=f"{base}:{extension}:queue",
        )

    @staticmethod
    def schema_key(room_id: str, category: str) -> str:
        """Get the schema hash key for a room and category.

        This is the Redis hash that stores all extension schemas for a category.

        Args:
            room_id: The room identifier
            category: The extension category

        Returns:
            Redis key for the schema hash
        """
        return f"room:{room_id}:extensions:{category}"

    @staticmethod
    def user_extensions_key(room_id: str, category: str, sid: str) -> str:
        """Get the reverse mapping key for a worker's registered extensions.

        This key maps from a session ID to the set of extension names
        that the worker can handle.

        Args:
            room_id: The room identifier
            category: The extension category
            sid: The session ID of the worker

        Returns:
            Redis key for the user's extensions set
        """
        return f"room:{room_id}:extensions:{category}:{sid}"
