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

    @classmethod
    def for_global_extension(cls, category: str, extension: str) -> "ExtensionKeys":
        """Create ExtensionKeys for a global extension (not room-specific).

        Args:
            category: The extension category (e.g., 'modifiers', 'selections')
            extension: The extension name

        Returns:
            ExtensionKeys instance with all relevant Redis keys for global extensions
        """
        base = f"global:extensions:{category}"
        return cls(
            schema=base,
            idle_workers=f"{base}:{extension}:idle_workers",
            progressing_workers=f"{base}:{extension}:progressing_workers",
            queue=f"{base}:{extension}:queue",
        )

    @staticmethod
    def global_schema_key(category: str) -> str:
        """Get the schema hash key for global extensions in a category.

        This is the Redis hash that stores all global extension schemas for a category.

        Args:
            category: The extension category

        Returns:
            Redis key for the global schema hash
        """
        return f"global:extensions:{category}"

    @staticmethod
    def global_user_extensions_key(category: str, sid: str) -> str:
        """Get the reverse mapping key for a worker's registered global extensions.

        This key maps from a session ID to the set of global extension names
        that the worker can handle.

        Args:
            category: The extension category
            sid: The session ID of the worker

        Returns:
            Redis key for the user's global extensions set
        """
        return f"global:extensions:{category}:{sid}"


@dataclass
class FilesystemKeys:
    """Redis keys for a filesystem provider.

    Centralizes key construction for filesystem-related Redis keys.
    """

    metadata: str
    worker: str

    @classmethod
    def for_filesystem(cls, room_id: str, fs_name: str) -> "FilesystemKeys":
        """Create FilesystemKeys for a specific filesystem.

        Parameters
        ----------
        room_id : str
            The room identifier
        fs_name : str
            The filesystem name

        Returns
        -------
        FilesystemKeys
            Instance with all relevant Redis keys
        """
        base = f"room:{room_id}:filesystems:{fs_name}"
        return cls(
            metadata=base,
            worker=f"{base}:worker",
        )

    @classmethod
    def for_global_filesystem(cls, fs_name: str) -> "FilesystemKeys":
        """Create FilesystemKeys for a global filesystem (not room-specific).

        Parameters
        ----------
        fs_name : str
            The filesystem name

        Returns
        -------
        FilesystemKeys
            Instance with all relevant Redis keys for global filesystems
        """
        base = f"global:filesystems:{fs_name}"
        return cls(
            metadata=base,
            worker=f"{base}:worker",
        )

    @staticmethod
    def user_filesystems_key(room_id: str, sid: str) -> str:
        """Get the reverse mapping key for a worker's registered filesystems.

        This key maps from a session ID to the set of filesystem names
        that the worker provides.

        Parameters
        ----------
        room_id : str
            The room identifier
        sid : str
            The session ID of the worker

        Returns
        -------
        str
            Redis key for the user's filesystems set
        """
        return f"room:{room_id}:filesystems:{sid}"

    @staticmethod
    def global_user_filesystems_key(sid: str) -> str:
        """Get the reverse mapping key for a worker's registered global filesystems.

        This key maps from a session ID to the set of global filesystem names
        that the worker provides.

        Parameters
        ----------
        sid : str
            The session ID of the worker

        Returns
        -------
        str
            Redis key for the user's global filesystems set
        """
        return f"global:filesystems:{sid}"
