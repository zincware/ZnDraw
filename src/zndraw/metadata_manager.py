"""Client-side interface for room metadata."""

import typing as t
from collections.abc import MutableMapping

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class RoomMetadata(MutableMapping):
    """Client-side interface for room metadata.

    Provides dict-like access to room metadata with automatic
    synchronization to the server.

    Parameters
    ----------
    zndraw_instance : ZnDraw
        The ZnDraw instance this metadata manager belongs to.

    Examples
    --------
    >>> vis.metadata["file"] = "data.xyz"
    >>> file = vis.metadata["file"]
    >>> del vis.metadata["file"]
    >>> for key in vis.metadata:
    ...     print(key, vis.metadata[key])

    Notes
    -----
    All metadata values must be strings (Redis hash requirement).
    Metadata is lazily loaded from the server on first access.
    Write operations require an active connection to the server.
    """

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        """Initialize the metadata manager.

        Parameters
        ----------
        zndraw_instance : ZnDraw
            The ZnDraw instance this metadata manager belongs to.
        """
        self.vis = zndraw_instance
        self._cache: dict[str, str] = {}
        self._loaded = False

    def _ensure_loaded(self) -> None:
        """Lazy load metadata from server."""
        if not self._loaded:
            self._cache = self.vis.api.get_metadata()
            self._loaded = True

    def __getitem__(self, key: str) -> str:
        """Get metadata value by key.

        Parameters
        ----------
        key : str
            The metadata field name.

        Returns
        -------
        str
            The metadata value.

        Raises
        ------
        KeyError
            If the key does not exist.
        """
        self._ensure_loaded()
        return self._cache[key]

    def __setitem__(self, key: str, value: str) -> None:
        """Set metadata field.

        Parameters
        ----------
        key : str
            The metadata field name.
        value : str
            The metadata value (must be string).

        Raises
        ------
        RuntimeError
            If client is not connected to the server.
        TypeError
            If value is not a string.

        Notes
        -----
        Requires connection. The backend will check if the room is locked.
        """
        if not isinstance(value, str):
            raise TypeError(
                f"Metadata values must be strings, got {type(value).__name__}"
            )

        if not self.vis.socket.sio.connected:
            raise RuntimeError("Client is not connected.")

        self.vis.api.set_metadata({key: value})
        self._cache[key] = value

    def __delitem__(self, key: str) -> None:
        """Delete metadata field.

        Parameters
        ----------
        key : str
            The metadata field name to delete.

        Raises
        ------
        RuntimeError
            If client is not connected to the server.
        KeyError
            If the key does not exist.

        Notes
        -----
        Requires connection. The backend will check if the room is locked.
        """
        if not self.vis.socket.sio.connected:
            raise RuntimeError("Client is not connected.")

        self.vis.api.delete_metadata_field(key)
        if key in self._cache:
            del self._cache[key]

    def __iter__(self):
        """Iterate over metadata keys.

        Yields
        ------
        str
            Metadata field names.
        """
        self._ensure_loaded()
        return iter(self._cache)

    def __len__(self) -> int:
        """Get number of metadata fields.

        Returns
        -------
        int
            Number of metadata fields.
        """
        self._ensure_loaded()
        return len(self._cache)

    def __repr__(self) -> str:
        """String representation of metadata.

        Returns
        -------
        str
            String representation showing all metadata fields.
        """
        self._ensure_loaded()
        return f"RoomMetadata({self._cache!r})"

    def refresh(self) -> None:
        """Refresh metadata from server.

        Forces a reload of metadata from the server, discarding
        the local cache.
        """
        self._loaded = False
        self._ensure_loaded()
