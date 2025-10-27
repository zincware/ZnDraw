import typing as t
from collections.abc import MutableMapping

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Bookmarks(MutableMapping):
    """Accessor for bookmarks using MutableMapping interface.

    Examples
    --------
    >>> vis.bookmarks[0] = "First Frame"
    >>> print(vis.bookmarks[0])
    >>> del vis.bookmarks[0]
    >>> len(vis.bookmarks)
    >>> list(vis.bookmarks)
    """

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def __getitem__(self, index: int) -> str:
        """Get bookmark label at a specific frame index.

        Parameters
        ----------
        index : int
            Frame index

        Returns
        -------
        str
            The bookmark label

        Raises
        ------
        KeyError
            If no bookmark exists at this index
        IndexError
            If index is out of range
        """
        if not isinstance(index, int):
            raise TypeError(
                f"Bookmark index must be an integer, not {type(index).__name__}"
            )

        # Check local cache first
        if index in self.vis._bookmarks:
            return self.vis._bookmarks[index]

        # Fetch from server
        response = self.vis.api.get_bookmark(index)
        label = response.get("label")
        if label is None:
            raise KeyError(f"Bookmark at index {index} does not exist")

        # Update cache
        self.vis._bookmarks[index] = label
        return label

    def __setitem__(self, index: int, label: str) -> None:
        """Set or update a bookmark at a specific frame index.

        Parameters
        ----------
        index : int
            Frame index
        label : str
            Bookmark label

        Raises
        ------
        TypeError
            If index is not an integer or label is not a string
        IndexError
            If index is out of range
        ValueError
            If label is empty
        """
        if not isinstance(index, int):
            raise TypeError(
                f"Bookmark index must be an integer, not {type(index).__name__}"
            )

        if not isinstance(label, str):
            raise TypeError("Bookmark label must be a string")

        if not label:
            raise ValueError("Bookmark label cannot be empty")

        # Validate index is in range
        if index < 0 or index >= len(self.vis):
            raise IndexError(
                f"Bookmark index {index} out of range (0-{len(self.vis) - 1})"
            )

        # Update via API
        self.vis.api.set_bookmark(index, label)

        # Update local cache
        self.vis._bookmarks[index] = label

    def __delitem__(self, index: int) -> None:
        """Delete a bookmark at a specific frame index.

        Parameters
        ----------
        index : int
            Frame index

        Raises
        ------
        KeyError
            If no bookmark exists at this index
        """
        if not isinstance(index, int):
            raise TypeError(
                f"Bookmark index must be an integer, not {type(index).__name__}"
            )

        # Delete via API
        self.vis.api.delete_bookmark(index)

        # Update local cache
        self.vis._bookmarks.pop(index, None)

    def __iter__(self):
        """Iterate over bookmark indices.

        Yields
        ------
        int
            Frame indices that have bookmarks
        """
        # Fetch from server to get current state
        bookmarks = self.vis.api.get_all_bookmarks()
        return iter(bookmarks.keys())

    def __len__(self) -> int:
        """Return the number of bookmarks.

        Returns
        -------
        int
            Number of bookmarks
        """
        # Fetch from server to get current count
        bookmarks = self.vis.api.get_all_bookmarks()
        return len(bookmarks)

    def update(self, other=(), /, **kwds) -> None:
        """Update bookmarks with key-value pairs.

        Parameters
        ----------
        other : dict or iterable
            A dictionary or iterable of (index, label) pairs
        **kwds
            Additional keyword arguments (index=label)

        Examples
        --------
        >>> vis.bookmarks.update({0: "First", 1: "Second"})
        >>> vis.bookmarks.update([(2, "Third"), (3, "Fourth")])
        >>> vis.bookmarks.update({4: "Fifth"}, sixth=5)

        Raises
        ------
        TypeError
            If indices are not integers or labels are not strings
        IndexError
            If any index is out of range
        ValueError
            If any label is empty
        """
        # Handle dict-like object
        if hasattr(other, "items"):
            for index, label in other.items():
                self[index] = label
        # Handle iterable of (index, label) pairs
        elif hasattr(other, "__iter__"):
            for index, label in other:
                self[index] = label

        # Handle keyword arguments
        for index, label in kwds.items():
            # Convert string keys to integers if they look like integers
            try:
                idx = int(index)
            except ValueError:
                raise TypeError(
                    f"Bookmark index must be an integer, not {type(index).__name__}"
                )
            self[idx] = label

    def clear(self) -> None:
        """Remove all bookmarks.

        Examples
        --------
        >>> vis.bookmarks.clear()
        >>> len(vis.bookmarks)
        0
        """
        # Get all current bookmark indices
        indices = list(self)

        # Delete each bookmark
        for index in indices:
            try:
                del self[index]
            except KeyError:
                # Bookmark was already deleted (race condition)
                pass

    def pop(self, index: int, default=None) -> str:
        """Remove and return bookmark label at index.

        Parameters
        ----------
        index : int
            Frame index
        default : str, optional
            Value to return if bookmark doesn't exist

        Returns
        -------
        str
            The bookmark label that was removed

        Raises
        ------
        KeyError
            If bookmark doesn't exist and no default provided

        Examples
        --------
        >>> label = vis.bookmarks.pop(0)
        >>> label = vis.bookmarks.pop(1, "Default")
        """
        try:
            # Get label from cache first (don't fetch from server)
            if index in self.vis._bookmarks:
                label = self.vis._bookmarks[index]
            else:
                # Not in cache, but might exist on server
                response = self.vis.api.get_bookmark(index)
                label = response.get("label")
                if label is None:
                    if default is None:
                        raise KeyError(f"Bookmark at index {index} does not exist")
                    return default

            # Delete the bookmark
            del self[index]
            return label
        except (KeyError, IndexError):
            if default is None:
                raise KeyError(f"Bookmark at index {index} does not exist")
            return default

    def __repr__(self) -> str:
        """Return string representation of bookmarks."""
        return f"Bookmarks({self.vis._bookmarks!r})"

    def __str__(self) -> str:
        """Return string representation of bookmarks."""
        return f"Bookmarks({self.vis._bookmarks!r})"
