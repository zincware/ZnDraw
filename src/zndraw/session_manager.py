"""Session management for frontend browser windows.

Provides Python API to access and control frontend sessions (browser windows).
Each frontend session has its own camera viewport state and rendering settings.

Note: Python clients do NOT appear in vis.sessions - only frontend browsers.
"""

from collections.abc import Mapping

from zndraw.geometries import Camera


class FrontendSession:
    """Represents a single browser window/tab.

    Camera and settings are accessed as properties that fetch from / sync to Redis.
    Each frontend session has independent camera and rendering settings.

    Parameters
    ----------
    vis : ZnDraw
        The parent ZnDraw instance.
    session_id : str
        The unique session identifier.

    Examples
    --------
    >>> session = vis.sessions["abc-123"]
    >>> # Get camera (fetches from Redis)
    >>> cam = session.camera
    >>> cam.position = (10, 10, 10)
    >>> # Save back (syncs to Redis and frontend)
    >>> session.camera = cam
    """

    def __init__(self, vis: "ZnDraw", session_id: str):
        self._vis = vis
        self.session_id = session_id

    @property
    def alias(self) -> str | None:
        """User-defined alias for this session.

        Aliases provide stable access to sessions across reconnects.
        Set from frontend via URL parameter (?alias=projector) or UI.

        Returns
        -------
        str or None
            The alias, or None if not set.
        """
        return self._vis.api.get_session_alias(self.session_id)

    @alias.setter
    def alias(self, value: str | None) -> None:
        """Set alias for stable access.

        Parameters
        ----------
        value : str or None
            The alias to set, or None to remove.
        """
        self._vis.api.set_session_alias(self.session_id, value)

    @property
    def camera(self) -> Camera:
        """Get camera state (fetches from Redis).

        Returns the Camera Pydantic model directly. Modify and set back to sync.

        Returns
        -------
        Camera
            The camera model with position, target, fov, etc.

        Examples
        --------
        >>> cam = session.camera
        >>> cam.position = (10, 10, 10)
        >>> session.camera = cam  # Sync to frontend
        """
        return self._vis.api.get_session_camera(self.session_id)

    @camera.setter
    def camera(self, value: Camera) -> None:
        """Set camera state (syncs to Redis and notifies frontend).

        Parameters
        ----------
        value : Camera
            The camera model to set.
        """
        self._vis.api.set_session_camera(self.session_id, value)

    @property
    def settings(self) -> dict:
        """Get session-scoped rendering settings.

        Returns
        -------
        dict
            Settings dictionary with camera_type, lighting, etc.
        """
        return self._vis.api.get_session_settings(self.session_id)

    @settings.setter
    def settings(self, value: dict) -> None:
        """Set session-scoped rendering settings.

        Parameters
        ----------
        value : dict
            Settings dictionary to set.
        """
        self._vis.api.set_session_settings(self.session_id, value)

    def __repr__(self) -> str:
        alias_str = f" (alias={self.alias!r})" if self.alias else ""
        return f"FrontendSession({self.session_id!r}{alias_str})"


class FrontendSessions(Mapping):
    """Collection of frontend sessions, keyed by session_id.

    Implements the Mapping interface for dict-like access.
    Only frontend (browser) sessions appear here, not Python clients.

    Supports:
    - Dict access: vis.sessions["abc-123"]
    - Alias lookup: vis.sessions.get(alias="projector")
    - Iteration: for session in vis.sessions.values()
    - Length: len(vis.sessions)

    Parameters
    ----------
    vis : ZnDraw
        The parent ZnDraw instance.

    Examples
    --------
    >>> # List all sessions
    >>> for sid, session in vis.sessions.items():
    ...     print(f"{sid}: {session.camera.position}")

    >>> # Access by alias
    >>> projector = vis.sessions.get(alias="projector")
    >>> if projector:
    ...     projector.camera.position = (10, 10, 10)
    """

    def __init__(self, vis: "ZnDraw"):
        self._vis = vis

    def _get_frontend_session_ids(self) -> list[str]:
        """Get session IDs that have frontend cameras (excludes Python clients).

        Returns
        -------
        list[str]
            List of frontend session IDs.
        """
        return self._vis.api.list_frontend_sessions()

    def __getitem__(self, session_id: str) -> FrontendSession:
        """Access session by session_id.

        Parameters
        ----------
        session_id : str
            The session identifier.

        Returns
        -------
        FrontendSession
            The session object.

        Raises
        ------
        KeyError
            If session_id not found.
        """
        session_ids = self._get_frontend_session_ids()
        if session_id not in session_ids:
            raise KeyError(f"Session '{session_id}' not found")
        return FrontendSession(self._vis, session_id)

    def __len__(self) -> int:
        """Return number of frontend sessions."""
        return len(self._get_frontend_session_ids())

    def __iter__(self):
        """Iterate over session IDs."""
        return iter(self._get_frontend_session_ids())

    def get(
        self, session_id: str | None = None, *, alias: str | None = None
    ) -> FrontendSession | None:
        """Get session by ID or alias.

        Parameters
        ----------
        session_id : str, optional
            Direct session ID lookup.
        alias : str, optional
            Lookup by user-defined alias.

        Returns
        -------
        FrontendSession or None
            The session, or None if not found.

        Examples
        --------
        >>> # Lookup by session ID
        >>> session = vis.sessions.get("abc-123")

        >>> # Lookup by alias
        >>> projector = vis.sessions.get(alias="projector")
        """
        if alias is not None:
            session_id = self._vis.api.get_session_by_alias(alias)
            if session_id is None:
                return None
        if session_id is None:
            return None
        try:
            return self[session_id]
        except KeyError:
            return None

    def values(self):
        """Iterate over FrontendSession objects.

        Yields
        ------
        FrontendSession
            Each frontend session.
        """
        for session_id in self:
            yield FrontendSession(self._vis, session_id)

    def items(self):
        """Iterate over (session_id, FrontendSession) pairs.

        Yields
        ------
        tuple[str, FrontendSession]
            Session ID and session object pairs.
        """
        for session_id in self:
            yield session_id, FrontendSession(self._vis, session_id)

    def __repr__(self) -> str:
        session_ids = self._get_frontend_session_ids()
        return f"FrontendSessions({session_ids})"
