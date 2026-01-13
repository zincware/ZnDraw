"""Session management for frontend browser windows.

Provides Python API to access and control frontend sessions (browser windows).
Each frontend session has its own camera viewport state and rendering settings.

Note: Python clients do NOT appear in vis.sessions - only frontend browsers.
Note: Session cameras are stored as geometries with key 'cam:session:{session_id}'.
"""

from collections.abc import Mapping

from zndraw.geometries import Camera
from zndraw.settings import RoomConfig


def get_session_camera_key(session_id: str) -> str:
    """Get the geometry key for a session camera."""
    return f"cam:session:{session_id}"


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
    >>> # Get camera via geometries
    >>> cam = session.camera
    >>> cam.position = (10, 10, 10)
    >>> # Save back (syncs to Redis and frontend via geometry system)
    >>> session.camera = cam
    >>>
    >>> # Get settings (fetches from Redis)
    >>> settings = session.settings
    >>> settings.studio_lighting.key_light = 1.5
    >>> # Auto-saves on attribute change via callback
    """

    def __init__(self, vis: "ZnDraw", session_id: str):
        self._vis = vis
        self.session_id = session_id

    @property
    def camera_key(self) -> str:
        """Get the geometry key for this session's camera."""
        return get_session_camera_key(self.session_id)

    @property
    def camera(self) -> Camera:
        """Get the currently active camera.

        Returns the Camera that the session is currently viewing through
        (as set by active_camera). Modify and set back to sync.

        Returns
        -------
        Camera
            The camera model with position, target, fov, etc.

        Examples
        --------
        >>> cam = session.camera
        >>> cam.position = (10, 10, 10)
        >>> session.camera = cam  # Updates the active camera
        """
        return self._vis.geometries[self.active_camera]

    @camera.setter
    def camera(self, value: Camera) -> None:
        """Set the currently active camera.

        Updates the camera geometry that this session is viewing through.
        Syncs to Redis and broadcasts INVALIDATE_GEOMETRY to notify frontend.

        Parameters
        ----------
        value : Camera
            The camera model to set.
        """
        self._vis.geometries[self.active_camera] = value

    def _make_settings_callback(self, category: str):
        """Create a callback for auto-saving settings on attribute change."""

        def callback(data: dict) -> None:
            self._vis.api.set_session_settings(self.session_id, {category: data})

        return callback

    @property
    def settings(self) -> RoomConfig:
        """Get session-scoped rendering settings.

        Returns a RoomConfig Pydantic model with validation. Modifying attributes
        automatically syncs to the backend via callbacks.

        Returns
        -------
        RoomConfig
            Settings model with studio_lighting, pathtracing, property_inspector.

        Examples
        --------
        >>> settings = session.settings
        >>> settings.studio_lighting.key_light = 1.5  # Auto-saves
        >>> settings.pathtracing.enabled = True  # Auto-saves
        """
        data = self._vis.api.get_session_settings(self.session_id)
        config = RoomConfig.model_validate(data)

        # Attach callbacks to each sub-setting for auto-save
        config.studio_lighting.callback = self._make_settings_callback(
            "studio_lighting"
        )
        config.pathtracing.callback = self._make_settings_callback("pathtracing")
        config.property_inspector.callback = self._make_settings_callback(
            "property_inspector"
        )

        return config

    @settings.setter
    def settings(self, value: RoomConfig) -> None:
        """Set session-scoped rendering settings.

        Parameters
        ----------
        value : RoomConfig
            Settings model to set (validated via Pydantic).
        """
        self._vis.api.set_session_settings(self.session_id, value.model_dump())

    @property
    def active_camera(self) -> str:
        """Get the camera key this session is viewing through.

        Returns
        -------
        str
            The geometry key of the camera being viewed through.
        """
        return self._vis.api.get_active_camera(self.session_id)

    @active_camera.setter
    def active_camera(self, value: str) -> None:
        """Set the camera key this session should view through.

        Parameters
        ----------
        value : str
            The geometry key of a camera to view through.

        Raises
        ------
        KeyError
            If the camera key doesn't exist in geometries.
        TypeError
            If the geometry exists but is not a Camera.
        """
        # Validate camera exists
        if value not in self._vis.geometries:
            raise KeyError(f"Camera '{value}' not found in geometries")

        # Validate it's actually a Camera type
        geom = self._vis.geometries[value]
        if not isinstance(geom, Camera):
            raise TypeError(
                f"Geometry '{value}' is not a Camera (got {type(geom).__name__})"
            )

        self._vis.api.set_active_camera(self.session_id, value)

    def __repr__(self) -> str:
        return f"FrontendSession({self.session_id!r})"


class FrontendSessions(Mapping):
    """Collection of frontend sessions, keyed by session_id.

    Implements the Mapping interface for dict-like access.
    Only frontend (browser) sessions appear here, not Python clients.

    Supports:
    - Dict access: vis.sessions["abc-123"]
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

    def get(self, session_id: str | None = None) -> FrontendSession | None:
        """Get session by ID.

        Parameters
        ----------
        session_id : str, optional
            Session ID lookup.

        Returns
        -------
        FrontendSession or None
            The session, or None if not found.

        Examples
        --------
        >>> session = vis.sessions.get("abc-123")
        """
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
