import contextlib
import dataclasses
import functools
import logging
import time
import typing as t
import uuid
from collections.abc import MutableSequence

import ase
import msgpack
import numpy as np
import requests
from pydantic import BaseModel, Field

from zndraw.api_manager import APIManager
from zndraw.bookmarks_manager import Bookmarks
from zndraw.exceptions import LockError
from zndraw.extensions import Extension, Category
from zndraw.figures_manager import Figures
from zndraw.frame_cache import FrameCache
from zndraw.metadata_manager import RoomMetadata
from zndraw.scene_manager import Geometries
from zndraw.server_manager import get_server_status
from zndraw.settings import RoomConfig, settings
from zndraw.socket_manager import SocketManager
from zndraw.lock import ZnDrawLock
from asebytes import encode, decode
from zndraw.connectivity import add_connectivity
from zndraw.utils import update_colors_and_radii
from zndraw.version_utils import validate_server_version

from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
)

log = logging.getLogger(__name__)


def _is_celery_extension(extension_name: str, category: str) -> bool:
    """Check if an extension is a server-side Celery extension.

    Parameters
    ----------
    extension_name : str
        The extension class name
    category : str
        The extension category

    Returns
    -------
    bool
        True if the extension is a Celery (server-side) extension
    """
    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections

    category_map = {
        "modifiers": modifiers,
        "selections": selections,
        "analysis": analysis,
    }

    return category in category_map and extension_name in category_map[category]


class _GeometryStore(t.TypedDict):
    type: str
    data: dict


class _ExtensionStore(t.TypedDict):
    public: bool
    run_kwargs: dict | None
    extension: t.Type[Extension]


class Selections:
    """Accessor for per-geometry selections."""

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def __getitem__(self, geometry: str) -> tuple[int]:
        response = self.vis.api.get_selection(geometry)
        return tuple(response.get("selection", []))

    def __setitem__(self, geometry: str, indices: t.Iterable[int]) -> None:
        # Validate indices
        if not hasattr(indices, "__iter__") or isinstance(indices, (str, bytes)):
            raise ValueError("Selection must be an iterable of integers.")

        indices_list = list(indices)

        # Check all elements are integers
        if not all(isinstance(idx, int) for idx in indices_list):
            raise ValueError("Selection must be an iterable of integers.")

        # For particles geometry, validate indices are in bounds
        if geometry == "particles" and len(self.vis) > 0:
            num_atoms = len(self.vis.atoms)
            for idx in indices_list:
                if idx < 0 or idx >= num_atoms:
                    raise ValueError(
                        f"Selection index {idx} out of range [0, {num_atoms})."
                    )

        self.vis.api.update_selection(geometry, indices_list)

    def __delitem__(self, geometry: str) -> None:
        self.vis.api.update_selection(geometry, [])

    def __iter__(self):
        data = self.vis.api.get_all_selections()
        return iter(data["selections"].keys())

    def __len__(self) -> int:
        data = self.vis.api.get_all_selections()
        return len(data["selections"])


class SelectionGroups:
    """Accessor for named selection groups."""

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def __getitem__(self, group_name: str) -> dict[str, list[int]]:
        response = self.vis.api.get_selection_group(group_name)
        return response["group"]

    def __setitem__(
        self, group_name: str, group_data: dict[str, t.Iterable[int]]
    ) -> None:
        # Convert iterables to lists
        data = {k: list(v) for k, v in group_data.items()}
        self.vis.api.create_update_selection_group(group_name, data)

    def __delitem__(self, group_name: str) -> None:
        self.vis.api.delete_selection_group(group_name)

    def __iter__(self):
        data = self.vis.api.get_all_selections()
        return iter(data["groups"].keys())

    def __len__(self) -> int:
        data = self.vis.api.get_all_selections()
        return len(data["groups"])


class Screenshots:
    """Accessor for room screenshots."""

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def list(self, limit: int = 20, offset: int = 0) -> list[dict]:
        """List all screenshots for this room.

        Parameters
        ----------
        limit : int
            Maximum number of screenshots to return (default: 20, max: 100).
        offset : int
            Number of screenshots to skip (default: 0).

        Returns
        -------
        list[dict]
            List of screenshot metadata dictionaries.
        """
        result = self.vis.api.list_screenshots(limit, offset)
        return result.get("screenshots", [])

    def get(self, screenshot_id: int) -> bytes:
        """Download screenshot data by ID.

        Parameters
        ----------
        screenshot_id : int
            Screenshot identifier.

        Returns
        -------
        bytes
            The image file data.
        """
        return self.vis.api.download_screenshot(screenshot_id)

    def metadata(self, screenshot_id: int) -> dict:
        """Get metadata for a screenshot.

        Parameters
        ----------
        screenshot_id : int
            Screenshot identifier.

        Returns
        -------
        dict
            Screenshot metadata (id, format, size, dimensions).
        """
        return self.vis.api.get_screenshot_metadata(screenshot_id)

    def delete(self, screenshot_id: int) -> bool:
        """Delete a screenshot.

        Parameters
        ----------
        screenshot_id : int
            Screenshot identifier.

        Returns
        -------
        bool
            True if deleted successfully.
        """
        return self.vis.api.delete_screenshot(screenshot_id)

    @property
    def latest(self) -> dict | None:
        """Get most recent screenshot metadata.

        Returns
        -------
        dict | None
            Latest screenshot metadata, or None if no screenshots exist.
        """
        results = self.list(limit=1)
        return results[0] if results else None

    def __len__(self) -> int:
        """Get total number of screenshots.

        Returns
        -------
        int
            Total screenshot count.
        """
        result = self.vis.api.list_screenshots(limit=1)
        return result.get("total", 0)


class LocalSettings(BaseModel):
    """Local client settings for ZnDraw operations.

    These settings control client-side behavior for operations like uploading frames.
    """

    target_chunk_size_bytes: int = Field(
        default=500_000,
        gt=0,
        le=100_000_000,
        description="Target size in bytes for each upload chunk. Default is 500KB. Maximum is 100MB.",
    )

    show_progress: bool = Field(
        default=True,
        description="Show progress bar during chunked uploads",
    )

    max_retries: int = Field(
        default=3,
        ge=0,
        description="Maximum number of retries for failed chunk uploads",
    )

    retry_delay: float = Field(
        default=1.0,
        ge=0,
        description="Delay in seconds between retries",
    )


class ProgressTracker:
    """Context manager for tracking progress of long-running operations.

    Parameters
    ----------
    vis : ZnDraw
        The ZnDraw instance
    description : str
        Initial progress description

    Examples
    --------
    >>> with vis.progress_tracker("Loading data...") as tracker:
    ...     # Do some work
    ...     tracker.update(progress=50)
    ...     # Do more work
    ...     tracker.update(description="Processing...", progress=100)
    """

    def __init__(self, vis: "ZnDraw", description: str):
        self.vis = vis
        self.description = description
        self.progress_id = str(uuid.uuid4())

    def __enter__(self):
        """Start tracking progress."""
        if not self.vis.socket.sio.connected:
            raise RuntimeError("Client is not connected.")

        # Use REST API to start progress tracking (more reliable than socket emit)
        self.vis.api.progress_start(self.progress_id, self.description)
        return self

    def update(self, description: str | None = None, progress: float | None = None):
        """Update the progress description and/or percentage.

        Parameters
        ----------
        description : str | None
            New progress description
        progress : float | None
            Progress percentage (0-100)
        """
        if not self.vis.socket.sio.connected:
            raise RuntimeError("Client is not connected.")

        # Use REST API to update progress (more reliable than socket emit)
        self.vis.api.progress_update(self.progress_id, description, progress)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Complete progress tracking."""
        if not self.vis.socket.sio.connected:
            # If disconnected, just return without error
            return False

        # Use REST API to complete progress tracking (more reliable than socket emit)
        try:
            self.vis.api.progress_complete(self.progress_id)
        except Exception as e:
            # Log but don't raise - progress completion is best-effort
            import logging
            logging.getLogger(__name__).warning(f"Failed to complete progress {self.progress_id}: {e}")
        return False


@dataclasses.dataclass
class ZnDraw(MutableSequence):
    """A client for interacting with the ZnDraw server.

    Parameters
    ----------
    url
        URL of the ZnDraw server. If None, will attempt to auto-discover
        a running local server (similar to zndraw CLI behavior).
    room
        Name of the room to connect to.
    user
        Username for authentication. If None, server will assign a guest username.
        If provided, will attempt to login/register with that username.
        After login, this field is updated with the actual username from the server.
    password
        Optional password for authentication. Required for registered users.
        For admin accounts, use admin credentials configured on the server.
    auto_pickup_jobs
        Whether to automatically pick up extension jobs.
    description
        Optional description for the room.
    copy_from
        Optional room name to copy initial state from.
    connectivity_threshold
        Maximum number of atoms for automatic connectivity calculation.
        When atoms are added via append/extend/insert/setitem, connectivity
        will be automatically computed if the structure has fewer atoms than
        this threshold and connectivity is not already present. Default: 1000.

    """

    url: str | None = None
    room: str = "default"
    user: str | None = None
    password: str | None = None
    auto_pickup_jobs: bool = True
    description: str | None = None
    copy_from: str | None = None
    connectivity_threshold: int = 1000

    _step: int = 0
    _len: int = 0
    _settings: dict = dataclasses.field(default_factory=dict, init=False)
    _public_extensions: dict[str, _ExtensionStore] = dataclasses.field(
        default_factory=dict, init=False
    )
    _private_extensions: dict[str, _ExtensionStore] = dataclasses.field(
        default_factory=dict, init=False
    )
    _filesystems: dict[str, dict] = dataclasses.field(
        default_factory=dict, init=False
    )
    role: str = dataclasses.field(default="guest", init=False)
    _worker_id: str | None = dataclasses.field(
        default=None, init=False
    )  # Server's sid for worker identification
    _selection: frozenset[int] = frozenset()
    _frame_selection: frozenset[int] = frozenset()
    _bookmarks: dict[int, str] = dataclasses.field(default_factory=dict, init=False)
    _geometries: dict[str, _GeometryStore] = dataclasses.field(
        default_factory=dict, init=False
    )
    _figures: dict[str, dict] = dataclasses.field(default_factory=dict, init=False)
    _metadata: RoomMetadata | None = dataclasses.field(default=None, init=False)

    def __post_init__(self):
        # Auto-discover local server if url is None
        if self.url is None:
            is_running, server_info, status_message = get_server_status()

            if is_running and server_info is not None:
                self.url = f"http://127.0.0.1:{server_info.port}"
                log.info(f"Auto-discovered local ZnDraw server: {status_message}")
            else:
                raise RuntimeError(
                    "No local ZnDraw server found. Please start a server with 'zndraw' "
                    "or provide an explicit URL."
                )

        self.url = self.url.rstrip("/")

        # Create APIManager (user_name will be set after login)
        self.api = APIManager(url=self.url, room=self.room)
        self.cache: FrameCache | None = FrameCache(maxsize=100)

        # Initialize local settings
        self.local = LocalSettings()

        # Validate server version compatibility before connecting
        import zndraw

        validate_server_version(self.api, zndraw.__version__)

        # Step 1: Login to get JWT token and userName
        # If user is None, server will assign a guest username
        # If user is provided, we try to login/register with that name
        login_data = self.api.login(user_name=self.user, password=self.password)

        # Update self.user with the actual username from server
        # (may be different if user was None and server assigned a guest name)
        self.user = login_data["userName"]
        self.role = login_data.get("role", "guest")
        log.info(f"Logged in as {self.user} (role: {self.role})")

        # Step 2: Join room (authenticated with JWT)
        response_data = self.api.join_room(
            description=self.description,
            copy_from=self.copy_from,
        )

        # Create socket manager and connect (with JWT)
        self.socket = SocketManager(zndraw_instance=self)
        self.connect()

        # Fetch room data separately via REST endpoints
        # (join response is now minimal - only returns status, sessionId, userName, roomId, created)

        # Get frame count
        try:
            room_info = self.api.get_room_info()
            self._len = room_info["frameCount"]
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                log.debug(f"Room {self.room} not found or has no frames")
                self._len = 0
            else:
                log.error(f"Failed to fetch room info: {e}")
                self._len = 0
        except Exception as e:
            log.error(f"Unexpected error loading room info: {e}")
            self._len = 0

        # Get frame selection
        try:
            frame_selection = self.api.get_frame_selection()
            if frame_selection is not None:
                self._frame_selection = frozenset(frame_selection)
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                log.debug("No frame selection for room")
            else:
                log.warning(f"Failed to fetch frame selection: {e}")
        except Exception as e:
            log.error(f"Unexpected error loading frame selection: {e}")

        # Get bookmarks
        try:
            bookmarks = self.api.get_all_bookmarks()
            if bookmarks:
                self._bookmarks = bookmarks
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                log.debug("No bookmarks for room")
            else:
                log.warning(f"Failed to fetch bookmarks: {e}")
        except Exception as e:
            log.error(f"Unexpected error loading bookmarks: {e}")

        # Get current step
        try:
            step_data = self.api.get_step()
            if step_data.get("step") is not None:
                self._step = int(step_data["step"])
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                log.debug("No current step set for room")
            else:
                log.warning(f"Failed to fetch current step: {e}")
        except Exception as e:
            log.error(f"Unexpected error loading current step: {e}")

        # Get geometries (returns dict of {key: {type, data}})
        try:
            geometries_response = self.api.list_geometries()
            # list_geometries now returns full geometry data as dict
            if isinstance(geometries_response, dict):
                self._geometries = geometries_response
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                log.debug("No geometries for room")
            else:
                log.warning(f"Failed to fetch geometries: {e}")
        except Exception as e:
            log.error(f"Unexpected error loading geometries: {e}")

    @classmethod
    def for_job_execution(
        cls,
        url: str,
        room: str,
        user: str | None = None,
        password: str | None = None,
    ) -> "ZnDraw":
        """Create a ZnDraw instance for job execution in a specific room.

        This factory method creates a fresh ZnDraw instance connected to the
        specified room with proper authentication. Used when an extension worker
        needs to execute a job in a different room than the one it registered in.

        Parameters
        ----------
        url : str
            URL of the ZnDraw server
        room : str
            Name of the room to connect to (the room where the job was triggered)
        user : str | None
            Username for authentication. If None, server will assign a guest username.
        password : str | None
            Optional password for authentication

        Returns
        -------
        ZnDraw
            A new ZnDraw instance connected to the specified room with auto_pickup_jobs=False
        """
        return cls(
            url=url,
            room=room,
            user=user,
            password=password,
            auto_pickup_jobs=False,
        )

    def get_lock(self, msg: str | None = None, target: str = "trajectory:meta") -> ZnDrawLock:
        """Get a SocketIOLock instance for distributed locking.

        Parameters
        ----------
        msg : str | None
            Optional message describing the lock purpose
        target : str
            Lock target identifier (default: "trajectory:meta")

        Returns
        -------
        SocketIOLock
            Lock instance ready to use as context manager

        Raises
        ------
        RuntimeError
            If the client is not connected

        Examples
        --------
        >>> with vis.get_lock(msg="Uploading trajectory") as lock:
        ...     vis.extend(frames)
        ...     lock.update_msg("Upload 50% complete")
        ...     vis.extend(more_frames)
        """
        if not self.socket.sio.connected:
            raise RuntimeError("Lock requires an active connection. Ensure client is connected.")
        return ZnDrawLock(self.api, target=target, msg=msg)

    @property
    def geometries(self) -> Geometries:
        return Geometries(self)

    @property
    def figures(self) -> Figures:
        return Figures(self)

    @property
    def metadata(self) -> RoomMetadata:
        """Access room metadata as a dict-like object.

        Returns
        -------
        RoomMetadata
            A MutableMapping interface to room metadata.

        Examples
        --------
        >>> vis.metadata["file"] = "data.xyz"
        >>> print(vis.metadata["file"])
        >>> del vis.metadata["file"]
        """
        if self._metadata is None:
            self._metadata = RoomMetadata(self)
        return self._metadata

    @property
    def selections(self) -> Selections:
        """Access selections by geometry name.

        Returns
        -------
        Selections
            A dict-like interface to per-geometry selections.

        Examples
        --------
        >>> vis.selections["particles"] = [1, 2, 3]
        >>> print(vis.selections["particles"])
        >>> del vis.selections["particles"]
        """
        if not hasattr(self, "_selections_accessor"):
            self._selections_accessor = Selections(self)
        return self._selections_accessor

    @property
    def selection_groups(self) -> SelectionGroups:
        """Access named selection groups.

        Returns
        -------
        SelectionGroups
            A dict-like interface to selection groups.

        Examples
        --------
        >>> vis.selection_groups["group1"] = {"particles": [1, 2], "forces": [3]}
        >>> print(vis.selection_groups["group1"])
        >>> del vis.selection_groups["group1"]
        """
        if not hasattr(self, "_selection_groups_accessor"):
            self._selection_groups_accessor = SelectionGroups(self)
        return self._selection_groups_accessor

    @property
    def active_selection_group(self) -> str | None:
        """Get the currently active selection group name.

        Returns
        -------
        str | None
            The name of the active selection group, or None if no group is active.
        """
        data = self.api.get_all_selections()
        return data.get("activeGroup")

    def load_selection_group(self, group_name: str) -> None:
        """Load a selection group (apply it to current selections).

        Parameters
        ----------
        group_name : str
            Name of the group to load.
        """
        self.api.load_selection_group(group_name)

    @property
    def sid(self) -> str | None:
        """Return the worker ID assigned by the server.

        The server assigns a worker ID (its request.sid) during extension registration.
        This ID is used consistently for both registration and disconnect cleanup.

        Returns
        -------
        str | None
            The worker ID assigned by server, client's socket.sio.sid if not yet registered,
            or None if not connected.
        """
        return self._worker_id if self._worker_id else self.socket.sio.sid

    @property
    def is_admin(self) -> bool:
        """Check if current user has admin privileges.

        Returns True if the user's role is "admin".

        In local mode (no admin credentials configured on server),
        all authenticated users have admin role. In deployment mode
        (admin credentials configured), only users who logged in
        with correct admin credentials have admin role.

        Returns
        -------
        bool
            True if user has admin role, False otherwise
        """
        return self.role == "admin"

    @property
    def step(self) -> int:
        return self._step

    @step.setter
    def step(self, value: int):
        if not isinstance(value, int) or value < 0:
            raise ValueError("Step must be a non-negative integer.")
        if value >= self._len:
            raise ValueError(
                f"Step {value} is out of bounds. Current number of frames: {self._len}."
            )
        self._step = value

        # Use REST API with auto-lock pattern: acquire lock → update step → release lock
        with self.get_lock(target="step"):
            self.api.update_step(value)

    @property
    def points(self) -> np.ndarray:
        """Get the current frame as an `ase.Atoms` object."""
        from zndraw.geometries import Curve

        curve: Curve | None = self.geometries.get("curve")
        if curve is not None:
            if isinstance(curve.position, str):
                raise ValueError(
                    "Curve position is string; cannot retrieve static points."
                )
            return np.array(curve.position)
        return np.empty((0, 3))

    @property
    def atoms(self) -> ase.Atoms:
        """Get the current frame as an `ase.Atoms` object."""
        return self[self.step]

    @atoms.setter
    def atoms(self, value: ase.Atoms):
        """Set the current frame from an `ase.Atoms` object."""
        self[self.step] = value

    @property
    def selection(self) -> tuple[int]:
        """Get selection for 'particles' geometry.

        Returns
        -------
        frozenset[int]
            The current selection indices for particles.
        """
        return self.selections["particles"]

    @selection.setter
    def selection(self, value: t.Iterable[int] | None):
        """Set selection for 'particles' geometry.

        Parameters
        ----------
        value : t.Iterable[int] | None
            The selection indices to set, or None to clear the selection.
        """
        if value is None:
            value = []
        self.selections["particles"] = value

    @property
    def frame_selection(self) -> tuple[int, ...]:
        return tuple(sorted(self._frame_selection))

    @frame_selection.setter
    def frame_selection(self, value: t.Iterable[int] | None):
        indices = [] if value is None else list(value)
        if not all(isinstance(idx, int) and 0 <= idx < len(self) for idx in indices):
            raise ValueError("Selection must be an iterable of valid frame indices.")

        if self.socket.sio.connected:
            response = self.socket.sio.call(
                "frame_selection:set", {"indices": indices}, timeout=5
            )
            if response and not response.get("success", False):
                raise RuntimeError(
                    response.get("message", "Failed to set frame selection")
                )
            self._frame_selection = frozenset(indices)
        else:
            raise RuntimeError("Client is not connected.")

    @property
    def bookmarks(self) -> Bookmarks:
        """Access bookmarks using MutableMapping interface.

        Returns
        -------
        Bookmarks
            A dict-like interface to frame bookmarks.

        Examples
        --------
        >>> vis.bookmarks[0] = "First Frame"
        >>> print(vis.bookmarks[0])
        >>> del vis.bookmarks[0]
        >>> len(vis.bookmarks)
        >>> list(vis.bookmarks)
        """
        if not hasattr(self, "_bookmarks_accessor"):
            self._bookmarks_accessor = Bookmarks(self)
        return self._bookmarks_accessor

    @property
    def screenshots(self) -> Screenshots:
        """Access screenshots for this room.

        Returns
        -------
        Screenshots
            An interface to manage room screenshots.

        Examples
        --------
        >>> # List all screenshots
        >>> screenshots = vis.screenshots.list()
        >>>
        >>> # Get latest screenshot
        >>> latest = vis.screenshots.latest
        >>>
        >>> # Download screenshot data
        >>> data = vis.screenshots.get(screenshot_id)
        >>>
        >>> # Delete a screenshot
        >>> vis.screenshots.delete(screenshot_id)
        """
        if not hasattr(self, "_screenshots_accessor"):
            self._screenshots_accessor = Screenshots(self)
        return self._screenshots_accessor

    def connect(self):
        self.socket.connect()

    def disconnect(self):
        self.socket.disconnect()

    def _upload_frames(self, action: str, data, **kwargs):
        """Internal upload method - does NOT acquire lock.

        Server validates locks, so only public methods need to acquire.
        """
        if not self.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        return self.api.upload_frames(action, data, **kwargs)

    def _append_frame(self, data: dict):
        """Internal append - does NOT acquire lock."""
        self._upload_frames("append", data)

    def _extend_frames(self, data: list[dict]):
        """Internal extend - does NOT acquire lock."""
        result = self._upload_frames("extend", data)
        return result.get("new_indices", [])

    def _replace_frame(self, frame_id: int, data: dict):
        """Internal replace - does NOT acquire lock."""
        self._upload_frames("replace", data, frame_id=frame_id)

    def _insert_frame(self, index: int, data: dict):
        """Internal insert - does NOT acquire lock."""
        self._upload_frames("insert", data, insert_position=index)

    def __len__(self) -> int:
        return self._len

    @t.overload
    def get(self, index: int, keys: list[str] | None = None) -> dict[str, t.Any]: ...
    @t.overload
    def get(self, index: list[int], keys: list[str] | None = None) -> list[dict[str, t.Any]]: ...
    @t.overload
    def get(self, index: slice, keys: list[str] | None = None) -> list[dict[str, t.Any]]: ...
    def _decode_frame_dict(self, frame_data: dict[bytes, bytes]) -> dict[str, t.Any]:
        """Decode dict[bytes, bytes] to user-friendly dict[str, Any] format.

        Converts byte keys to strings and unpacks msgpack-encoded byte values.
        """
        import msgpack
        import msgpack_numpy as m

        result = {}
        for k, v in frame_data.items():
            # Decode key from bytes to string
            key = k.decode() if isinstance(k, bytes) else k

            # Unpack value from msgpack bytes
            if isinstance(v, bytes):
                try:
                    value = msgpack.unpackb(v, object_hook=m.decode, strict_map_key=False)
                except:
                    value = v
            else:
                value = v

            result[key] = value

        return result

    def get(self, index, keys: list[str] | None = None):
        """Get frame data as dictionaries with decoded values.

        Internally caches raw dict[bytes, bytes] format for efficiency, but decodes
        to user-friendly dict[str, Any] format before returning.

        For ase.Atoms objects, use __getitem__ (vis[index]) instead.
        """
        if isinstance(index, np.ndarray):
            index = index.tolist() if index.ndim > 0 else int(index.item())

        length = len(self)

        if isinstance(index, int):
            normalized_index = index if index >= 0 else length + index
            if not (0 <= normalized_index < length):
                raise IndexError("Index out of range")

            # Only use cache if keys is None (full frame)
            if (
                keys is None
                and self.cache is not None
                and normalized_index in self.cache
            ):
                # Cache stores dict[bytes, bytes], decode to dict[str, Any]
                raw_data = self.cache.get(normalized_index)
                return self._decode_frame_dict(raw_data)

            # Fetch from API (returns dict[bytes, bytes])
            frame_data = self.api.get_frames([normalized_index], keys=keys)[0]

            # Cache raw bytes format
            if keys is None and self.cache is not None:
                self.cache.set(normalized_index, frame_data)

            # Return decoded format
            return self._decode_frame_dict(frame_data)

        elif isinstance(index, slice):
            # Pass slice directly to API for efficiency
            indices = list(range(*index.indices(length)))
            if not indices:
                return []

            # Fetch frames using slice (returns list[dict[bytes, bytes]])
            fetched_data = self.api.get_frames(index, keys=keys)

            # Cache individual frames if keys is None (cache stores raw bytes)
            if keys is None and self.cache is not None:
                for idx, frame in zip(indices, fetched_data):
                    self.cache.set(idx, frame)

            # Return decoded format
            return [self._decode_frame_dict(frame) for frame in fetched_data]

        elif isinstance(index, list):
            # Validate all indices are integers
            for idx in index:
                if not isinstance(idx, (int, np.integer)):
                    raise TypeError(
                        f"List indices must be integers, not {type(idx).__name__}"
                    )

            if not index:
                return []

            # Normalize negative indices
            normalized_indices = [idx if idx >= 0 else length + idx for idx in index]

            # Only use cache if keys is None (full frames)
            results_dict, misses = {}, []
            for idx in normalized_indices:
                if keys is None and self.cache is not None and idx in self.cache:
                    # Get from cache (raw bytes) and decode
                    raw_data = self.cache.get(idx)
                    results_dict[idx] = self._decode_frame_dict(raw_data)
                elif idx not in misses:
                    misses.append(idx)

            if misses:
                # Fetch from API (returns list[dict[bytes, bytes]])
                fetched_data = self.api.get_frames(misses, keys=keys)
                for idx, frame in zip(misses, fetched_data):
                    # Cache raw bytes format
                    if keys is None and self.cache is not None:
                        self.cache.set(idx, frame)
                    # Store decoded format in results
                    results_dict[idx] = self._decode_frame_dict(frame)

            return [results_dict[idx] for idx in normalized_indices]

        raise TypeError(
            f"Index must be int, slice, or list, not {type(index).__name__}"
        )

    @t.overload
    def __getitem__(self, index: int) -> ase.Atoms: ...
    @t.overload
    def __getitem__(self, index: slice) -> list[ase.Atoms]: ...
    @t.overload
    def __getitem__(self, index: list[int]) -> list[ase.Atoms]: ...
    @t.overload
    def __getitem__(self, index: np.ndarray) -> ase.Atoms | list[ase.Atoms]: ...
    def __getitem__(self, index):
        """Get frame(s) as ase.Atoms object(s).

        Retrieves raw dict[bytes, bytes] from cache/server and decodes
        on-demand to ase.Atoms.
        """
        # We need to get the raw dict[bytes, bytes] format, not the decoded format from get()
        # So we'll fetch directly from cache/API and bypass get()'s decoding
        if isinstance(index, np.ndarray):
            index = index.tolist() if index.ndim > 0 else int(index.item())

        length = len(self)

        if isinstance(index, int):
            normalized_index = index if index >= 0 else length + index
            if not (0 <= normalized_index < length):
                raise IndexError("Index out of range")

            # Check cache first
            if self.cache is not None and normalized_index in self.cache:
                raw_data = self.cache.get(normalized_index)
            else:
                # Fetch from API
                raw_data = self.api.get_frames([normalized_index], keys=None)[0]
                if self.cache is not None:
                    self.cache.set(normalized_index, raw_data)

            return decode(raw_data)

        elif isinstance(index, slice):
            indices = list(range(*index.indices(length)))
            if not indices:
                return []

            # Fetch frames (check cache for each)
            frames = []
            misses = []
            miss_indices = []

            for idx in indices:
                if self.cache is not None and idx in self.cache:
                    frames.append(self.cache.get(idx))
                else:
                    misses.append(idx)
                    miss_indices.append(len(frames))
                    frames.append(None)  # Placeholder

            if misses:
                fetched_data = self.api.get_frames(misses, keys=None)
                for i, (idx, frame) in enumerate(zip(misses, fetched_data)):
                    if self.cache is not None:
                        self.cache.set(idx, frame)
                    frames[miss_indices[i]] = frame

            return [decode(f) for f in frames]

        elif isinstance(index, list):
            for idx in index:
                if not isinstance(idx, (int, np.integer)):
                    raise TypeError(
                        f"List indices must be integers, not {type(idx).__name__}"
                    )

            if not index:
                return []

            normalized_indices = [idx if idx >= 0 else length + idx for idx in index]

            # Fetch frames (check cache for each)
            results_dict, misses = {}, []
            for idx in normalized_indices:
                if self.cache is not None and idx in self.cache:
                    results_dict[idx] = self.cache.get(idx)
                elif idx not in misses:
                    misses.append(idx)

            if misses:
                fetched_data = self.api.get_frames(misses, keys=None)
                for idx, frame in zip(misses, fetched_data):
                    if self.cache is not None:
                        self.cache.set(idx, frame)
                    results_dict[idx] = frame

            return [decode(results_dict[idx]) for idx in normalized_indices]

        raise TypeError(
            f"Index must be int, slice, or list, not {type(index).__name__}"
        )

    @t.overload
    def __setitem__(self, index: int, atoms: ase.Atoms) -> None: ...
    @t.overload
    def __setitem__(self, index: slice, atoms: list[ase.Atoms]) -> None: ...
    @t.overload
    def __setitem__(self, index: list[int], atoms: list[ase.Atoms]) -> None: ...
    @t.overload
    def __setitem__(
        self, index: np.ndarray, atoms: list[ase.Atoms] | ase.Atoms
    ) -> None: ...
    def __setitem__(self, index, atoms):
        if isinstance(atoms, list):
            if not all(isinstance(a, ase.Atoms) for a in atoms):
                raise TypeError("All elements must be ase.Atoms objects")
            dicts = []
            for atom in atoms:
                self._prepare_atoms(atom)
                dicts.append(encode(atom))
            self.set_frames(index, dicts)
        elif isinstance(atoms, ase.Atoms):
            self._prepare_atoms(atoms)
            self.set_frames(index, encode(atoms))
        else:
            raise TypeError("Only ase.Atoms or list of ase.Atoms are supported.")

    def set_frames(self, index, value):
        if not self.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        if isinstance(index, np.ndarray):
            index = index.tolist() if index.ndim > 0 else int(index.item())

        length = len(self)

        if isinstance(index, int):
            if index < 0:
                index += length
            # Single frame replacement - acquire lock
            with self.get_lock(msg=f"Replacing frame at index {index}"):
                self._replace_frame(index, value)
        elif isinstance(index, (slice, list)):
            with self.get_lock(msg=f"Replacing frames in index {index}"):
                if isinstance(index, slice):
                    start, stop, step = index.indices(length)
                    if step == 1:
                        self.api.bulk_patch_frames(value, start=start, stop=stop)
                    else:
                        indices = list(range(start, stop, step))
                        if len(value) != len(indices):
                            raise ValueError(
                                f"attempt to assign sequence of size {len(value)} to extended slice of size {len(indices)}"
                            )
                        self.api.bulk_patch_frames(value, indices=indices)
                else:  # list
                    if len(value) != len(index):
                        raise ValueError("Attempt to assign sequence of wrong size.")
                    # Validate all indices are integers
                    for idx in index:
                        if not isinstance(idx, (int, np.integer)):
                            raise TypeError(
                                f"List indices must be integers, not {type(idx).__name__}"
                            )
                    # Normalize negative indices
                    normalized_indices = [
                        idx if idx >= 0 else length + idx for idx in index
                    ]
                    self.api.bulk_patch_frames(value, indices=normalized_indices)
        else:
            raise TypeError(
                f"Index must be int, slice, or list, not {type(index).__name__}"
            )

    def __delitem__(self, index: int | slice | list[int] | np.ndarray):
        if isinstance(index, np.ndarray):
            index = index.tolist() if index.ndim > 0 else int(index.item())

        # Validate index type
        if not isinstance(index, (int, slice, list)):
            raise TypeError(
                f"Indices must be integers, slices, or lists, not {type(index).__name__}"
            )

        length = len(self)

        # Normalize negative indices and validate types
        if isinstance(index, int):
            index = index if index >= 0 else length + index
        elif isinstance(index, list):
            # Validate all indices are integers
            for idx in index:
                if not isinstance(idx, (int, np.integer)):
                    raise TypeError(
                        f"List indices must be integers, not {type(idx).__name__}"
                    )
            index = [idx if idx >= 0 else length + idx for idx in index]

        with self.get_lock(msg="Deleting frames"):
            self.api.delete_frames(index)

    def _prepare_atoms(self, atoms: ase.Atoms) -> None:
        """Prepare atoms for upload: add connectivity if needed, update colors and radii.

        Parameters
        ----------
        atoms
            The atoms object to prepare.

        """
        if len(atoms) < self.connectivity_threshold and "connectivity" not in atoms.info:
            add_connectivity(atoms)
        update_colors_and_radii(atoms)

    def insert(self, index: int, atoms: ase.Atoms):
        if not isinstance(atoms, ase.Atoms):
            raise TypeError("Only ase.Atoms objects are supported")
        self._prepare_atoms(atoms)
        value = encode(atoms)
        if index < 0:
            index = max(0, len(self) + index + 1)
        elif index > len(self):
            index = len(self)
        # Public API - acquire lock
        with self.get_lock(msg="Inserting frame"):
            self._insert_frame(index, value)

    def append(self, atoms: ase.Atoms):
        if not isinstance(atoms, ase.Atoms):
            raise TypeError("Only ase.Atoms objects are supported")
        self._prepare_atoms(atoms)
        # Public API - acquire lock
        with self.get_lock(msg="Appending frame"):
            self._append_frame(encode(atoms))

    def _calculate_chunk_boundaries(
        self, dicts: list[dict[bytes, bytes]]
    ) -> tuple[list[list[dict[bytes, bytes]]], list[int]]:
        """Calculate chunk boundaries based on exact sizes.

        Parameters
        ----------
        dicts : list[dict[bytes, bytes]]
            List of frame data as msgpack bytes dictionaries

        Returns
        -------
        chunks : list[list[dict[bytes, bytes]]]
            List of chunks (each chunk is a list of msgpack dicts)
        chunk_sizes : list[int]
            Size in bytes of each chunk
        """
        # Get exact sizes by packing the msgpack frames
        packed_frames = [msgpack.packb(frame) for frame in dicts]
        frame_sizes = [len(p) for p in packed_frames]

        # Calculate chunk boundaries
        target_bytes = self.local.target_chunk_size_bytes
        chunks = []
        chunk_sizes = []

        current_chunk = []
        current_size = 0

        for i, (frame_dict, size) in enumerate(zip(dicts, frame_sizes)):
            # Check if adding this frame exceeds target
            if current_size > 0 and current_size + size > target_bytes:
                # Finalize current chunk
                chunks.append(current_chunk)
                chunk_sizes.append(current_size)

                # Start new chunk
                current_chunk = [frame_dict]
                current_size = size
            else:
                # Add to current chunk
                current_chunk.append(frame_dict)
                current_size += size

        # Add last chunk
        if current_chunk:
            chunks.append(current_chunk)
            chunk_sizes.append(current_size)

        return chunks, chunk_sizes

    @contextlib.contextmanager
    def _progress_bar(
        self, total_bytes: int, total_frames: int, num_chunks: int
    ):
        """Create progress bar for chunked upload.

        Yields
        ------
        update_fn : callable
            Function to call with (bytes_uploaded, frames_uploaded) to update progress
        """
        if not self.local.show_progress or num_chunks == 1:
            yield lambda bytes_uploaded, frames_uploaded: None
            return

        # Rich progress bar with full details
        # Note: Using decimal MB (1 MB = 1,000,000 bytes) per SI standard,
        # not binary MiB (1 MiB = 1,048,576 bytes)
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("{task.fields[mb_completed]:.2f}/{task.fields[mb_total]:.2f} MB"),
            TextColumn("({task.percentage:>3.0f}%)"),
            TextColumn("|"),
            TextColumn("{task.fields[frames_done]}/{task.fields[frames_total]} frames"),
            TextColumn("|"),
            TimeElapsedColumn(),
            transient=False,
        ) as progress:
            task = progress.add_task(
                "Uploading",
                total=total_bytes,
                frames_done=0,
                frames_total=total_frames,
                mb_completed=0.0,
                mb_total=total_bytes / 1_000_000,  # Decimal MB
            )

            def update_progress(bytes_uploaded, frames_uploaded):
                progress.update(
                    task,
                    completed=bytes_uploaded,
                    frames_done=frames_uploaded,
                    mb_completed=bytes_uploaded / 1_000_000,
                )

            yield update_progress

    def extend(self, atoms_list: t.Iterable[ase.Atoms]):
        """Extend trajectory with automatic chunked uploads.

        Chunking is always performed based on self.local.target_chunk_size_bytes.
        Progress bar is shown if self.local.show_progress is True and multiple
        chunks are created.

        Warning
        -------
        Uploads are performed in chunks. If an upload fails after retries,
        previously uploaded chunks remain on the server (no automatic rollback).

        Parameters
        ----------
        atoms_list : t.Iterable[ase.Atoms]
            Atoms objects to append to trajectory

        Example
        -------
        >>> vis = ZnDraw(url="...")
        >>>
        >>> # Default: 500KB chunks with progress bar
        >>> vis.extend([atoms for _ in range(1000)])
        >>>
        >>> # Customize chunk size
        >>> vis.local.target_chunk_size_bytes = 1_000_000  # 1MB
        >>> vis.extend([atoms for _ in range(1000)])
        >>>
        >>> # Disable progress bar
        >>> vis.local.show_progress = False
        >>> vis.extend([atoms for _ in range(1000)])
        """
        # Convert all atoms to dicts
        dicts = []
        for atoms in atoms_list:
            if not isinstance(atoms, ase.Atoms):
                raise TypeError("Only ase.Atoms objects are supported")
            self._prepare_atoms(atoms)
            dicts.append(encode(atoms))

        if not dicts:
            log.warning("No frames to upload")
            return

        # Calculate chunk boundaries based on exact sizes
        chunks, chunk_sizes = self._calculate_chunk_boundaries(dicts)

        total_bytes = sum(chunk_sizes)
        total_frames = len(dicts)
        num_chunks = len(chunks)

        log.info(
            f"Uploading {total_frames} frames ({total_bytes / 1_000_000:.2f} MB) "
            f"in {num_chunks} chunk(s)"
        )

        # Upload chunks with progress bar - acquire lock once for entire extend
        with self.get_lock(msg="Uploading frames"):
            with self._progress_bar(total_bytes, total_frames, num_chunks) as update_progress:
                bytes_uploaded = 0
                frames_uploaded = 0

                for chunk_idx, (chunk, chunk_size) in enumerate(zip(chunks, chunk_sizes)):
                    # Upload this chunk with retry logic
                    for attempt in range(self.local.max_retries + 1):
                        try:
                            self._extend_frames(chunk)
                            break  # Success
                        except (ConnectionError, IOError, TimeoutError, OSError) as e:
                            # Only retry network-related exceptions
                            if attempt == self.local.max_retries:
                                # Final attempt failed
                                log.error(
                                    f"Failed to upload chunk {chunk_idx + 1}/{num_chunks} "
                                    f"after {self.local.max_retries} retries: {e}"
                                )
                                raise
                            else:
                                # Retry with backoff
                                delay = self.local.retry_delay * (2**attempt)
                                log.warning(
                                    f"Chunk {chunk_idx + 1}/{num_chunks} failed (attempt {attempt + 1}), "
                                    f"retrying in {delay:.1f}s: {e}"
                                )
                                time.sleep(delay)

                    # Update progress
                    bytes_uploaded += chunk_size
                    frames_uploaded += len(chunk)
                    update_progress(bytes_uploaded, frames_uploaded)

                    log.debug(
                        f"Uploaded chunk {chunk_idx + 1}/{num_chunks}: "
                        f"{len(chunk)} frames, {chunk_size / 1024:.1f} KB"
                    )

        log.info(f"Successfully uploaded {total_frames} frames")

    @property
    def settings(self) -> RoomConfig:
        def callback_fn(data, extension: str):
            self.api.submit_extension_settings(extension, data)

        for key in settings:
            if key not in self._settings:
                data = self.api.get_extension_settings(key)
                self._settings[key] = settings[key](**(data.get("data") or {}))
                self._settings[key].callback = functools.partial(
                    callback_fn, extension=key
                )

        return RoomConfig(**self._settings)

    def register_extension(
        self,
        extension: t.Type[Extension],
        public: bool = False,
        run_kwargs: dict | None = None,
    ):
        name = extension.__name__

        # Select the appropriate dictionary based on public flag
        extensions_dict = self._public_extensions if public else self._private_extensions

        # Check if this specific extension is already registered in this namespace
        if name in extensions_dict:
            namespace = "public" if public else "private"
            raise ValueError(
                f"Extension '{name}' is already registered in {namespace} namespace."
            )

        if not hasattr(extension, "category") or extension.category not in [
            cat.value for cat in Category
        ]:
            raise ValueError("Extension must have a valid 'category' attribute.")

        # Validate that public extensions require admin privileges
        if public and not self.is_admin:
            raise PermissionError(
                "Only admin users can register public extensions. "
                "Please authenticate with admin credentials to register global extensions."
            )

        extensions_dict[name] = {
            "public": public,
            "run_kwargs": run_kwargs,
            "extension": extension,
        }
        print(f"Registered extension '{name}' of category '{extension.category}'.")

        scope = "global" if public else self.room
        print(f"Registering {'global' if public else 'room-scoped'} extension '{name}'...")

        worker_id = self.api.register_extension(
            name=name,
            category=extension.category,
            schema=extension.model_json_schema(),
            socket_manager=self.socket,
            public=public,
        )
        # Store the worker_id assigned by server (server's request.sid)
        if worker_id:
            self._worker_id = worker_id
        print(
            f"Extension '{name}' registered with {scope} (worker_id: {self._worker_id})."
        )

    def run(self, extension: Extension, public: bool | None = None):
        """Run an extension by submitting a job to the server.

        Parameters
        ----------
        extension : Extension
            The extension instance to run
        public : bool | None
            Whether to use the public/global namespace (True) or room-scoped namespace (False).
            If None (default), tries room-scoped first, then auto-retries with public if the
            extension is a server-side extension (Celery-based).

        Returns
        -------
        Job
            Job object for tracking progress and retrieving results
        """
        from zndraw.job import Job

        extension_name = extension.__class__.__name__
        category = extension.category.value

        # If public is None, we'll try private first and let APIManager auto-retry
        # If public is explicitly True/False, we pass that through
        response = self.api.run_extension(
            category=category,
            name=extension_name,
            data=extension.model_dump(),
            public=public if public is not None else False,
            auto_retry=public is None,  # Only auto-retry if public wasn't explicitly set
        )

        # Extract jobId from response and create Job object
        job_id = response.get("jobId")
        if not job_id:
            raise RuntimeError(f"Server response missing jobId: {response}")

        assert self.url is not None, "URL must be set after initialization"
        return Job(job_id=job_id, url=self.url, room=self.room, api=self.api, socket=self.socket)

    def register_filesystem(
        self,
        fs,
        name: str,
        public: bool = False,
    ):
        """Register a filesystem for remote file access.

        Parameters
        ----------
        fs : fsspec.AbstractFileSystem
            An fsspec filesystem instance (e.g., LocalFileSystem, S3FileSystem)
        name : str
            Unique name for this filesystem instance
        public : bool
            If True, register as global filesystem accessible to all rooms.
            Requires admin privileges. Default is False (room-scoped).

        Raises
        ------
        ValueError
            If filesystem name is already registered.
        PermissionError
            If trying to register public filesystem without admin privileges.
        RuntimeError
            If registration with server fails.

        Examples
        --------
        >>> from fsspec.implementations.local import LocalFileSystem
        >>> fs = LocalFileSystem()
        >>> vis.register_filesystem(fs, name="local-data")
        >>>
        >>> # S3 example
        >>> import s3fs
        >>> s3 = s3fs.S3FileSystem(anon=False, key='...', secret='...')
        >>> vis.register_filesystem(s3, name="s3-data", public=True)
        """
        # Use composite key to separate global and room-scoped namespaces
        fs_key = f"global:{name}" if public else f"room:{self.room}:{name}"

        if fs_key in self._filesystems:
            scope = "global" if public else f"room '{self.room}'"
            raise ValueError(f"Filesystem '{name}' is already registered in {scope}.")

        # Validate that public filesystems require admin privileges
        if public and not self.is_admin:
            raise PermissionError(
                "Only admin users can register public filesystems. "
                "Please authenticate with admin credentials to register global filesystems."
            )

        # Store locally with composite key
        self._filesystems[fs_key] = {
            "fs": fs,
            "public": public,
            "name": name,
        }
        print(f"Registered filesystem '{name}' (type: {fs.__class__.__name__}).")

        scope = "global" if public else self.room
        print(f"Registering {'global' if public else 'room-scoped'} filesystem '{name}'...")

        # Register with server via Socket.IO
        worker_id = self.api.register_filesystem(
            name=name,
            fs_type=fs.__class__.__name__,
            socket_manager=self.socket,
            public=public,
        )

        # Store the worker_id assigned by server (server's request.sid)
        if worker_id:
            self._worker_id = worker_id

        print(
            f"Filesystem '{name}' registered with {scope} (worker_id: {self._worker_id})."
        )

    def log(self, message: str) -> dict | None:
        if not self.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        return self.socket.sio.call(
            "chat:message:create", {"content": message}, timeout=5
        )

    def edit_message(self, message_id: str, new_content: str) -> dict | None:
        if not self.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        return self.socket.sio.call(
            "chat:message:edit",
            {"messageId": message_id, "content": new_content},
            timeout=5,
        )

    def get_messages(
        self, limit: int = 30, before: int | None = None, after: int | None = None
    ) -> dict:
        return self.api.get_messages(limit=limit, before=before, after=after)

    def progress_tracker(self, description: str) -> ProgressTracker:
        """Create a progress tracker context manager for tracking long-running operations.

        Parameters
        ----------
        description : str
            Initial progress description

        Returns
        -------
        ProgressTracker
            A context manager that tracks operation progress

        Examples
        --------
        >>> with vis.progress_tracker("Loading data...") as tracker:
        ...     # Do some work
        ...     tracker.update(progress=50)
        ...     # Do more work
        ...     tracker.update(description="Processing...", progress=100)
        """
        return ProgressTracker(self, description)

    def _repr_html_(self):
        """Get an HTML representation for embedding the viewer in Jupyter notebooks.

        Returns
        -------
        IPython.display.IFrame
            An IFrame object displaying the ZnDraw viewer.

        """
        try:
            from IPython.display import IFrame
        except ImportError:
            raise ImportError(
                "IPython is required for viewer display. Install with: pip install ipython"
            )

        viewer_url = f"{self.url}/rooms/{self.room}/{self.user}"
        return IFrame(src=viewer_url, width="100%", height=600)._repr_html_()
