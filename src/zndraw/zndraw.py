import dataclasses
import functools
import logging
import typing as t
from collections.abc import MutableSequence

import ase
import numpy as np

from zndraw.api_manager import APIManager
from zndraw.bookmarks_manager import Bookmarks
from zndraw.exceptions import LockError
from zndraw.extensions import Extension, ExtensionType
from zndraw.figures_manager import Figures
from zndraw.frame_cache import FrameCache
from zndraw.metadata_manager import RoomMetadata
from zndraw.scene_manager import Geometries
from zndraw.server_manager import get_server_status
from zndraw.settings import RoomConfig, settings
from zndraw.socket_manager import SocketIOLock, SocketManager
from zndraw.utils import atoms_from_dict, atoms_to_dict, update_colors_and_radii
from zndraw.version_utils import validate_server_version

log = logging.getLogger(__name__)


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

    """

    url: str | None = None
    room: str = "default"
    user: str | None = None
    password: str | None = None
    auto_pickup_jobs: bool = True
    description: str | None = None
    copy_from: str | None = None

    _step: int = 0
    _len: int = 0
    _settings: dict = dataclasses.field(default_factory=dict, init=False)
    _extensions: dict[str, _ExtensionStore] = dataclasses.field(
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
    _lock: SocketIOLock | None = dataclasses.field(default=None, init=False)
    _metadata: RoomMetadata | None = dataclasses.field(default=None, init=False)

    def __post_init__(self):
        # Auto-discover local server if url is None
        if self.url is None:
            is_running, server_info, status_message = get_server_status()

            if is_running and server_info is not None:
                self.url = f"http://localhost:{server_info.port}"
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

        # Initialize the lock after socket is connected
        # TTL of 60 seconds with automatic refresh every 30 seconds
        self._lock = SocketIOLock(self.socket.sio, target="trajectory:meta", ttl=60)

        if response_data["frame_selection"] is not None:
            self._frame_selection = frozenset(response_data["frame_selection"])
        if response_data.get("bookmarks") is not None:
            self._bookmarks = {int(k): v for k, v in response_data["bookmarks"].items()}
        if response_data.get("step") is not None:
            self._step = int(response_data["step"])
        if response_data.get("geometries") is not None:
            self._geometries = response_data["geometries"]
        self._len = response_data["frameCount"]

    @property
    def lock(self) -> SocketIOLock:
        if self._lock is None:
            raise RuntimeError("Lock not initialized. Ensure client is connected.")
        return self._lock

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
        if self.socket.sio.connected:
            response = self.socket.sio.call(
                "set_frame_atomic", {"frame": value}, timeout=5
            )
            if response and not response.get("success", False):
                error_type = response.get("error")
                error_msg = response.get("message", "Failed to set frame")
                if error_type == "LockError":
                    raise LockError(error_msg)
                else:
                    raise RuntimeError(error_msg)
        else:
            raise RuntimeError("Client is not connected.")

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

    def _perform_locked_upload(self, action: str, data, **kwargs):
        if not self.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        with self.lock:
            return self.api.upload_frames(action, data, **kwargs)

    def append_frame(self, data: dict):
        self._perform_locked_upload("append", data)

    def extend_frames(self, data: list[dict]):
        result = self._perform_locked_upload("extend", data)
        return result.get("new_indices", [])

    def replace_frame(self, frame_id: int, data: dict):
        self._perform_locked_upload("replace", data, frame_id=frame_id)

    def insert_frame(self, index: int, data: dict):
        self._perform_locked_upload("insert", data, insert_position=index)

    def __len__(self) -> int:
        return self._len

    @t.overload
    def get(self, index: int, keys: list[str] | None = None) -> dict: ...
    @t.overload
    def get(self, index: list[int], keys: list[str] | None = None) -> list[dict]: ...
    @t.overload
    def get(self, index: slice, keys: list[str] | None = None) -> list[dict]: ...
    def get(self, index, keys: list[str] | None = None):
        """Get frame data as dictionaries, optionally with specific keys."""
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
                return self.cache.get(normalized_index)

            frame_data = self.api.get_frames([normalized_index], keys=keys)[0]
            if keys is None and self.cache is not None:
                self.cache.set(normalized_index, frame_data)
            return frame_data

        elif isinstance(index, slice):
            # Pass slice directly to API for efficiency
            indices = list(range(*index.indices(length)))
            if not indices:
                return []

            # Fetch frames using slice
            fetched_data = self.api.get_frames(index, keys=keys)

            # Cache individual frames if keys is None
            if keys is None and self.cache is not None:
                for idx, frame in zip(indices, fetched_data):
                    self.cache.set(idx, frame)

            return fetched_data

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
                    results_dict[idx] = self.cache.get(idx)
                elif idx not in misses:
                    misses.append(idx)

            if misses:
                fetched_data = self.api.get_frames(misses, keys=keys)
                for idx, frame in zip(misses, fetched_data):
                    if keys is None and self.cache is not None:
                        self.cache.set(idx, frame)
                    results_dict[idx] = frame

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
        data = self.get(index, keys=None)
        if isinstance(data, list):
            return [atoms_from_dict(d) for d in data]
        return atoms_from_dict(data)

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
                update_colors_and_radii(atom)
                dicts.append(atoms_to_dict(atom))
            self.set_frames(index, dicts)
        elif isinstance(atoms, ase.Atoms):
            update_colors_and_radii(atoms)
            self.set_frames(index, atoms_to_dict(atoms))
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
            self.replace_frame(index, value)
        elif isinstance(index, (slice, list)):
            with self.lock:
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

        with self.lock:
            self.api.delete_frames(index)

    def insert(self, index: int, atoms: ase.Atoms):
        if not isinstance(atoms, ase.Atoms):
            raise TypeError("Only ase.Atoms objects are supported")
        update_colors_and_radii(atoms)
        value = atoms_to_dict(atoms)
        if index < 0:
            index = max(0, len(self) + index + 1)
        elif index > len(self):
            index = len(self)
        self.insert_frame(index, value)

    def append(self, atoms: ase.Atoms):
        if not isinstance(atoms, ase.Atoms):
            raise TypeError("Only ase.Atoms objects are supported")
        update_colors_and_radii(atoms)
        self.append_frame(atoms_to_dict(atoms))

    def extend(self, atoms_list: t.Iterable[ase.Atoms]):
        dicts = []
        for atoms in atoms_list:
            if not isinstance(atoms, ase.Atoms):
                raise TypeError("Only ase.Atoms objects are supported")
            update_colors_and_radii(atoms)
            dicts.append(atoms_to_dict(atoms))
        self.extend_frames(dicts)

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
        if name in self._extensions:
            raise ValueError(f"Extension '{name}' is already registered.")
        if not hasattr(extension, "category") or extension.category not in [
            cat.value for cat in ExtensionType
        ]:
            raise ValueError("Extension must have a valid 'category' attribute.")

        self._extensions[name] = {
            "public": public,
            "run_kwargs": run_kwargs,
            "extension": extension,
        }
        print(f"Registered extension '{name}' of category '{extension.category}'.")

        if public:
            raise NotImplementedError("Public extensions are not supported yet.")

        worker_id = self.api.register_extension(
            name=name,
            category=extension.category,
            schema=extension.model_json_schema(),
            socket_manager=self.socket,
        )
        # Store the worker_id assigned by server (server's request.sid)
        if worker_id:
            self._worker_id = worker_id
        print(
            f"Extension '{name}' registered with room '{self.room}' (worker_id: {self._worker_id})."
        )
        self.socket._on_queue_update({})

    def run(self, extension: Extension) -> dict:
        return self.api.run_extension(
            category=extension.category.value,
            name=extension.__class__.__name__,
            data=extension.model_dump(),
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
