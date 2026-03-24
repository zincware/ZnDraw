"""ZnDraw main client class implementing MutableSequence[ase.Atoms]."""

from __future__ import annotations

import contextlib
import logging
import uuid
import warnings
from collections.abc import Generator, Iterable, MutableSequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Annotated, Any, Literal, cast, overload

import ase
import msgpack
import typing_extensions
from pydantic import SecretStr
from zndraw_joblib.client import ClaimedTask, Extension as JoblibExtension, JobManager

from zndraw.accessors import (
    Bookmarks,
    ChatMessages,
    Extensions,
    Figures,
    Geometries,
    Presets,
    RoomMetadata,
    Screenshots,
    SelectionGroups,
    Selections,
    Sessions,
    TaskHandle,
    Tasks,
)
from zndraw.client.api import APIManager
from zndraw.client.exceptions import NotConnectedError
from zndraw.client.lock import ZnDrawLock
from zndraw.client.serialization import (
    _MAX_CHUNK_FRAMES,
    _TARGET_CHUNK_BYTES,
    _estimate_frame_size,
    atoms_to_json_dict,
    raw_frame_to_atoms,
)
from zndraw.client.socket import SocketManager
from zndraw.geometries.camera import Camera

if TYPE_CHECKING:
    from zndraw.extensions.abc import Extension
    from zndraw.providers.frame_source import FrameSource
    from zndraw.tqdm import ZnDrawTqdm

log = logging.getLogger(__name__)


# =============================================================================
# ZnDraw - Main Client Class
# =============================================================================


@dataclass
class ZnDraw(MutableSequence[ase.Atoms]):
    """A synchronous client for interacting with the ZnDraw server.

    Implements MutableSequence for frame operations, providing list-like access
    to trajectory frames as ase.Atoms objects.

    Parameters
    ----------
    url : str | None
        URL of the ZnDraw server. If None, auto-discovers via PID file.
    room : str | None
        Room ID to connect to. If None, generates a random UUID.
    user : str | None
        User email for authentication. If None, creates a guest session.
    password : SecretStr | str | None
        Password for login. Accepts ``str`` (auto-wrapped to ``SecretStr``)
        or ``SecretStr``. Required when ``user`` is provided.
    copy_from : str | None
        Room ID to copy frames from when creating a new room, or an @-prefixed
        preset (``@empty`` for one empty frame, ``@none`` for zero frames).
        If None, uses the server's default template room. Default is ``@none``.

    Notes
    -----
    The Socket.IO connection is established lazily on the first call to
    ``mount()``, ``register_job()``, ``register_fs()``, or ``wait()``.
    All frame CRUD and metadata operations use REST and work without a socket.
    Call ``connect()`` explicitly if you need the socket immediately.

    Examples
    --------
    >>> import ase
    >>> vis = ZnDraw(url="http://localhost:8000")
    >>> len(vis)  # Number of frames
    >>> vis.append(ase.Atoms("H2O", positions=[[0,0,0], [1,0,0], [0,1,0]]))
    >>> atoms = vis[0]  # Get first frame as ase.Atoms
    >>> vis[0] = ase.Atoms("CO2")  # Update frame
    >>> del vis[0]  # Delete frame
    """

    url: str | None = None
    room: str | None = None
    user: str | None = None
    password: SecretStr | str | None = None
    token: str | None = None
    copy_from: str | None = "@none"
    create_if_missing: bool = True
    auto_pickup: bool = True
    polling_interval: float = 5.0
    heartbeat_interval: float = 30.0

    # Internal state
    api: APIManager = field(init=False)
    socket: SocketManager = field(init=False)
    _jobs: JobManager = field(init=False)
    cached_length: int | None = field(default=None, init=False, repr=False)
    mounted_source: FrameSource | None = field(default=None, init=False, repr=False)
    mounted_source_name: str | None = field(default=None, init=False, repr=False)

    # Accessors (lazy initialized)
    _selections: Selections | None = field(default=None, init=False)
    _selection_groups: SelectionGroups | None = field(default=None, init=False)
    _bookmarks: Bookmarks | None = field(default=None, init=False)
    _geometries: Geometries | None = field(default=None, init=False)
    _figures: Figures | None = field(default=None, init=False)
    _metadata: RoomMetadata | None = field(default=None, init=False)
    _presets: Presets | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        """Initialize the client (REST-only, socket connects lazily)."""
        import atexit

        from zndraw.auth_utils import resolve_token

        # Normalize password to SecretStr
        if isinstance(self.password, str):
            self.password = SecretStr(self.password)

        # Generate room ID if not provided
        if self.room is None:
            self.room = str(uuid.uuid4())

        # Resolve URL (auto-discover if None)
        self.url = self._resolve_url(self.url)

        # Create API manager (token may be None initially)
        self.api = APIManager(url=self.url, room_id=self.room, token=self.token)

        # Authenticate via shared resolve_token
        self.api.token = resolve_token(
            self.url,
            token=self.token,
            user=self.user,
            password=self.password,
        )

        # Populate self.user for guest/stored-token sessions
        # Skip when token was explicitly provided (caller already knows identity)
        if self.user is None and self.token is None:
            resp = self.api.http.get(
                "/v1/auth/users/me",
                headers={"Authorization": f"Bearer {self.api.token}"},
            )
            if resp.status_code == 200:
                self.user = resp.json().get("email")

        # Create socket manager (no connection yet -- connects lazily)
        self.socket = SocketManager(zndraw=self)

        # Create job manager (zero-cost until first register())
        self._jobs = JobManager(
            api=self.api,
            tsio=self.socket.tsio,
            execute=self._execute_task if self.auto_pickup else None,
            heartbeat_interval=self.heartbeat_interval,
            polling_interval=self.polling_interval,
        )

        # Verify/create room via REST and seed frame count cache
        try:
            info = self.api.get_room_info()
            self.cached_length = info.get("frame_count", 0)
        except KeyError:
            if not self.create_if_missing:
                raise
            self.api.create_room(copy_from=self.copy_from)
            self.cached_length = 0

        # Ensure cleanup on interpreter exit
        atexit.register(self.disconnect)

    # -------------------------------------------------------------------------
    # Connection Management
    # -------------------------------------------------------------------------

    def _ensure_socket_connected(self) -> None:
        """Connect the socket lazily if not already connected."""
        if not self.socket.connected:
            self.socket.connect()

    def connect(self) -> None:
        """Explicitly connect the Socket.IO connection.

        Most operations work over REST without a socket. Call this only
        if you need the socket immediately (e.g. to receive broadcasts).
        Socket-dependent methods like ``mount()``, ``register_job()``,
        and ``wait()`` call this automatically.
        """
        self._ensure_socket_connected()

    def disconnect(self) -> None:
        """Disconnect from the server. Idempotent.

        Worker cleanup is handled server-side by the ``on_disconnect``
        handler when the socket disconnects (clears providers, frame
        counts, etc.).  If the socket was never connected, there is no
        worker to clean up.
        """
        self.socket.disconnect()
        self.api.close()

    def wait(self) -> None:
        """Block until disconnected."""
        self._ensure_socket_connected()
        self.socket.wait()

    @property
    def connected(self) -> bool:
        """Check if the socket is connected."""
        return self.socket.connected

    def __enter__(self) -> ZnDraw:
        """Context manager entry."""
        return self

    def __exit__(self, *args: Any) -> None:
        """Context manager exit."""
        self.disconnect()

    # -------------------------------------------------------------------------
    # Source Mount
    # -------------------------------------------------------------------------

    def mount(self, source: FrameSource) -> None:
        """Mount a virtual frame source on this room.

        The room must be empty (``len(vis) == 0``) and have no existing mount.
        After mounting, the room becomes read-only -- frames are served on demand
        from the source via the provider system. Connects the socket if needed.

        Registers two providers:
        - ``FrameSourceRead`` (category "frames") -- serves individual frames.
        - ``FrameSourceLength`` (category "frames_meta") -- serves source length.

        Then PATCHes the room with the frame count so ``get_length()`` returns
        the correct value immediately.

        Parameters
        ----------
        source
            Any object satisfying the ``FrameSource`` protocol
            (e.g. ``asebytes.ASEIO``, ``list[ase.Atoms]``).
        """
        from zndraw.providers.frame_source import FrameSourceLength, FrameSourceRead

        try:
            length = len(source)
        except TypeError:
            raise TypeError(
                "Frame source does not support len(). "
                "For file-based backends (e.g. asebytes.ASEIO with XYZ), "
                "call db._backend.count_frames() before mounting."
            ) from None

        if self.mounted_source is not None:
            raise RuntimeError("Room already has a mount")
        if len(self) != 0:
            raise RuntimeError("Room must be empty before mounting")
        if self.room is None:
            raise NotConnectedError("Cannot mount: no room set")

        self._ensure_socket_connected()
        self.mounted_source = source
        # UUID ensures each mount gets fresh provider cache keys --
        # prevents stale results from a previous mount being served.
        mount_name = f"mount-{uuid.uuid4().hex[:12]}"
        self.mounted_source_name = mount_name
        self.register_provider(
            FrameSourceRead, name=mount_name, handler=source, room=self.room
        )
        self.register_provider(
            FrameSourceLength, name=mount_name, handler=source, room=self.room
        )
        self.api.update_room({"frame_count": length})
        self.cached_length = length

    def unmount(self) -> None:
        """Unmount the current source. Room returns to empty.

        Unregisters both providers and clears the external frame count.
        """
        if self.mounted_source is None:
            raise RuntimeError("No mount to remove")

        if self.room is None:
            raise NotConnectedError("Cannot unmount: no room set")
        mount_name = self.mounted_source_name
        self.jobs.unregister_provider(f"{self.room}:frames:{mount_name}")
        self.jobs.unregister_provider(f"{self.room}:frames_meta:{mount_name}")
        self.api.update_room({"frame_count": 0})
        self.mounted_source = None
        self.mounted_source_name = None
        self.cached_length = 0

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

    @property
    def step(self) -> int:
        """Current frame index."""
        return self.api.get_step()["step"]

    @step.setter
    def step(self, value: int) -> None:
        """Set current frame index."""
        if value < 0:
            value = len(self) + value
        self.api.update_step(value)

    @property
    def selections(self) -> Selections:
        """Access selections by geometry name."""
        if self._selections is None:
            self._selections = Selections(self.api)
        return self._selections

    @property
    def selection_groups(self) -> SelectionGroups:
        """Access named selection groups."""
        if self._selection_groups is None:
            self._selection_groups = SelectionGroups(self.api)
        return self._selection_groups

    @property
    def selection(self) -> tuple[int, ...]:
        """Get selection for 'particles' geometry (convenience property)."""
        return self.selections.get("particles", ())

    @selection.setter
    def selection(self, value: Iterable[int] | None) -> None:
        """Set selection for 'particles' geometry."""
        self.selections["particles"] = [] if value is None else list(value)

    @property
    def frame_selection(self) -> tuple[int, ...]:
        """Selected frame indices (sorted)."""
        indices = self.api.get_frame_selection()
        return tuple(sorted(indices)) if indices else ()

    @frame_selection.setter
    def frame_selection(self, value: Iterable[int] | None) -> None:
        """Set selected frame indices."""
        self.api.update_frame_selection(sorted(value) if value else [])

    @property
    def atoms(self) -> ase.Atoms:
        """Get the current frame as an ase.Atoms object."""
        return self[self.step]

    @atoms.setter
    def atoms(self, value: ase.Atoms) -> None:
        """Set the current frame from an ase.Atoms object."""
        self[self.step] = value

    @property
    def bookmarks(self) -> Bookmarks:
        """Access frame bookmarks."""
        if self._bookmarks is None:
            self._bookmarks = Bookmarks(self.api)
        return self._bookmarks

    @property
    def geometries(self) -> Geometries:
        """Access custom geometries."""
        if self._geometries is None:
            self._geometries = Geometries(self.api)
        return self._geometries

    @property
    def default_camera(self) -> str | None:
        """Default camera geometry key for new sessions."""
        return self.api.get_default_camera()

    @default_camera.setter
    def default_camera(self, value: str | None) -> None:
        """Set default camera. Validates key exists and is a Camera."""
        if value is not None:
            if value not in self.geometries:
                raise KeyError(f"Camera '{value}' not found in geometries")
            geom = self.geometries[value]
            if not isinstance(geom, Camera):
                raise TypeError(f"Geometry '{value}' is not a Camera")
        self.api.set_default_camera(value)

    @property
    def figures(self) -> Figures:
        """Access Plotly figures."""
        if self._figures is None:
            self._figures = Figures(self.api)
        return self._figures

    @property
    def presets(self) -> Presets:
        """Visual presets for this room."""
        if self._presets is None:
            self._presets = Presets(self.api)
        return self._presets

    @property
    def metadata(self) -> RoomMetadata:
        """Access room metadata."""
        if self._metadata is None:
            self._metadata = RoomMetadata(self.api)
        return self._metadata

    @property
    def sessions(self) -> Sessions:
        """Active frontend browser sessions in this room."""
        return Sessions(_api=self.api)

    @property
    def chat(self) -> ChatMessages:
        """Chat messages for this room.

        Read-only Sequence with ``.send()`` for posting messages.
        """
        return ChatMessages(self.api)

    @property
    def screenshots(self) -> Screenshots:
        """Screenshots for this room.

        Read-only Mapping where keys are screenshot IDs.
        """
        return Screenshots(self.api)

    @property
    def extensions(self) -> Extensions:
        """Available extensions for this room.

        Read-only Mapping where keys are full_names (e.g.
        ``'@internal:modifiers:Delete'``).
        """
        return Extensions(self.api)

    @property
    def tasks(self) -> Tasks:
        """Tasks submitted to this room.

        Read-only Mapping where keys are task IDs.
        Call as ``vis.tasks(status='running')`` for filtered views.
        """
        return Tasks(self.api)

    @property
    def locked(self) -> bool:
        """Whether the room is locked."""
        info = self.api.get_room_info()
        return info.get("locked", False)

    @locked.setter
    def locked(self, value: bool) -> None:
        self.api.update_room({"locked": value})

    def log(self, message: str) -> None:
        """Send a chat message to the room.

        .. deprecated::
            Use ``vis.chat.send(message)`` instead.
        """
        warnings.warn(
            "vis.log() is deprecated, use vis.chat.send() instead",
            DeprecationWarning,
            stacklevel=2,
        )
        self.chat.send(message)

    # -------------------------------------------------------------------------
    # Class Methods (server-level, no room needed)
    # -------------------------------------------------------------------------

    @staticmethod
    def _resolve_url(url: str | None) -> str:
        """Resolve the server URL from argument or PID file auto-discovery.

        Parameters
        ----------
        url
            Explicit URL, or None to auto-discover via PID file.

        Raises
        ------
        ConnectionError
            If no running server is found.
        """
        if url is not None:
            return url.rstrip("/")
        from zndraw.server_manager import find_running_server

        server_info = find_running_server()
        if server_info is not None:
            return f"http://localhost:{server_info.port}"
        raise ConnectionError(
            "No running zndraw server found."
            " Start one with `uv run zndraw` or pass `url`."
        )

    @classmethod
    def list_rooms(
        cls,
        url: str | None = None,
        *,
        token: str | None = None,
        search: str | None = None,
    ) -> list[dict[str, Any]]:
        """List all rooms on the server.

        Parameters
        ----------
        url
            Server URL. If None, auto-discovers via PID file.
        token
            JWT token. If None, uses stored token or creates guest session.
        search
            Optional search filter.
        """
        from zndraw.auth_utils import resolve_token

        resolved = cls._resolve_url(url)
        resolved_token = resolve_token(resolved, token=token)
        api = APIManager(url=resolved, room_id="", token=resolved_token)
        try:
            return api.list_rooms(search=search)
        finally:
            api.close()

    @classmethod
    def login(
        cls,
        url: str | None = None,
        username: str = "",
        password: str = "",
    ) -> str:
        """Authenticate and return a JWT token.

        Parameters
        ----------
        url
            Server URL. If None, auto-discovers via PID file.
        username
            User email.
        password
            User password.

        Returns
        -------
        str
            JWT access token.
        """
        from zndraw.auth_utils import resolve_token

        resolved = cls._resolve_url(url)
        return resolve_token(resolved, user=username, password=password)

    @property
    def jobs(self) -> JobManager:
        """Access the job manager for registering jobs and submitting tasks.

        Lazily connects the socket, since all job operations require it.
        """
        self._ensure_socket_connected()
        return self._jobs

    def _resolve_room(self, room: str | None) -> str:
        """Resolve room argument, defaulting to self.room."""
        resolved = room or self.room
        if resolved is None:
            raise NotConnectedError("No room set")
        return resolved

    def register_job(
        self,
        cls: type,
        *,
        room: Literal["@global"] | str | None = None,
        public: Annotated[
            bool | None,
            typing_extensions.deprecated("Use room='@global' instead of public=True"),
        ] = None,
    ) -> None:
        """Register an extension as a job. Connects the socket if needed.

        Parameters
        ----------
        cls
            Extension subclass to register.
        room
            Room scope. Use ``"@global"`` for global registration (admin-only).
            Defaults to ``self.room``.
        public
            .. deprecated::
                Use ``room='@global'`` instead.
        """
        if public is not None and room is not None:
            raise ValueError("Cannot specify both 'room' and 'public'")
        if public:
            warnings.warn(
                "public=True is deprecated, use room='@global' instead",
                DeprecationWarning,
                stacklevel=2,
            )
            room = "@global"
        elif room != "@global":
            room = self._resolve_room(room)

        self._ensure_socket_connected()
        self.jobs.register(cls, room=room)

    @typing_extensions.deprecated(
        "Use register_job(cls, room='@global') for global, "
        "or register_job(cls) for room-scoped"
    )
    def register_extension(
        self, cls: type, *, public: bool = False, **kwargs: Any
    ) -> None:
        """Register an extension.

        .. deprecated::
            Use :meth:`register_job` instead.
        """
        room = kwargs.pop("room", None)
        if kwargs:
            unexpected = ", ".join(sorted(kwargs))
            raise TypeError(f"Unexpected keyword argument(s): {unexpected}")
        if public and room is not None:
            raise ValueError("Cannot specify both 'room' and 'public'")
        room = "@global" if public else room
        self.register_job(cls, room=room)

    def register_provider(
        self,
        provider_cls: type,
        *,
        name: str,
        handler: Any,
        room: str | None = None,
    ) -> uuid.UUID:
        """Register a provider for serving read requests. Connects the socket if needed.

        Parameters
        ----------
        provider_cls
            Provider subclass defining category and read schema.
        name
            Unique name for this provider instance.
        handler
            Object passed to ``provider.read(handler)`` on dispatch.
        room
            Room scope. Defaults to ``self.room``.
        """
        self._ensure_socket_connected()
        return self.jobs.register_provider(
            provider_cls, name=name, handler=handler, room=self._resolve_room(room)
        )

    def register_fs(self, fs: Any, *, name: str, room: str | None = None) -> None:
        """Register an fsspec filesystem as a provider with file loading.

        Parameters
        ----------
        fs
            An fsspec filesystem instance (e.g. ``fsspec.filesystem("file")``).
        name
            Unique name for this filesystem (e.g. "local", "s3-bucket").
        room
            Room scope. Defaults to ``self.room``.
        """
        from zndraw.extensions.filesystem import LoadFile
        from zndraw.providers.filesystem import FilesystemRead

        r = self._resolve_room(room)
        self.register_provider(FilesystemRead, name=name, handler=fs, room=r)
        self.register_job(LoadFile, room=r)

    # -------------------------------------------------------------------------
    # Locking
    # -------------------------------------------------------------------------

    def get_lock(self, msg: str | None = None) -> ZnDrawLock:
        """Get an edit lock context manager.

        Parameters
        ----------
        msg : str | None
            Optional message describing the lock purpose.

        Returns
        -------
        ZnDrawLock
            Lock context manager.

        Examples
        --------
        >>> with vis.get_lock(msg="Uploading frames"):
        ...     vis.extend(frames)
        """
        return ZnDrawLock(api=self.api, msg=msg)

    # -------------------------------------------------------------------------
    # MutableSequence Implementation (Frame Operations)
    # -------------------------------------------------------------------------

    def __len__(self) -> int:
        """Return number of frames.

        When the socket is connected, uses a cached value (updated by
        ``FramesInvalidate`` events and local mutations).  Without socket,
        always queries the server to avoid stale reads.
        """
        if self.cached_length is not None and self.socket.connected:
            return self.cached_length
        info = self.api.get_room_info()
        length = info.get("frame_count", 0)
        self.cached_length = length
        return length

    @overload
    def __getitem__(self, index: int) -> ase.Atoms: ...

    @overload
    def __getitem__(self, index: slice) -> list[ase.Atoms]: ...

    def __getitem__(self, index: int | slice) -> ase.Atoms | list[ase.Atoms]:
        """Get frame(s) by index or slice as ase.Atoms objects."""
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            frame = self.api.get_frame(index)
            return raw_frame_to_atoms(frame)

        if isinstance(index, slice):
            length = len(self)
            start, stop, step = index.indices(length)
            if step == 1:
                frames = self.api.get_frames(start=start, stop=stop)
            else:
                indices = list(range(start, stop, step))
                frames = self.api.get_frames(indices=indices) if indices else []
            return [raw_frame_to_atoms(f) for f in frames]

        raise TypeError(f"Index must be int or slice, not {type(index).__name__}")

    @overload
    def __setitem__(self, index: int, value: ase.Atoms) -> None: ...

    @overload
    def __setitem__(self, index: slice, value: Iterable[ase.Atoms]) -> None: ...

    def __setitem__(  # type: ignore[override]
        self, index: int | slice, value: ase.Atoms | Iterable[ase.Atoms]
    ) -> None:
        """Set frame(s) at given index from ase.Atoms objects."""
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            if not isinstance(value, ase.Atoms):
                raise TypeError("Value must be an ase.Atoms object")
            self.api.update_frame(index, atoms_to_json_dict(value))

        elif isinstance(index, slice):
            if isinstance(value, ase.Atoms):
                raise TypeError(
                    "Value must be iterable of ase.Atoms for slice assignment"
                )
            value_list = list(value)
            if not all(isinstance(v, ase.Atoms) for v in value_list):
                raise TypeError("All values must be ase.Atoms objects")
            length = len(self)
            start, stop, step = index.indices(length)
            indices = list(range(start, stop, step))

            if step != 1 and len(value_list) != len(indices):
                raise ValueError(
                    f"attempt to assign sequence of size {len(value_list)} "
                    f"to extended slice of size {len(indices)}"
                )

            for i, idx in enumerate(indices):
                if i < len(value_list):
                    self.api.update_frame(idx, atoms_to_json_dict(value_list[i]))

        else:
            raise TypeError(f"Index must be int or slice, not {type(index).__name__}")

    def __delitem__(self, index: int | slice) -> None:
        """Delete frame(s) at given index."""
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            self.api.delete_frame(index)
            self.cached_length = None

        elif isinstance(index, slice):
            length = len(self)
            start, stop, step = index.indices(length)
            indices = sorted(range(start, stop, step), reverse=True)
            for idx in indices:
                self.api.delete_frame(idx)
            self.cached_length = None

        else:
            raise TypeError(f"Index must be int or slice, not {type(index).__name__}")

    def insert(self, _index: int, value: ase.Atoms) -> None:
        """Insert an ase.Atoms frame at the given index.

        Note: Due to API limitations, this appends and then reorders.
        """
        if not isinstance(value, ase.Atoms):
            raise TypeError("Value must be an ase.Atoms object")
        # For simplicity, just append (true insert requires server support)
        result = self.api.append_frames([atoms_to_json_dict(value)])
        self.cached_length = result.get("total")

    def append(self, value: ase.Atoms) -> None:
        """Append an ase.Atoms frame to the end."""
        if not isinstance(value, ase.Atoms):
            raise TypeError("Value must be an ase.Atoms object")
        result = self.api.append_frames([atoms_to_json_dict(value)])
        self.cached_length = result.get("total")

    @contextlib.contextmanager
    @typing_extensions.deprecated(
        "Use ZnDrawTqdm directly: "
        "for x in ZnDrawTqdm(items, vis=vis, description='...', unit='it'): ..."
    )
    def progress_bar(
        self,
        iterable: Iterable[Any] | None = None,
        *,
        total: int | None = None,
        description: str = "Processing...",
        unit: str = "it",
    ) -> Generator[ZnDrawTqdm, None, None]:
        """Context manager yielding a ZnDrawTqdm progress bar.

        .. deprecated::
            Use :class:`ZnDrawTqdm` directly instead.

        Parameters
        ----------
        iterable : Iterable[Any] | None
            Optional iterable to wrap.
        total : int | None
            Total expected iterations.
        description : str
            Label shown in the UI.
        unit : str
            Unit label (e.g. ``"frames"``).

        Yields
        ------
        ZnDrawTqdm
            A tqdm-compatible progress bar.
        """
        from zndraw.tqdm import ZnDrawTqdm

        pbar = ZnDrawTqdm(
            iterable, total=total, vis=self, description=description, unit=unit
        )
        try:
            yield pbar
        finally:
            pbar.close()

    def extend(self, values: Iterable[ase.Atoms]) -> None:
        """Extend with multiple ase.Atoms frames.

        Streams frames in size-targeted chunks (~2 MB, max 1000 frames each).
        Shows a tqdm progress bar in the terminal and broadcasts progress to
        the ZnDraw UI via ``ZnDrawTqdm``. Accepts any iterable including
        generators.
        """
        from zndraw.tqdm import ZnDrawTqdm

        try:
            total_frames: int | None = len(values)  # type: ignore[arg-type]
        except TypeError:
            total_frames = None

        chunk: list[dict[str, Any]] = []
        chunk_size = 0
        result: dict[str, Any] = {}

        progress = ZnDrawTqdm(
            total=total_frames,
            vis=self,
            description="Uploading frames",
            unit="frames",
        )

        try:
            for atoms in values:
                if not isinstance(atoms, ase.Atoms):
                    raise TypeError("All values must be ase.Atoms objects")
                frame = atoms_to_json_dict(atoms)
                frame_size = _estimate_frame_size(frame)

                if chunk_size > 0 and (
                    chunk_size + frame_size > _TARGET_CHUNK_BYTES
                    or len(chunk) >= _MAX_CHUNK_FRAMES
                ):
                    result = self.api.append_frames(chunk)
                    progress.update(len(chunk))
                    chunk = []
                    chunk_size = 0

                chunk.append(frame)
                chunk_size += frame_size

            # Flush remaining frames
            if chunk:
                result = self.api.append_frames(chunk)
                progress.update(len(chunk))
        finally:
            progress.close()

        if result:
            self.cached_length = result.get("total")

    @staticmethod
    def _decode_raw_frame(frame: dict[bytes, bytes]) -> dict[str, Any]:
        """Decode a raw msgpack frame to a dict with string keys and decoded values."""
        import msgpack_numpy

        return {
            k.decode(): msgpack.unpackb(v, object_hook=msgpack_numpy.decode)
            for k, v in frame.items()
        }

    def get(
        self,
        index: int | list[int] | slice,
        keys: list[str] | None = None,
    ) -> dict[str, Any] | list[dict[str, Any]]:
        """Get decoded frame data with optional key filtering.

        Returns decoded dictionaries with string keys and Python/numpy values,
        transferring only the requested keys from the server.

        Parameters
        ----------
        index : int | list[int] | slice
            Frame index, list of indices, or slice.
        keys : list[str] | None
            If provided, only return these keys from each frame
            (e.g. ``["info.energy", "arrays.positions"]``).

        Returns
        -------
        dict[str, Any] | list[dict[str, Any]]
            Single dict for int index, list of dicts for slice/list index.
        """
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            frames = self.api.get_frames(indices=[index], keys=keys)
            return self._decode_raw_frame(frames[0]) if frames else {}

        if isinstance(index, list):
            length = len(self)
            normalized = [i if i >= 0 else length + i for i in index]
            frames = self.api.get_frames(indices=normalized, keys=keys)
        elif isinstance(index, slice):
            length = len(self)
            start, stop, step = index.indices(length)
            if step == 1:
                frames = self.api.get_frames(start=start, stop=stop, keys=keys)
            else:
                indices = list(range(start, stop, step))
                frames = self.api.get_frames(indices=indices, keys=keys)
        else:
            raise TypeError(
                f"Index must be int, list, or slice, not {type(index).__name__}"
            )

        return [self._decode_raw_frame(f) for f in frames]

    def set_frames(
        self,
        index: int | slice | list[int],
        value: ase.Atoms | list[ase.Atoms],
    ) -> None:
        """Set frame(s) at given index from ase.Atoms objects.

        This is an alternative to __setitem__ that supports list indices.
        """
        if isinstance(index, int):
            if not isinstance(value, ase.Atoms):
                raise TypeError("Value must be ase.Atoms for single index")
            self[index] = value
        elif isinstance(index, slice):
            if isinstance(value, ase.Atoms):
                raise TypeError("Value must be list of ase.Atoms for slice index")
            self[index] = value
        elif isinstance(index, list):
            if not isinstance(value, list):
                raise TypeError("Value must be list of ase.Atoms for list index")
            if not all(isinstance(v, ase.Atoms) for v in value):
                raise TypeError("All values must be ase.Atoms objects")
            # For list indices, update each frame individually
            if len(index) != len(value):
                raise ValueError("Index and value lists must have same length")
            for idx, atoms in zip(index, value, strict=True):
                self.api.update_frame(idx, atoms_to_json_dict(atoms))
        else:
            raise TypeError(
                f"Index must be int, slice, or list, not {type(index).__name__}"
            )

    # -------------------------------------------------------------------------
    # Extension Execution
    # -------------------------------------------------------------------------

    @overload
    def run(self, extension: str, **kwargs: Any) -> TaskHandle: ...

    @overload
    def run(
        self, extension: Extension, *, job_room: str | None = None
    ) -> TaskHandle: ...

    def run(
        self,
        extension: str | Extension,
        *,
        job_room: str | None = None,
        **kwargs: Any,
    ) -> TaskHandle:
        """Submit an extension for execution.

        Parameters
        ----------
        extension : str | Extension
            Either a fully-qualified name (e.g. ``"@internal:modifiers:Delete"``)
            for pure-REST dispatch, or an ``Extension`` instance for the
            existing job-system path.
        job_room : str | None
            Only used with ``Extension`` instances.  If *None*, auto-detects:
            ``@internal`` for built-in extensions, ``@global`` otherwise.
        **kwargs : Any
            Only used with string dispatch -- passed as the task payload.

        Returns
        -------
        TaskHandle
            Handle for tracking task progress.
        """
        if isinstance(extension, str):
            if self.room is None:
                raise NotConnectedError("Cannot run: no room set")
            task_id = self.api.submit_task(extension, kwargs)["id"]
            return TaskHandle(id=task_id, _api=self.api)

        # Extension instance path
        if job_room is None:
            module = extension.__class__.__module__
            job_room = (
                "@internal" if module.startswith("zndraw.extensions") else "@global"
            )

        if self.room is None:
            raise NotConnectedError("Cannot run: no room set")
        task_id = self.jobs.submit(
            cast("JoblibExtension", extension), room=self.room, job_room=job_room
        )
        return TaskHandle(id=task_id, _api=self.api)

    # -------------------------------------------------------------------------
    # Worker Serve Loop
    # -------------------------------------------------------------------------

    def _execute_task(self, task: ClaimedTask) -> None:
        """Execute a claimed task's extension logic.

        Lifecycle management (start/complete/fail) is handled by
        ``JobManager._claim_loop`` -- this callback only runs the extension.
        """
        task_vis = ZnDraw(
            url=self.url,
            room=task.room_id,
            token=self.api.token,
        )
        try:
            task.extension.run(task_vis, providers=self.jobs.handlers)
        finally:
            task_vis.disconnect()
