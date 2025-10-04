import dataclasses
import functools
import logging
import typing as t
from collections.abc import MutableSequence

import ase
import numpy as np

from zndraw.api_manager import APIManager
from zndraw.exceptions import LockError
from zndraw.extensions import Extension, ExtensionType
from zndraw.frame_cache import FrameCache
from zndraw.scene_manager import Geometries
from zndraw.settings import RoomConfig, settings
from zndraw.socket_manager import SocketIOLock, SocketManager
from zndraw.utils import atoms_from_dict, atoms_to_dict, update_colors_and_radii

log = logging.getLogger(__name__)


class _GeometryStore(t.TypedDict):
    type: str
    data: dict


class _ExtensionStore(t.TypedDict):
    public: bool
    run_kwargs: dict | None
    extension: t.Type[Extension]


class _TemplateValue:
    """Sentinel value for template parameter."""

    pass


@dataclasses.dataclass
class ZnDraw(MutableSequence):
    """A client for interacting with the ZnDraw server."""

    url: str = "http://localhost:5000"
    room: str = "default"
    user: str = "guest"
    auto_pickup_jobs: bool = True
    template: str | None | t.Type[_TemplateValue] = _TemplateValue

    _step: int = 0
    _len: int = 0
    _settings: dict = dataclasses.field(default_factory=dict, init=False)
    _extensions: dict[str, _ExtensionStore] = dataclasses.field(
        default_factory=dict, init=False
    )
    _client_id: str | None = dataclasses.field(default=None, init=False)
    _selection: frozenset[int] = frozenset()
    _frame_selection: frozenset[int] = frozenset()
    _bookmarks: dict[int, str] = dataclasses.field(default_factory=dict, init=False)
    _geometries: dict[str, _GeometryStore] = dataclasses.field(
        default_factory=dict, init=False
    )

    def __post_init__(self):
        self.api = APIManager(url=self.url, room=self.room, client_id=self._client_id)
        self.cache: FrameCache | None = FrameCache(maxsize=100)

        # Call join_room FIRST to get join token and prepare room
        response_data = self.api.join_room(template=self.template, user_id=self.user)

        # Update client_id if server assigned a new one
        if "clientId" in response_data:
            self._client_id = response_data["clientId"]
            self.api.client_id = self._client_id

        # Get join token for socket authentication
        join_token = response_data.get("joinToken")
        if not join_token:
            raise RuntimeError("Server did not provide join token")

        # Now create socket manager and connect WITH the token
        self.socket = SocketManager(zndraw_instance=self, join_token=join_token)
        self.connect()

        if response_data["selection"] is not None:
            self._selection = frozenset(response_data["selection"])
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
        return SocketIOLock(self.socket.sio, target="trajectory:meta")

    @property
    def geometries(self) -> Geometries:
        return Geometries(self)

    @property
    def sid(self) -> str:
        return self._client_id

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
    def selection(self) -> frozenset[int]:
        return self._selection

    @selection.setter
    def selection(self, value: t.Iterable[int] | None):
        indices = [] if value is None else list(value)
        if not all(
            isinstance(idx, int) and 0 <= idx < len(self[self.step]) for idx in indices
        ):
            raise ValueError("Selection must be an iterable of valid atom indices.")

        if self.socket.sio.connected:
            response = self.socket.sio.call(
                "selection:set", {"indices": indices}, timeout=5
            )
            if response and not response.get("success", False):
                raise RuntimeError(response.get("message", "Failed to set selection"))
            self._selection = frozenset(indices)
        else:
            raise RuntimeError("Client is not connected.")

    @property
    def frame_selection(self) -> frozenset[int]:
        return self._frame_selection

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
    def bookmarks(self) -> dict[int, str]:
        return self._bookmarks.copy()

    @bookmarks.setter
    def bookmarks(self, value: dict[int, str] | None):
        bookmarks = {} if value is None else value
        if not isinstance(bookmarks, dict):
            raise TypeError("Bookmarks must be a dictionary.")

        for idx, label in bookmarks.items():
            if not (isinstance(idx, int) and 0 <= idx < len(self)):
                raise IndexError(f"Bookmark index {idx} out of range.")
            if not isinstance(label, str):
                raise TypeError("Bookmark label must be a string.")

        if self.socket.sio.connected:
            bookmarks_json = {str(k): v for k, v in bookmarks.items()}
            response = self.socket.sio.call(
                "bookmarks:set", {"bookmarks": bookmarks_json}, timeout=5
            )
            if response and not response.get("success", False):
                raise RuntimeError(response.get("message", "Failed to set bookmarks"))
            self._bookmarks = bookmarks
        else:
            raise RuntimeError("Client is not connected.")

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
            self.api.submit_extension_settings(extension, self.user, data)

        for key in settings:
            if key not in self._settings:
                data = self.api.get_extension_settings(key, self.user)
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

        self.api.register_extension(
            name=name,
            category=extension.category,
            schema=extension.model_json_schema(),
            client_id=self.sid,
        )
        print(f"Extension '{name}' registered with room '{self.room}'.")
        self.socket._on_queue_update({})

    def run(self, extension: Extension) -> dict:
        return self.api.run_extension(
            category=extension.category.value,
            name=extension.__class__.__name__,
            user_id=self.user,
            data=extension.model_dump(),
        )

    def log(self, message: str) -> dict:
        if not self.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        return self.socket.sio.call(
            "chat:message:create", {"content": message}, timeout=5
        )

    def edit_message(self, message_id: str, new_content: str) -> dict:
        if not self.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        return self.socket.sio.call(
            "chat:message:edit",
            {"messageId": message_id, "content": new_content},
            timeout=5,
        )

    def get_messages(
        self, limit: int = 30, before: int = None, after: int = None
    ) -> dict:
        return self.api.get_messages(limit=limit, before=before, after=after)
