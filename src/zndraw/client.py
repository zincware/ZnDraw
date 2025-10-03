import dataclasses
import functools
import logging
import typing as t
import uuid
import warnings
from collections.abc import MutableSequence

import msgpack
import numpy as np
import requests
import socketio
from cachetools import LRUCache
from tqdm import tqdm

from zndraw.exceptions import LockError
from zndraw.extensions import Extension, ExtensionType
from zndraw.settings import RoomConfig, settings
from zndraw.storage import decode_data, encode_data

log = logging.getLogger(__name__)


@dataclasses.dataclass
class SocketIOLock:
    """A client-side context manager for a distributed lock via Socket.IO."""

    sio: socketio.Client
    target: str

    def acquire(self, timeout: float = 60) -> bool:
        """
        Acquire a lock for the specific target.
        Waits for the server's confirmation.
        """
        payload = {"target": self.target}
        # sio.call is inherently blocking, so it waits for the server's response.
        response = self.sio.call("lock:acquire", payload, timeout=timeout)
        return response and response.get("success", False)

    def release(self) -> bool:
        """Release the lock."""
        payload = {"target": self.target}
        response = self.sio.call("lock:release", payload, timeout=10)
        return response and response.get("success", False)

    def __enter__(self):
        if not self.acquire():
            raise RuntimeError(f"Failed to acquire lock for target '{self.target}'")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.release():
            # Use a warning here because an exception in release could hide
            # the original exception that occurred inside the 'with' block.
            warnings.warn(
                f"Failed to release lock for target '{self.target}'. It may have expired."
            )


class _ExtensionStore(t.TypedDict):
    public: bool
    run_kwargs: dict | None
    extension: t.Type[Extension]


class _TemplateValue:
    """Sentinel value for template parameter."""

    pass


@dataclasses.dataclass
class Client(MutableSequence):
    """A client for interacting with the ZnDraw server. Implements MutableSequence for frame operations."""

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
    _client_id: str = dataclasses.field(
        default_factory=lambda: str(uuid.uuid4()), init=False
    )
    _selection: frozenset[int] = frozenset()
    _frame_selection: frozenset[int] = frozenset()
    _bookmarks: dict[int, str] = dataclasses.field(default_factory=dict, init=False)
    _lock: SocketIOLock = dataclasses.field(init=False)
    _cache: LRUCache[int, dict] = dataclasses.field(
        default_factory=lambda: LRUCache(maxsize=100), init=False
    )

    def __post_init__(self):
        self.sio = socketio.Client()
        self.sio.on("connect", self._on_connect)
        self.sio.on("frame_update", self._on_frame_update)
        self.sio.on("selection:update", self._on_selection_update)
        self.sio.on("len_frames", self._on_len_frames_update)
        self.sio.on("invalidate", self._on_invalidate)
        self.sio.on("queue:update", self._on_queue_update)
        self.sio.on("frame_selection:update", self._on_frame_selection_update)
        self.sio.on("bookmarks:update", self._on_bookmarks_update)
        self.sio.on("frames:invalidate", self._on_frames_invalidate)

        self._lock = SocketIOLock(self.sio, target="trajectory:meta")

        if self.template is _TemplateValue:
            response = requests.post(f"{self.url}/api/rooms/{self.room}/join", json={})
        else:
            response = requests.post(
                f"{self.url}/api/rooms/{self.room}/join",
                json={"template": self.template},
            )

        if response.status_code != 200:
            raise RuntimeError(
                f"Failed to join room '{self.room}': {response.status_code} {response.text}"
            )
        response_data = response.json()
        if response_data["selection"] is not None:
            self._selection = frozenset(response_data["selection"])
        if response_data["frame_selection"] is not None:
            self._frame_selection = frozenset(response_data["frame_selection"])
        if response_data.get("bookmarks") is not None:
            # Convert string keys to int keys
            self._bookmarks = {int(k): v for k, v in response_data["bookmarks"].items()}
        self._len = response_data["frameCount"]

    @property
    def lock(self) -> SocketIOLock:
        """Get a lock object for the trajectory metadata."""
        return self._lock

    @property
    def sid(self) -> str:
        """Get the client's unique identifier (distinct from Socket.IO SIDs)."""
        return self._client_id

    def _on_frame_update(self, data):
        """Internal callback for when a frame update is received."""
        if "frame" in data:
            self._step = data["frame"]

    def _on_len_frames_update(self, data):
        """Internal callback for when a len_frames update is received."""
        if "count" in data:
            self._len = data["count"]

    def _on_selection_update(self, data):
        """Internal callback for when a selection update is received."""
        if "indices" in data:
            self._selection = frozenset(data["indices"])

    def _on_frame_selection_update(self, data):
        """Internal callback for when a frame selection update is received."""
        if "indices" in data:
            self._frame_selection = frozenset(data["indices"])

    def _on_bookmarks_update(self, data):
        """Internal callback for when bookmarks are updated."""
        if "bookmarks" in data:
            # Convert string keys to int keys
            self._bookmarks = {int(k): v for k, v in data["bookmarks"].items()}

    def _on_queue_update(self, data: dict):
        print(f"Queue update received: {data}")
        if not self.auto_pickup_jobs:
            return
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/jobs/next", json={"workerId": self.sid}
        )
        if response.status_code == 200:
            data = response.json()
            if "jobId" in data:
                try:
                    self._on_task_run(
                        data=data.get("data"),
                        extension=data.get("extension"),
                        category=data.get("category"),
                    )
                    response = requests.put(
                        f"{self.url}/api/rooms/{self.room}/jobs/{data.get('jobId')}/status",
                        json={"status": "completed", "workerId": self.sid},
                    )
                    if response.status_code != 200:
                        log.error(
                            f"Failed to mark job {data.get('jobId')} as complete: {response.status_code} {response.text}"
                        )
                    # log as completed
                except Exception as e:
                    log.error(f"Error processing job {data.get('jobId')}: {e}")
                    response = requests.put(
                        f"{self.url}/api/rooms/{self.room}/jobs/{data.get('jobId')}/status",
                        json={
                            "status": "failed",
                            "error": str(e),
                            "workerId": self.sid,
                        },
                    )
                # we have finished now, so we check for more jobs
                self._on_queue_update({})
        elif response.status_code == 400:
            print("No jobs available.")
        else:
            log.error(
                f"Failed to fetch next job: {response.status_code} {response.text}"
            )

    def _on_task_run(self, data: dict, extension: str, category: str):
        """Internal callback for when a task is run."""

        ext = self._extensions[extension]["extension"]
        instance = ext(**(data))
        # TODO: if public=True, we need to create a new vis instance, connected to the correct room
        instance.run(self, **(self._extensions[extension]["run_kwargs"] or {}))

    def _on_invalidate(self, data: dict):
        """Internal callback for when settings are invalidated."""
        # "invalidate", {"userId": user_id, "category": category, "extension": extension, "roomId": room_id}, to=f"user:{user_id}"
        if data["category"] == "settings":
            self._settings.pop(data["extension"], None)

    def _on_frames_invalidate(self, data: dict):
        """
        Callback to invalidate the client-side cache based on server events.
        """
        log.debug(f"Received cache invalidation event: {data}")
        operation = data.get("operation")

        if operation == "replace":
            # This now correctly handles ONLY single-item replacements
            idx = data.get("affectedIndex")
            if idx is not None:
                self._cache.pop(idx, None)

        elif operation in ("insert", "delete", "bulk_replace"):
            # An insert, delete, or bulk slice operation has shifted subsequent frames.
            # Invalidate all cache entries from the point of change onwards.
            from_idx = data.get("affectedFrom")
            if from_idx is not None:
                new_cache = LRUCache(maxsize=self._cache.maxsize)
                for k, v in self._cache.items():
                    if k < from_idx:
                        new_cache[k] = v
                self._cache = new_cache


        elif operation == "clear_all":
            self._cache.clear()

        else:
            log.warning(
                f"Unknown or broad invalidation event received. Clearing entire frame cache."
            )
            self._cache.clear()
    @property
    def step(self) -> int:
        """Get the current step/frame index."""
        return self._step

    @step.setter
    def step(self, value: int):
        """Set the current step/frame index and notify the server."""
        if not isinstance(value, int) or value < 0:
            raise ValueError("Step must be a non-negative integer.")
        if value >= self._len:
            raise ValueError(
                f"Step {value} is out of bounds. Current number of frames: {self._len}."
            )
        self._step = value
        if self.sio.connected:
            response = self.sio.call("set_frame_atomic", {"frame": value}, timeout=5)
            if response and not response.get("success", False):
                error_type = response.get("error")
                error_msg = response.get("message", "Failed to set frame")
                if error_type == "LockError":
                    raise LockError(error_msg)
                else:
                    raise RuntimeError(error_msg)
        else:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

    @property
    def selection(self) -> frozenset[int]:
        """Get the current selection of frame indices."""
        return self._selection

    @selection.setter
    def selection(self, value: t.Iterable[int] | None):
        """Set the current selection of frame indices."""
        if value is None:
            indices = []
        else:
            indices = [x for x in value]
            if not all(
                isinstance(idx, int) and 0 <= idx < len(self[self.step])
                for idx in indices
            ):
                raise ValueError(
                    "Selection must be an iterable of valid frame indices."
                )
        if self.sio.connected:
            response = self.sio.call("selection:set", {"indices": indices}, timeout=5)
            if response and not response.get("success", False):
                error_type = response.get("error")
                error_msg = response.get("message", "Failed to set selection")
                if error_type == "LockError":
                    raise LockError(error_msg)
                else:
                    raise RuntimeError(error_msg)
            self._selection = frozenset(indices)
        else:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

    @property
    def frame_selection(self) -> frozenset[int]:
        """Get the current selection of frame indices."""
        return self._frame_selection

    @frame_selection.setter
    def frame_selection(self, value: t.Iterable[int] | None):
        """Set the current selection of frame indices."""
        if value is None:
            indices = []
        else:
            indices = [x for x in value]
            if not all(
                isinstance(idx, int) and 0 <= idx < len(self) for idx in indices
            ):
                raise ValueError(
                    "Selection must be an iterable of valid frame indices."
                )
        if self.sio.connected:
            response = self.sio.call(
                "frame_selection:set", {"indices": indices}, timeout=5
            )
            if response and not response.get("success", False):
                error_type = response.get("error")
                error_msg = response.get("message", "Failed to set selection")
                if error_type == "LockError":
                    raise LockError(error_msg)
                else:
                    raise RuntimeError(error_msg)
            self._frame_selection = frozenset(indices)
        else:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

    @property
    def bookmarks(self) -> dict[int, str]:
        """Get the current bookmarks mapping frame indices to labels."""
        return self._bookmarks.copy()

    @bookmarks.setter
    def bookmarks(self, value: dict[int, str] | None):
        """Set bookmarks for frames.

        Args:
            value: Dictionary mapping frame indices to bookmark labels, or None to clear all bookmarks
        """
        if value is None:
            bookmarks = {}
        else:
            if not isinstance(value, dict):
                raise TypeError("Bookmarks must be a dictionary")

            bookmarks = {}
            for idx, label in value.items():
                if not isinstance(idx, int) or idx < 0 or idx >= len(self):
                    raise IndexError(
                        f"Bookmark index {idx} out of range [0, {len(self) - 1}]"
                    )
                if not isinstance(label, str):
                    raise TypeError(
                        f"Bookmark label must be a string, got {type(label).__name__}"
                    )
                bookmarks[idx] = label

        if self.sio.connected:
            # Convert int keys to strings for JSON
            bookmarks_json = {str(k): v for k, v in bookmarks.items()}
            response = self.sio.call(
                "bookmarks:set", {"bookmarks": bookmarks_json}, timeout=5
            )
            if response and not response.get("success", False):
                error_type = response.get("error")
                error_msg = response.get("message", "Failed to set bookmarks")
                if error_type == "LockError":
                    raise LockError(error_msg)
                else:
                    raise RuntimeError(error_msg)
            self._bookmarks = bookmarks
        else:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

    def connect(self):
        """Establishes a connection to the server."""
        if self.sio.connected:
            print("Already connected.")
            return
        self.sio.connect(self.url)

    def disconnect(self):
        """Disconnects from the server."""
        if self.sio.connected:
            self.sio.disconnect()
            print("Disconnected.")

    def _on_connect(self):
        """Internal callback for when a connection is established."""
        # Prepare join_room data
        join_data = {
            "room": self.room,
            "userId": self.user,
            "clientId": self._client_id,
        }

        self.sio.emit("join_room", join_data)
        log.debug(f"Joined room: '{self.room}' with client ID: {self._client_id}")
        for name, ext in self._extensions.items():
            response = requests.post(
                f"{self.url}/api/rooms/{self.room}/extensions/register",
                json={
                    "name": name,
                    "category": ext["extension"].category,
                    "schema": ext["extension"].model_json_schema(),
                    "clientId": self.sid,
                },
            )
            response.raise_for_status()
        # check for jobs
        self._on_queue_update({})

    def get_frame(self, frame_id: int, keys: list[str] | None = None) -> dict:
        """Fetches a single frame's data from the server."""
        # Handle negative indices
        if frame_id < 0:
            frame_id = self.len_frames() + frame_id
        return self.get_frames([frame_id], keys=keys)[0]

    def _perform_locked_upload(self, action: str, data, **kwargs):
        """Perform a locked upload operation with common error handling."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            try:
                serialized_data = (
                    encode_data(data)
                    if isinstance(data, dict)
                    else [encode_data(frame) for frame in data]
                )
                packed_data = msgpack.packb(serialized_data)

                upload_url = f"{self.url}/api/rooms/{self.room}/frames"
                params = {"action": action}
                params.update(kwargs)

                http_response = requests.post(
                    upload_url, data=packed_data, params=params, timeout=30
                )

                # Check for errors
                if http_response.status_code == 404:
                    try:
                        error_data = http_response.json()
                        error_type = error_data.get("type", "")
                        error_msg = error_data.get("error", http_response.text)

                        if error_type == "IndexError":
                            raise IndexError(error_msg)
                    except ValueError:
                        pass

                http_response.raise_for_status()
                result = http_response.json()

                if not result.get("success"):
                    raise RuntimeError(
                        f"Server reported failure: {result.get('error')}"
                    )

                return result
            except requests.exceptions.RequestException as e:
                raise RuntimeError(f"Error uploading frame data: {e}") from e

    def append_frame(self, data: dict):
        """Appends a single frame to the trajectory."""
        self._perform_locked_upload("append", data)

    def extend_frames(self, data: list[dict]):
        """
        Extends the trajectory by adding multiple frames in a single operation.
        Uses a single lock for the entire operation to ensure atomicity.

        Args:
            data: List of dictionaries, each containing data for one frame
        """
        result = self._perform_locked_upload("extend", data)
        return result.get("new_indices", [])

    def get_frames(self, indices_or_slice, keys: list[str] | None = None) -> list[dict]:
        """
        Fetches multiple frames' data from the server in a single call.

        Args:
            indices_or_slice: Either a list of frame indices [0, 2, 5] or a slice object slice(start, stop, step)
            keys: Optional list of keys to retrieve for each frame

        Returns:
            List of dictionaries, each containing data for one frame
        """
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        # Prepare request payload based on input type
        if isinstance(indices_or_slice, list):
            # Direct list of indices - convert to comma-separated string
            payload = {"indices": ",".join(str(i) for i in indices_or_slice)}
        elif isinstance(indices_or_slice, slice):
            # Slice object - extract start, stop, step
            payload = {}
            if indices_or_slice.start is not None:
                payload["start"] = indices_or_slice.start
            if indices_or_slice.stop is not None:
                payload["stop"] = indices_or_slice.stop
            if indices_or_slice.step is not None:
                payload["step"] = indices_or_slice.step
        else:
            raise ValueError(
                "indices_or_slice must be either a list of integers or a slice object"
            )

        # Add keys parameter if specified (comma-separated)
        if keys is not None:
            payload["keys"] = ",".join(keys)

        full_url = f"{self.url}/api/rooms/{self.room}/frames"
        response = requests.get(full_url, params=payload, timeout=30)

        # Check for errors
        if response.status_code == 404:
            try:
                error_data = response.json()
                error_type = error_data.get("type", "")
                error_msg = error_data.get("error", response.text)

                if error_type == "KeyError":
                    raise KeyError(error_msg)
                elif error_type == "IndexError":
                    raise IndexError(error_msg)
            except ValueError:
                # JSON parsing failed, let raise_for_status handle it
                pass

        response.raise_for_status()

        serialized_frames = msgpack.unpackb(response.content, strict_map_key=False)
        return [decode_data(frame) for frame in serialized_frames]

    def len_frames(self) -> int:
        """Returns the number of frames in the current room."""
        return self._len

    def delete_frame(self, index: int | slice | list[int]):
        """Deletes frame(s) from the current room."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            delete_url = f"{self.url}/api/rooms/{self.room}/frames"
            params = {"action": "delete"}

            # Determine if this is a single frame or batch delete
            if isinstance(index, int):
                # Single frame deletion
                params["frame_id"] = index
            elif isinstance(index, list):
                # List of indices - convert to comma-separated string
                params["indices"] = ",".join(str(i) for i in index)
            elif isinstance(index, slice):
                # Slice object - extract start, stop, step
                if index.start is not None:
                    params["start"] = index.start
                if index.stop is not None:
                    params["stop"] = index.stop
                if index.step is not None:
                    params["step"] = index.step

            response = requests.delete(delete_url, params=params, timeout=30)

            # Check for errors
            if response.status_code == 404:
                try:
                    error_data = response.json()
                    error_type = error_data.get("type", "")
                    error_msg = error_data.get("error", response.text)

                    if error_type == "IndexError":
                        raise IndexError(error_msg)
                    elif error_type == "PermissionError":
                        raise PermissionError(error_msg)
                except ValueError:
                    pass

            response.raise_for_status()
            return response.json()

    def replace_frame(self, frame_id: int, data: dict):
        """
        Replaces an existing logical frame with new data.
        This is a locked, non-destructive operation that appends the new data
        and updates the logical-to-physical mapping.
        """
        self._perform_locked_upload("replace", data, frame_id=frame_id)

    def insert_frame(self, index: int, data: dict):
        """
        Inserts a frame at the specified logical position.
        All frames at position >= index will be shifted to the right.

        Args:
            index: Logical position where to insert the frame (0-based)
            data: Dictionary containing data for the frame
        """
        self._perform_locked_upload("insert", data, insert_position=index)

    # MutableSequence interface implementation
    def __len__(self) -> int:
        """Return the number of frames."""
        return self._len

    @t.overload
    def __getitem__(self, index: int) -> dict: ...
    @t.overload
    def __getitem__(self, index: slice) -> list[dict]: ...
    @t.overload
    def __getitem__(self, index: list[int]) -> list[dict]: ...
    @t.overload
    def __getitem__(self, index: np.ndarray) -> dict | list[dict]: ...
    def __getitem__(
        self, index: int | slice | list[int] | np.ndarray
    ) -> dict | list[dict]:
        """Get frame(s) by index or slice, utilizing a local cache.
        
        For slices and lists, this method checks for partial hits in the cache
        and only fetches the missing frames from the server.
        """
        # Handle numpy arrays first
        if isinstance(index, np.ndarray):
            index = index.tolist() if index.ndim > 0 else int(index.item())

        length = len(self)

        if isinstance(index, int):
            # Caching logic for single integer index
            normalized_index = index if index >= 0 else length + index
            if not (0 <= normalized_index < length):
                raise IndexError("Index out of range")
            
            if normalized_index in self._cache:
                return self._cache[normalized_index]
            
            frame_data = self.get_frame(normalized_index)
            self._cache[normalized_index] = frame_data
            return frame_data

        elif isinstance(index, slice):
            # Caching logic for slices with partial hits
            if index.step is not None:
                if not isinstance(index.step, int):
                    raise TypeError("Slice step must be an integer")
                if index.step == 0:
                    raise ValueError("Slice step cannot be zero")

            resolved_indices = list(range(*index.indices(length)))
            if not resolved_indices:
                return []

            results_dict = {}
            misses = []

            # 1. Check the cache for hits and identify misses
            for idx in resolved_indices:
                if idx in self._cache:
                    results_dict[idx] = self._cache[idx]
                else:
                    misses.append(idx)
            
            # 2. If there were misses, fetch only them from the server
            if misses:
                fetched_data = self.get_frames(misses)
                
                # 3. Populate cache and results with the new data
                for idx, frame in zip(misses, fetched_data):
                    self._cache[idx] = frame
                    results_dict[idx] = frame
            
            # 4. Reconstruct the final list in the correct order
            return [results_dict[idx] for idx in resolved_indices]
            
        elif isinstance(index, list):
            # Caching logic for list of indices with partial hits
            if not index:
                return []

            validated_indices = []
            for i in index:
                if not isinstance(i, (int, np.integer)):
                    raise TypeError(f"List indices must be integers, not {type(i).__name__}")
                normalized_i = i if i >= 0 else length + i
                if not (0 <= normalized_i < length):
                    raise IndexError("list index out of range")
                validated_indices.append(normalized_i)
            
            results_dict = {}
            misses = []

            # 1. Check the cache for hits and identify unique misses
            for idx in validated_indices:
                if idx in self._cache:
                    results_dict[idx] = self._cache[idx]
                else:
                    # Avoid requesting the same missing index multiple times in one call
                    if idx not in misses:
                        misses.append(idx)
            
            # 2. If there were misses, fetch them from the server
            if misses:
                fetched_data = self.get_frames(misses)
                
                # 3. Populate cache and results with the new data
                for idx, frame in zip(misses, fetched_data):
                    self._cache[idx] = frame
                    results_dict[idx] = frame

            # 4. Reconstruct the final list, preserving original order and duplicates
            return [results_dict[idx] for idx in validated_indices]
        else:
            raise TypeError(
                f"Index must be int, slice, or list, not {type(index).__name__}"
            )


    def __setitem__(self, index, value):
        """Replace frame(s) at index, slice, or list of indices."""
        # Handle numpy arrays
        if isinstance(index, np.ndarray):
            if index.ndim == 0:
                # 0-d array (scalar)
                index = int(index.item())
            else:
                # Multi-dimensional array - convert to list
                index = index.tolist()

        if isinstance(index, slice):
            self._setitem_slice(index, value)
        elif isinstance(index, list):
            self._setitem_list(index, value)
        elif isinstance(index, int):
            # Single index assignment
            if index < 0:
                index += len(self)
            self.replace_frame(index, value)
        else:
            raise TypeError(
                f"Index must be int, slice, or list, not {type(index).__name__}"
            )

    def _setitem_slice(self, slice_obj: slice, values):
        """Handle slice assignment like Python lists."""
        current_len = len(self)

        # Convert values to list if it's iterable
        if hasattr(values, "__iter__") and not isinstance(values, dict):
            values = list(values)
        else:
            # Single value for slice assignment
            values = [values]

        # Calculate slice indices
        start, stop, step = slice_obj.indices(current_len)

        if step == 1:
            # Simple slice assignment: replace range with new values
            self._simple_slice_assignment(start, stop, values)
        else:
            # Extended slice assignment: must have same length
            self._extended_slice_assignment(start, stop, step, values)

    def _simple_slice_assignment(self, start: int, stop: int, values: list):
        """Handle simple slice assignment like data[2:5] = [a, b, c] using bulk endpoint."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            # Use the new bulk endpoint for efficient slice replacement
            serialized_data = [encode_data(value) for value in values]
            packed_data = msgpack.packb(serialized_data)

            bulk_url = f"{self.url}/api/rooms/{self.room}/frames/bulk"
            params = {"start": start, "stop": stop}

            response = requests.patch(
                bulk_url, data=packed_data, params=params, timeout=30
            )

            # Check for errors
            if response.status_code == 404:
                try:
                    error_data = response.json()
                    error_type = error_data.get("type", "")
                    error_msg = error_data.get("error", response.text)

                    if error_type == "IndexError":
                        raise IndexError(error_msg)
                except ValueError:
                    pass

            response.raise_for_status()

    def _setitem_list(self, indices: list, values):
        """Handle list index assignment like data[[1,2,3]] = [a, b, c] using bulk endpoint."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        # Convert values to list if it's iterable
        if hasattr(values, "__iter__") and not isinstance(values, dict):
            values = list(values)
        else:
            # Single value for all indices
            values = [values] * len(indices)

        # Validate list length
        if len(values) != len(indices):
            raise ValueError(
                f"attempt to assign sequence of size {len(values)} to list of size {len(indices)}"
            )

        # Validate and convert indices
        length = len(self)
        validated_indices = []
        for i in indices:
            if not isinstance(i, (int, np.integer)):
                raise TypeError(
                    f"List indices must be integers, not {type(i).__name__}"
                )
            # Convert negative indices
            if i < 0:
                i += length
            # Check bounds
            if i < 0 or i >= length:
                raise IndexError("list index out of range")
            validated_indices.append(int(i))

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            # Use the new bulk endpoint for efficient list replacement
            serialized_data = [encode_data(value) for value in values]
            packed_data = msgpack.packb(serialized_data)

            bulk_url = f"{self.url}/api/rooms/{self.room}/frames/bulk"
            params = {"indices": ",".join(str(i) for i in validated_indices)}

            response = requests.patch(
                bulk_url, data=packed_data, params=params, timeout=30
            )

            # Check for errors
            if response.status_code == 404:
                try:
                    error_data = response.json()
                    error_type = error_data.get("type", "")
                    error_msg = error_data.get("error", response.text)

                    if error_type == "IndexError":
                        raise IndexError(error_msg)
                except ValueError:
                    pass

            response.raise_for_status()

    def _extended_slice_assignment(
        self, start: int, stop: int, step: int, values: list
    ):
        """Handle extended slice assignment like data[::2] = [a, b, c] using bulk endpoint."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        # Calculate the indices that would be affected
        indices = list(range(start, stop, step))

        if len(values) != len(indices):
            raise ValueError(
                f"attempt to assign sequence of size {len(values)} to extended slice of size {len(indices)}"
            )

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            # Use the new bulk endpoint for efficient extended slice replacement
            serialized_data = [encode_data(value) for value in values]
            packed_data = msgpack.packb(serialized_data)

            bulk_url = f"{self.url}/api/rooms/{self.room}/frames/bulk"
            params = {"indices": ",".join(str(i) for i in indices)}

            response = requests.patch(
                bulk_url, data=packed_data, params=params, timeout=30
            )

            # Check for errors
            if response.status_code == 404:
                try:
                    error_data = response.json()
                    error_type = error_data.get("type", "")
                    error_msg = error_data.get("error", response.text)

                    if error_type == "IndexError":
                        raise IndexError(error_msg)
                except ValueError:
                    pass

            response.raise_for_status()

    def __delitem__(self, index: int | slice | list[int] | np.ndarray):
        """Delete frame(s) at index, slice, or list of indices."""
        # Handle numpy arrays
        if isinstance(index, np.ndarray):
            if index.ndim == 0:
                # 0-d array (scalar)
                index = int(index.item())
            else:
                # Multi-dimensional array - convert to list
                index = index.tolist()

        length = len(self)

        if isinstance(index, slice):
            # Validate slice step
            if index.step is not None:
                if not isinstance(index.step, int):
                    raise TypeError("Slice step must be an integer")
                if index.step == 0:
                    raise ValueError("Slice step cannot be zero")
            self.delete_frame(index)
        elif isinstance(index, list):
            # Validate list elements and convert negative indices
            validated_indices = []
            for i in index:
                if not isinstance(i, (int, np.integer)):
                    raise TypeError(
                        f"List indices must be integers, not {type(i).__name__}"
                    )
                # Convert negative indices
                if i < 0:
                    i += length
                # Check bounds
                if i < 0 or i >= length:
                    raise IndexError(f"list index out of range")
                validated_indices.append(int(i))
            self.delete_frame(validated_indices)
        elif isinstance(index, int):
            # Single index - convert negative and check bounds
            if index < 0:
                index += length
            if index < 0 or index >= length:
                raise IndexError(f"list index out of range")
            self.delete_frame(index)
        else:
            raise TypeError(
                f"Index must be int, slice, or list, not {type(index).__name__}"
            )

    def insert(self, index: int, value: dict):
        """Insert frame at index."""
        # Handle negative indices and clamp to valid range
        if index < 0:
            index = max(0, len(self) + index + 1)
        elif index > len(self):
            index = len(self)

        self.insert_frame(index, value)

    def append(self, value: dict):
        """Append a frame."""
        self.append_frame(value)

    def extend(self, values: list[dict]):
        """Extend with multiple frames."""
        if hasattr(values, "__iter__"):
            values = list(values)
        self.extend_frames(values)

    @property
    def settings(self) -> RoomConfig:
        def callback_fn(data, extension: str):
            # print("Settings updated from ZnDraw")
            category = "settings"
            response = requests.post(
                f"{self.url}/api/rooms/{self.room}/extensions/{category}/{extension}/submit?userId={self.user}",
                json=data,
            )
            response.raise_for_status()

        for key in settings:
            if key not in self._settings:
                response = requests.get(
                    f"{self.url}/api/rooms/{self.room}/extensions/settings/{key}/data?userId={self.user}"
                )
                data = response.json()
                if data["data"] is None:
                    self._settings[key] = settings[key]()
                else:
                    self._settings[key] = settings[key](**data["data"])
                self._settings[key].callback = functools.partial(
                    callback_fn, extension=key
                )

        config = RoomConfig(**self._settings)
        # TODO: do not allow changing config.<setting> directly, only sub-fields
        return config

    def register_extension(
        self,
        extension: t.Type[Extension],
        public: bool = False,
        run_kwargs: dict | None = None,
    ):
        # A WARNING ABOUT RUN_KWARGS!
        # If multiple workers are registering the same extension, work load will be distributed among them.
        # We do check that the extension schema are the same, but run_kwargs are not validated and can
        # lead to unreproducible behavior if they are different among workers.
        if not hasattr(extension, "category"):
            raise ValueError(
                "Extension must have a 'category' attribute and inherit from Extension base class."
            )
        name = extension.__name__
        if name in self._extensions:
            raise ValueError(f"Extension '{name}' is already registered.")
        if extension.category not in (cat.value for cat in ExtensionType):
            raise ValueError(
                f"Extension category '{extension.category}' is not valid. Must be one of {[cat.value for cat in ExtensionType]}."
            )
        self._extensions[name] = {
            "public": public,
            "run_kwargs": run_kwargs,
            "extension": extension,
        }
        print(f"Registered extension '{name}' of category '{extension.category}'.")

        schema = extension.model_json_schema()
        if public:
            raise NotImplementedError("Public extensions are not supported yet.")
        else:
            # /api/rooms/<string:room_id>/extensions/register
            response = requests.post(
                f"{self.url}/api/rooms/{self.room}/extensions/register",
                json={
                    "name": name,
                    "category": extension.category,
                    "schema": schema,
                    "clientId": self.sid,
                },
            )
            if response.status_code != 200:
                raise RuntimeError(
                    f"Failed to register extension '{name}': {response.status_code} {response.text}"
                )
            print(f"Extension '{name}' registered with room '{self.room}'.")

    def run(self, extension: Extension) -> None | dict:
        print(f"Running extension: {extension}")
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/extensions/{extension.category.value}/{extension.__class__.__name__}/submit",
            json={"data": extension.model_dump(), "userId": self.user},
        )

        response.raise_for_status()
        return response.json()

    def log(self, message: str) -> dict:
        """
        Send a chat message to the room.

        Args:
            message: Chat message content (markdown supported)

        Returns:
            dict: Response with success status and message data
        """
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        response = self.sio.call("chat:message:create", {"content": message}, timeout=5)
        return response

    def edit_message(self, message_id: str, new_content: str) -> dict:
        """
        Edit a previously sent message.

        Args:
            message_id: ID of the message to edit
            new_content: New message content

        Returns:
            dict: Response with success status and updated message data
        """
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        response = self.sio.call(
            "chat:message:edit",
            {"messageId": message_id, "content": new_content},
            timeout=5,
        )
        return response

    def get_messages(
        self, limit: int = 30, before: int = None, after: int = None
    ) -> dict:
        """
        Fetch chat message history.

        Args:
            limit: Number of messages to fetch (1-100)
            before: Get messages before this timestamp
            after: Get messages after this timestamp

        Returns:
            dict: Response with messages and metadata
        """
        params = {"limit": limit}
        if before:
            params["before"] = before
        if after:
            params["after"] = after

        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/chat/messages", params=params
        )
        response.raise_for_status()
        return response.json()
