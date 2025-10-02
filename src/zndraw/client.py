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
    _lock: SocketIOLock = dataclasses.field(init=False)

    def __post_init__(self):
        self.sio = socketio.Client()
        self.sio.on("connect", self._on_connect)
        self.sio.on("frame_update", self._on_frame_update)
        self.sio.on("selection:update", self._on_selection_update)
        self.sio.on("len_frames", self._on_len_frames_update)
        self.sio.on("invalidate", self._on_invalidate)
        self.sio.on("queue:update", self._on_queue_update)
        self.sio.on("frame_selection:update", self._on_frame_selection_update)

        self._lock = SocketIOLock(self.sio, target="trajectory:meta")

        if self.template is _TemplateValue:
            response = requests.post(f"{self.url}/api/room/{self.room}/join", json={})
        else:
            response = requests.post(
                f"{self.url}/api/room/{self.room}/join",
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
                    response = requests.post(
                        f"{self.url}/api/rooms/{self.room}/jobs/{data.get('jobId')}/complete",
                        json={"workerId": self.sid},
                    )
                    if response.status_code != 200:
                        log.error(
                            f"Failed to mark job {data.get('jobId')} as complete: {response.status_code} {response.text}"
                        )
                    # log as completed
                except Exception as e:
                    log.error(f"Error processing job {data.get('jobId')}: {e}")
                    response = requests.post(
                        f"{self.url}/api/rooms/{self.room}/jobs/{data.get('jobId')}/fail",
                        json={"workerId": self.sid, "error": str(e)},
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

    def get_frame(self, frame_id: int, keys: list[str] | None = None) -> dict:
        """Fetches a single frame's data from the server."""
        # Handle negative indices
        if frame_id < 0:
            frame_id = self.len_frames() + frame_id
        return self.get_frames([frame_id], keys=keys)[0]

    def _prepare_upload_token(self, action: str, **kwargs) -> str:
        """Prepare upload token for frame operations."""
        request_data = {"action": action}
        request_data.update(kwargs)

        response = self.sio.call("upload:prepare", request_data)
        if not response or not response.get("success"):
            error_msg = response.get("error") if response else "No response"
            error_type = response.get("error_type") if response else None

            # Raise the appropriate error type based on server response
            if error_type == "IndexError":
                raise IndexError(error_msg)
            raise RuntimeError(f"Failed to prepare for upload: {error_msg}")
        return response["token"]

    def _upload_frame_data(self, token: str, serialized_data) -> dict:
        """Upload frame data using the provided token."""
        packed_data = msgpack.packb(serialized_data)

        upload_url = f"{self.url}/api/rooms/{self.room}/frames"
        headers = {
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/octet-stream",
        }

        http_response = requests.post(
            upload_url, data=packed_data, headers=headers, timeout=30
        )
        http_response.raise_for_status()

        result = http_response.json()
        if not result.get("success"):
            raise RuntimeError(f"Server reported failure: {result.get('error')}")

        return result

    def _perform_locked_upload(self, action: str, data, **kwargs):
        """Perform a locked upload operation with common error handling."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            try:
                token = self._prepare_upload_token(action, **kwargs)
                serialized_data = (
                    encode_data(data)
                    if isinstance(data, dict)
                    else [encode_data(frame) for frame in data]
                )
                return self._upload_frame_data(token, serialized_data)
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
            # Direct list of indices
            payload = {"indices": indices_or_slice}
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

        # Add keys parameter if specified
        if keys is not None:
            payload["keys"] = keys

        full_url = f"{self.url}/api/frames/{self.room}"
        response = requests.post(full_url, json=payload, timeout=30)

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

    def delete_frame(self, frame_id: int):
        """Deletes a frame from the current room."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            response = self.sio.call("frame:delete", {"frame_id": frame_id})

            if not response or not response.get("success"):
                error_msg = response.get("error") if response else "No response"
                error_type = response.get("error_type") if response else None

                # Raise the appropriate error type based on server response
                if error_type == "IndexError":
                    raise IndexError(error_msg)
                raise RuntimeError(f"Failed to delete frame: {error_msg}")

            return response

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

    def __getitem__(self, index) -> dict | list[dict]:
        """Get frame(s) by index or slice."""
        if isinstance(index, slice):
            return self.get_frames(index)
        elif isinstance(index, int):
            if index < 0:
                index += len(self)
            return self.get_frame(index)
        else:
            raise TypeError("Index must be int or slice")

    def __setitem__(self, index, value):
        """Replace frame(s) at index or slice."""
        if isinstance(index, slice):
            self._setitem_slice(index, value)
        else:
            # Single index assignment
            if index < 0:
                index += len(self)
            self.replace_frame(index, value)

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
        """Handle simple slice assignment like data[2:5] = [a, b, c]."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        # Number of positions being replaced
        old_count = max(0, stop - start)

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            # TODO: this is very chatty and should be improved in the server!!

            # First, delete the old range if it exists
            for _ in range(old_count):
                if start < len(self):
                    response = self.sio.call("frame:delete", {"frame_id": start})
                    if not response or not response.get("success"):
                        raise RuntimeError(
                            f"Failed to delete frame: {response.get('error') if response else 'No response'}"
                        )

            # Then insert the new values at the start position
            for i, value in enumerate(values):
                # Use the direct upload mechanism
                token = self._prepare_upload_token("insert", insert_position=start + i)
                serialized_data = encode_data(value)
                self._upload_frame_data(token, serialized_data)

    def _extended_slice_assignment(
        self, start: int, stop: int, step: int, values: list
    ):
        """Handle extended slice assignment like data[::2] = [a, b, c]."""
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
            # Replace each position individually using direct calls
            for i, value in zip(indices, values):
                if i < len(self):
                    token = self._prepare_upload_token("replace", frame_id=i)
                    serialized_data = encode_data(value)
                    self._upload_frame_data(token, serialized_data)

    def __delitem__(self, index: int):
        """Delete frame at index."""
        if isinstance(index, slice):
            raise NotImplementedError("Slice deletion not supported")
        if index < 0:
            index += len(self)
        self.delete_frame(index)

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
                f"{self.url}/api/rooms/{self.room}/extensions/{category}/{extension}?userId={self.user}",
                json=data,
            )
            response.raise_for_status()

        for key in settings:
            if key not in self._settings:
                response = requests.get(
                    f"{self.url}/api/rooms/{self.room}/extension-data/settings/{key}?userId={self.user}"
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
        if not public:
            response = self.sio.call(
                "register:extension",
                {
                    "name": name,
                    "category": extension.category,
                    "schema": schema,
                    "public": False,
                },
            )
            print(f"Extension '{name}' registered with room '{self.room}'.")
        else:
            response = self.sio.call(
                "register:extension",
                {
                    "name": name,
                    "category": extension.category,
                    "schema": schema,
                    "public": True,
                },
            )
            print(f"Extension '{name}' registered as public.")

        if response.get("status") != "success":
            raise RuntimeError(f"Failed to register extension '{name}': {response}")


if __name__ == "__main__":
    client = Client()
    client.connect()
    for idx in tqdm(range(50)):
        client.append_frame({"index": np.array([idx])})
    for idx in range(client.len_frames()):
        frame = client.get_frame(idx)
        print(f"Frame {idx} keys: {list(frame.keys())}, index: {frame['index']}")
    client.delete_frame(0)
    client.delete_frame(0)
    print("After deletion:")
    for idx in range(client.len_frames()):
        frame = client.get_frame(idx)
        print(f"Frame {idx} keys: {list(frame.keys())}, index: {frame['index']}")
    client.replace_frame(5, {"index": np.array([999])})
    for idx in range(10):
        frame = client.get_frame(idx)
        print(f"Frame {idx} keys: {list(frame.keys())}, index: {frame['index']}")

    data = [{"index": np.array([1000 + i])} for i in range(50)]
    for _ in tqdm(range(1)):
        client.extend_frames(data)
    for _ in tqdm(range(1)):
        client.extend_frames(data)
    for entry in tqdm(data):
        client.append_frame(entry)

    all_data = client.get_frames(slice(None, None, None))
    print(f"Fetched {len(all_data)} frames via get_frames.")
    print(f"Total frames: {client.len_frames()}")
