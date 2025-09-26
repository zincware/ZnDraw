import dataclasses
import json
import logging
import warnings
from collections.abc import MutableSequence

import msgpack
import numpy as np
import requests
import socketio
from tqdm import tqdm

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


@dataclasses.dataclass
class Client(MutableSequence):
    """A client for interacting with the ZnDraw server. Implements MutableSequence for frame operations."""

    room: str = "default"
    url: str = "http://localhost:5000"

    def __post_init__(self):
        self.sio = socketio.Client()
        self.sio.on("connect", self._on_connect)

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
        log.debug(f"Connected to {self.url} with session ID {self.sio.sid}")
        self.sio.emit("join_room", {"room": self.room})
        log.debug(f"Joined room: '{self.room}'")

    def get_frame(self, frame_id: int, keys: list[str] | None = None) -> dict:
        """Fetches a single frame's data from the server."""
        # Handle negative indices
        if frame_id < 0:
            frame_id = self.len_frames() + frame_id
        return self.get_frames([frame_id], keys=keys)[0]

    def _serialize_frame_data(self, data: dict) -> dict:
        """Convert nested dict with numpy arrays to server-compatible format."""
        flattened = self._flatten_data(data)
        # Convert to server format (data/shape/dtype for arrays)
        serialized = {}
        non_array_data = {}

        for key, value in flattened.items():
            if isinstance(value, np.ndarray):
                serialized[key] = {
                    "data": value.tobytes(),
                    "shape": value.shape,
                    "dtype": str(value.dtype),
                }
            else:
                # Collect non-array values for metadata
                non_array_data[key] = value

        # If we have non-array data, store it as a JSON-encoded numpy array
        if non_array_data:
            json_str = json.dumps(non_array_data)
            json_array = np.array([json_str], dtype="U")
            serialized["_metadata"] = {
                "data": json_array.tobytes(),
                "shape": json_array.shape,
                "dtype": str(json_array.dtype),
            }

        return serialized

    def _flatten_data(self, data, prefix=""):
        """Flatten nested dictionaries using dot notation for keys."""
        flattened = {}
        for key, value in data.items():
            full_key = f"{prefix}.{key}" if prefix else key
            if isinstance(value, dict):
                # Recursively flatten nested dicts
                flattened.update(self._flatten_data(value, full_key))
            else:
                flattened[full_key] = value
        return flattened

    def _unflatten_data(self, flattened):
        """Reconstruct nested dictionary from flattened dot notation."""
        result = {}
        for key, value in flattened.items():
            parts = key.split(".")
            current = result
            for part in parts[:-1]:
                if part not in current:
                    current[part] = {}
                current = current[part]
            current[parts[-1]] = value
        return result

    def _deserialize_frame_data(self, serialized_data):
        """Convert server format back to nested dict with numpy arrays."""
        # First convert arrays back from server format
        converted = {}
        metadata = {}

        for key, array_info in serialized_data.items():
            if key == "_metadata":
                # Handle metadata specially
                data_bytes = array_info["data"]
                shape = tuple(array_info["shape"])
                dtype = array_info["dtype"]

                # Reconstruct the JSON string array
                json_array = np.frombuffer(data_bytes, dtype=dtype).reshape(shape)
                json_str = str(json_array[0]) if len(json_array) > 0 else "{}"

                try:
                    metadata = json.loads(json_str)
                except (json.JSONDecodeError, TypeError):
                    metadata = {}
            else:
                # Regular array
                data_bytes = array_info["data"]
                shape = tuple(array_info["shape"])
                dtype = array_info["dtype"]
                converted[key] = np.frombuffer(data_bytes, dtype=dtype).reshape(shape)

        # Add metadata back
        converted.update(metadata)

        # Unflatten the structure
        return self._unflatten_data(converted)

    def _serialize_frames_data(self, frames: list[dict]) -> list[dict]:
        """Convert list of frame data to msgpack-compatible format."""
        return [self._serialize_frame_data(frame) for frame in frames]

    def _prepare_upload_token(self, action: str, **kwargs) -> str:
        """Prepare upload token for frame operations."""
        request_data = {"action": action}
        request_data.update(kwargs)

        response = self.sio.call("upload:prepare", request_data)
        if not response or not response.get("success"):
            raise RuntimeError(
                f"Failed to prepare for upload: {response.get('error') if response else 'No response'}"
            )
        return response["token"]

    def _upload_frame_data(self, token: str, serialized_data) -> dict:
        """Upload frame data using the provided token."""
        packed_data = msgpack.packb(serialized_data)

        upload_url = f"{self.url}/rooms/{self.room}/frames"
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
                    self._serialize_frame_data(data)
                    if isinstance(data, dict)
                    else self._serialize_frames_data(data)
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

        full_url = f"{self.url}/frames/{self.room}"
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

        # Unpack msgpack data
        serialized_frames = msgpack.unpackb(response.content, strict_map_key=False)

        # Deserialize frames
        return [self._deserialize_frame_data(frame) for frame in serialized_frames]

    def len_frames(self) -> int:
        """Returns the number of frames in the current room."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        response = self.sio.call("frames:count", {})

        if not response or not response.get("success"):
            raise RuntimeError(
                f"Failed to get frame count: {response.get('error') if response else 'No response'}"
            )

        return response["count"]

    def delete_frame(self, frame_id: int):
        """Deletes a frame from the current room."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            response = self.sio.call("frame:delete", {"frame_id": frame_id})

            if not response or not response.get("success"):
                raise RuntimeError(
                    f"Failed to delete frame: {response.get('error') if response else 'No response'}"
                )

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
        return self.len_frames()

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
                serialized_data = self._serialize_frame_data(value)
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
                    serialized_data = self._serialize_frame_data(value)
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
