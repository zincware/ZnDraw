import dataclasses
import socketio
import warnings
import numpy as np
import requests
import typing as t
from tqdm import tqdm
import msgpack
import logging

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
            warnings.warn(f"Failed to release lock for target '{self.target}'. It may have expired.")

@dataclasses.dataclass
class Client:
    """A client for interacting with the ZnDraw server."""
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

    def get_frame(self, frame_id: int) -> dict[str, np.ndarray]:
        """Fetches a single frame's data from the server."""
        full_url = f"{self.url}/frame/{self.room}/{frame_id}"
        response = requests.get(full_url, timeout=10)
        response.raise_for_status()

        # Unpack msgpack data
        serialized_data = msgpack.unpackb(response.content, strict_map_key=False)

        # Reconstruct numpy arrays from bytes, shape, and dtype
        result = {}
        for key, array_info in serialized_data.items():
            data_bytes = array_info['data']
            shape = tuple(array_info['shape'])
            dtype = array_info['dtype']

            result[key] = np.frombuffer(data_bytes, dtype=dtype).reshape(shape)

        return result

    def append_frame(self, data: dict[str, np.ndarray]):
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        lock = SocketIOLock(self.sio, target="trajectory:meta")

        with lock:
            response = self.sio.call("upload:prepare", {})
            
            if not response or not response.get("success"):
                raise RuntimeError(f"Failed to prepare for upload: {response.get('error')}")
            
            token = response["token"]            
            try:
                # Pack the data using msgpack with bytes and shape info
                serialized_data = {}
                for key, array in data.items():
                    serialized_data[key] = {
                        'data': array.tobytes(),
                        'shape': array.shape,
                        'dtype': str(array.dtype)
                    }

                packed_data = msgpack.packb(serialized_data)

                upload_url = f"{self.url}/rooms/{self.room}/frames"
                headers = {
                    "Authorization": f"Bearer {token}",
                    "Content-Type": "application/octet-stream"
                }

                http_response = requests.post(
                    upload_url,
                    data=packed_data,
                    headers=headers,
                    timeout=30
                )
                http_response.raise_for_status()

                result = http_response.json()
                if not result.get("success"):
                    raise RuntimeError(f"Server reported failure: {result.get('error')}")
                
            except requests.exceptions.RequestException as e:
                # Wrap the HTTP error in a RuntimeError
                raise RuntimeError(f"Error uploading frame data: {e}") from e

    def len_frames(self) -> int:
        """Returns the number of frames in the current room."""
        if not self.sio.connected:
            raise RuntimeError("Client is not connected. Please call .connect() first.")

        response = self.sio.call("frames:count", {})

        if not response or not response.get("success"):
            raise RuntimeError(f"Failed to get frame count: {response.get('error') if response else 'No response'}")

        return response["count"]

if __name__ == '__main__':
    client = Client()
    client.connect()
    for idx in tqdm(range(50)):
        client.append_frame({"index": np.array([idx])})
    for idx in range(50):
        frame = client.get_frame(idx)
        print(f"Frame {idx} keys: {list(frame.keys())}, index: {frame['index']}")
    
    print(f"Total frames: {client.len_frames()}")
