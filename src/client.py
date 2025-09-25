import dataclasses
import socketio
import warnings
import numpy as np
import requests
from tqdm import tqdm

@dataclasses.dataclass
class Client:
    room: str = "default"
    url: str = "http://localhost:5000"


    def __post_init__(self):
        self.sio = socketio.Client()
        self.lock = SocketIOLock(self.sio, target="trajectory:meta")
        self.sio.on("connect", self.on_connect)
        self.sio.connect(self.url)

    def on_connect(self):
        print(f"Connected to {self.url}")
        self.sio.emit("join_room", {"room": self.room})

    def get_frame(self, frame_id: int) -> np.ndarray:
        response = requests.get(f"{self.url}/frame/{self.room}/{frame_id}", timeout=10)
        response.raise_for_status()
        positions = np.frombuffer(response.content, dtype=np.float32)
        return positions

    def append_frame(self, data: np.ndarray[tuple[int, int], np.float32]):
        with self.lock:
            response = self.sio.call("upload:prepare", {})
            if not response.get("success"):
                print("Failed to prepare for upload")
                return
            token = response["token"]
            try:
                headers = {
                    "Authorization": f"Bearer {token}",
                    "Content-Type": "application/octet-stream"
                }
                response = requests.post(
                    f"{self.url}/upload",
                    data=data.tobytes(),
                    headers=headers,
                    timeout=30
                )
                response.raise_for_status() # Raise an exception for 4xx/5xx errors
                
                result = response.json()
                if not result.get("success"):
                    print(f"Server reported failure: {result.get('error')}")
                    return

            except requests.exceptions.RequestException as e:
                print(f"Error uploading frame data: {e}")




@dataclasses.dataclass
class SocketIOLock:
    sio: socketio.Client
    target: str

    def acquire(self, blocking: bool = True, timeout: float = 60):
        """Acquire a lock."""
        return self.sio.call("lock:acquire", {"target": self.target}, timeout=int(timeout))

    def release(self):
        """Release a lock."""
        return self.sio.call("lock:release", {"target": self.target})

    def __enter__(self):
        if not self.acquire():
            raise RuntimeError("Failed to acquire lock")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.release():
            warnings.warn("Failed to release lock")

if __name__ == '__main__':
    import znh5md
    io = znh5md.IO("/Users/fzills/tools/zndraw-communication-testing/structures.h5")
    positions = np.array([frame.positions for frame in io[:]])
    client = Client()
    data = np.random.rand(100, 3).astype(np.float32)
    client.append_frame(positions)
    for frame in tqdm(io[:]):
        client.append_frame(frame.positions)
    # for _ in tqdm(range(1000)):
    #     data = np.random.rand(10_000, 3).astype(np.float32)
    #     client.append_frame(data)

    print(client.get_frame(0))
    client.sio.wait()



# import requests
# import numpy as np
# from tqdm import tqdm

# def get_frame(frame_id: int, server_url: str = "http://localhost:5000") -> np.ndarray:
#     """
#     Fetches binary frame data from the Flask server and returns it as a NumPy array.

#     Args:
#         frame_id: The integer ID of the frame to retrieve.
#         server_url: The base URL of the Flask server.

#     Returns:
#         A NumPy array of shape (N, 3) containing the atomic positions,
#         or None if an error occurs.
#     """
#     url = f"{server_url}/frame/{frame_id}"    
#     response = requests.get(url, timeout=10)
#     response.raise_for_status()
#     positions = np.frombuffer(response.content, dtype=np.float32)
#     positions_reshaped = positions.reshape(-1, 3)
#     return positions_reshaped


# # --- Example Usage ---
# if __name__ == '__main__':
#     frames = []
#     print("Round 1")
#     for idx in tqdm(range(50)):
#         frame_data = get_frame(idx)
#         frames.append(frame_data)
#     print("Round 2")
#     for idx in tqdm(range(50)):
#         frame_data = get_frame(idx)
#         frames.append(frame_data)
#     print("Done")
#     print(np.array(frames).shape)

# import socketio

# sio = socketio.Client()

# @sio.on("set_frame")
# def on_set_frame(data):
#     frame_id = data["frame"]
#     print(f"Switched to frame {frame_id}")

# sio.connect("http://localhost:5000")

# # Example: set frame from this client
# sio.emit("set_frame", {"frame": 10})
# sio.wait()  # Keep the client running to listen for events