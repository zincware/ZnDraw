import dataclasses
import logging
from abc import abstractmethod
from collections.abc import MutableSequence

import numpy as np
import socketio
import splines

from zndraw.data import CeleryTaskData

log = logging.getLogger(__name__)


@dataclasses.dataclass
class ZnDrawBase(MutableSequence):
    url: str
    token: str

    socket: socketio.Client = dataclasses.field(default_factory=socketio.Client)

    def __post_init__(self):
        self.url = self.url.replace("http", "ws")
        log.critical(f"Connecting to {self.url}")
        self.socket.connect(self.url, wait_timeout=1)
        self.socket.emit("join", {"token": str(self.token), "auth_token": None})

    def reconnect(self):
        self.socket.connect(self.url)
        self.socket.emit("join", {"token": str(self.token), "auth_token": None})

    @abstractmethod
    def log(self, message: str):
        pass

    @property
    @abstractmethod
    def bookmarks(self) -> dict[int, str]:
        pass

    @bookmarks.setter
    @abstractmethod
    def bookmarks(self, value: dict[int, str]):
        pass

    @property
    @abstractmethod
    def step(self) -> int:
        pass

    @step.setter
    @abstractmethod
    def step(self, value: int):
        pass

    @property
    @abstractmethod
    def selection(self) -> list[int]:
        pass

    @selection.setter
    @abstractmethod
    def selection(self, value: list[int]):
        pass

    @property
    @abstractmethod
    def points(self) -> np.ndarray:
        pass

    @points.setter
    @abstractmethod
    def points(self, value: np.ndarray):
        pass

    @property
    @abstractmethod
    def segments(self) -> np.ndarray:
        pass

    @property
    def figure(self):
        raise NotImplementedError("Gathering figure from webclient not implemented yet")

    @figure.setter
    def figure(self, fig: str):
        data = {"figure": fig, "token": self.token}
        self.socket.emit("analysis:figure", data)

    @property
    def atoms(self):
        return self[self.step]

    @staticmethod
    def calculate_segments(points: np.ndarray) -> np.ndarray:
        if points.shape[0] <= 1:
            return points
        t = np.linspace(0, len(points) - 1, len(points) * 50)
        return splines.CatmullRom(points).evaluate(t)

    @property
    def camera(self):
        raise NotImplementedError("Getting camera from webclient not implemented yet")

    @camera.setter
    def camera(self, camera: dict):
        """Set the camera position and orientation

        camera: dict
            A dictionary with the following
            - position: list[float]
                The position of the camera
            - target: list[float]
                The target of the camera
        """
        if set(camera) != {"position", "target"}:
            raise ValueError("camera must have keys 'position' and 'target'")
        msg = CeleryTaskData(
            target=str(self.token),
            event="camera:update",
            data=camera,
        )
        self.socket.emit("celery:task:emit", msg.to_dict())

    def get_screenshot(self, filename):
        
        def _receive_screenshot(data):
            data = data.replace("data:image/png;base64,", "")
            print(data)
            import base64
            with open(filename, "wb") as fh:
                fh.write(base64.decodebytes(bytes(data, "utf-8")))
            self.socket.on("screenshot", lambda x: None)

        self.socket.on("screenshot", _receive_screenshot)
        self.socket.emit("screenshot")
