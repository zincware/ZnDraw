import dataclasses
from abc import abstractmethod
from collections.abc import MutableSequence

import numpy as np
import socketio


@dataclasses.dataclass
class ZnDrawBase(MutableSequence):
    token: str
    url: str

    socket: socketio.Client = dataclasses.field(default_factory=socketio.Client)

    def __post_init__(self):
        self.url = self.url.replace("http", "ws")
        print(f"Connecting to {self.url}")
        self.socket.connect(self.url, wait_timeout=1)
        self.socket.emit("join", str(self.token))

    def reconnect(self):
        self.socket.connect(self.url)
        self.socket.emit("join", str(self.token))

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
    def wrap_and_check_index(index: int | slice | list[int], length: int) -> list[int]:
        is_slice = isinstance(index, slice)
        if is_slice:
            index = list(range(*index.indices(length)))
        index = [index] if isinstance(index, int) else index
        index = [i if i >= 0 else length + i for i in index]
        # check if index is out of range
        for i in index:
            if i >= length:
                raise IndexError(f"Index {i} out of range for length {length}")
        return index
