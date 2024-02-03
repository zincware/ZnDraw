import dataclasses
from abc import abstractmethod
import collections
import socketio
import typing
import numpy as np

@dataclasses.dataclass
class ZnDrawBase(collections.abc.MutableSequence):
    socket: socketio.Client
    token: str

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

    


    