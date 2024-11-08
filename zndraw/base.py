import dataclasses
import logging
from abc import abstractmethod
from collections.abc import MutableSequence

import numpy as np
import splines
from pydantic import BaseModel

log = logging.getLogger(__name__)


class Extension(BaseModel):
    def run(self, vis, **kwargs) -> None:
        raise NotImplementedError("run method must be implemented in subclass")


@dataclasses.dataclass  # TODO: move to a separate file, so it can be imported in other files
class FileIO:
    name: str | None = None
    start: int = 0
    stop: int | None = None
    step: int = 1
    remote: str | None = None
    rev: str | None = None
    convert_nan: bool = False

    def to_dict(self):
        return dataclasses.asdict(self)


class ZnDrawBase(MutableSequence):
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
    def segments(self) -> np.ndarray:
        points = self.points
        if points.shape[0] <= 1:
            return points
        t = np.linspace(0, len(points) - 1, len(points) * 50)
        return splines.CatmullRom(points).evaluate(t)

    @property
    @abstractmethod
    def figures(self):
        pass

    @figures.setter
    @abstractmethod
    def figures(self, fig: str):
        pass

    @property
    def atoms(self):
        return self[self.step]

    @property
    @abstractmethod
    def camera(self):
        pass

    @camera.setter
    @abstractmethod
    def camera(self, camera: dict):
        pass

    @property
    @abstractmethod
    def locked(self) -> bool:
        pass

    @locked.setter
    @abstractmethod
    def locked(self, value: bool) -> None:
        pass
