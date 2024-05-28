import dataclasses
import logging
import typing as t
from abc import abstractmethod
from collections.abc import MutableSequence

import numpy as np
import socketio
import splines
from pydantic import BaseModel, Field, create_model

from zndraw.data import CeleryTaskData

log = logging.getLogger(__name__)


class Extension(BaseModel):
    # TODO: can I hide the discriminator field in the model json schema automatically here?
    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        # Automatically add the discriminator field
        cls.__annotations__["discriminator"] = t.Literal[cls.__name__]
        setattr(cls, "discriminator", cls.__name__)
    
    def run(self, vis, **kwargs) -> None:
        raise NotImplementedError("run method must be implemented in subclass")


class MethodsCollection(BaseModel):
    """Base class for collections of methods for modification, analysis, etc."""
    method: t.Type[Extension] = Field(
        ..., description="Select a method."
    )

    def run(self, vis, **kwargs) -> None:
        self.method.run(vis, **kwargs)

    @classmethod
    def updated_schema(cls, extensions: list[t.Type[Extension]]|None = None) -> dict:
        methods = cls.__annotations__["method"]
        if extensions is not None:
            extensions_types = t.Union[tuple(extensions)]
            extended_methods = t.Union[methods, extensions_types]
        else:
            extended_methods = methods

        # get the description of the cls.method field
        method_description = cls.model_fields['method'].description

        extended_cls = create_model(
            cls.__name__,
            __base__=cls,
            method=(extended_methods, Field(
                ..., description=method_description, discriminator="discriminator"
            ))
        )
        schema = extended_cls.model_json_schema()
        for prop in [x.__name__ for x in t.get_args(methods)]:
            schema["$defs"][prop]["properties"]["discriminator"]["options"] = {
                "hidden": True
            }
        return schema


@dataclasses.dataclass  # TODO: move to a separate file, so it can be imported in other files
class FileIO:
    name: str = None
    start: int = 0
    stop: int = None
    step: int = 1
    remote: str = None
    rev: str = None

    def to_dict(self):
        return dataclasses.asdict(self)


@dataclasses.dataclass
class ZnDrawBase(MutableSequence):
    url: str
    token: str
    auth_token: str = None

    socket: socketio.Client = dataclasses.field(default_factory=socketio.Client)

    def __post_init__(self):
        self.url = self.url.replace("http", "ws")
        self.socket.on("connect", self._on_connect)
        self.socket.connect(self.url, wait_timeout=1)

    def reconnect(self):
        self.socket.connect(self.url)

    def _on_connect(self):
        self.socket.emit(
            "join",
            {
                "token": str(self.token),
                "auth_token": self.auth_token,
            },
        )

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
