import dataclasses
import logging
import typing as t
from abc import abstractmethod
from collections.abc import MutableSequence

import ase
import numpy as np
import splines
import znjson
import znsocket
from flask import current_app, session
from pydantic import BaseModel, Field, create_model
from redis import Redis

from zndraw.utils import ASEConverter

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

    @staticmethod
    def get_atoms() -> ase.Atoms:
        """Get the ase atoms object at the current position in the room"""
        room = session["token"]
        r: Redis = current_app.extensions["redis"]
        step = r.get(f"room:{room}:step")
        key = (
            f"room:{room}:frames"
            if r.exists(f"room:{room}:frames")
            else "room:default:frames"
        )
        lst = znsocket.List(r, key)
        try:
            frame_json = lst[int(step)]
            return znjson.loads(
                frame_json, cls=znjson.ZnDecoder.from_converters([ASEConverter])
            )
        except TypeError:
            # step is None
            return ase.Atoms()
        except IndexError:
            return ase.Atoms()


class MethodsCollection(BaseModel):
    """Base class for collections of methods for modification, analysis, etc."""

    method: t.Type[Extension] = Field(..., description="Select a method.")

    def run(self, vis, **kwargs) -> None:
        self.method.run(vis, **kwargs)

    @classmethod
    def updated_schema(
        cls, extensions: t.Optional[t.List[t.Type[Extension]]] = None
    ) -> dict:
        methods = cls.__annotations__["method"]
        if extensions is not None and len(extensions) > 0:
            extensions_types = t.Union[tuple(extensions)]
            extended_methods = t.Union[extensions_types, methods]
        else:
            extended_methods = methods

        # get the description of the cls.method field
        method_field = cls.model_fields["method"]

        field_kwargs = {
            "description": method_field.description,
            "discriminator": "discriminator",
        }
        field_kwargs["default"] = method_field.default
        field_kwargs["default_factory"] = method_field.default_factory

        extended_cls = create_model(
            cls.__name__,
            __base__=cls,
            method=(
                extended_methods,
                Field(**field_kwargs),
            ),
        )
        schema = extended_cls.model_json_schema()
        # TODO: iterate through all fields that have the
        for prop in [x.__name__ for x in t.get_args(extended_methods)]:
            schema["$defs"][prop]["properties"]["discriminator"]["options"] = {
                "hidden": True
            }

        return schema


@dataclasses.dataclass  # TODO: move to a separate file, so it can be imported in other files
class FileIO:
    name: str | None = None
    start: int = 0
    stop: int | None = None
    step: int = 1
    remote: str | None = None
    rev: str | None = None

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
    def figure(self):
        pass

    @figure.setter
    @abstractmethod
    def figure(self, fig: str):
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
