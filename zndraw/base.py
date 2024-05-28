import dataclasses
import logging
import typing as t
from abc import abstractmethod
from collections.abc import MutableSequence

import ase
import numpy as np
import splines
import znframe
from flask import current_app, session
from pydantic import BaseModel, Field, create_model
from redis import Redis


log = logging.getLogger(__name__)


class RedisList(MutableSequence):
    def __init__(self, redis: Redis, key: str):
        self.redis = redis
        self.key = key
    
    def __len__(self) -> int:
        return int(self.redis.llen(self.key))
    
    def __getitem__(self, index: int|list|slice):
        single_item = isinstance(index, int)
        if single_item:
            index = [index]
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))

        items = []
        for i in index:
            item = self.redis.lindex(self.key, i)
            if item is None:
                raise IndexError("list index out of range")
            items.append(item)
        return items[0] if single_item else items
    
    def __setitem__(self, index: int|list|slice, value: str|list[str]):
        single_item = isinstance(index, int)
        if single_item:
            index = [index]
            assert isinstance(value, str), "single index requires single value"
            value = [value]
        
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        
        index = [int(i) for i in index]
        
        if len(index) != len(value):
            raise ValueError(f"attempt to assign sequence of size {len(value)} to extended slice of size {len(index)}")
        
        for i, v in zip(index, value):
            # TODO: this prohibits appending to the list, right?
            if i >= self.__len__() or i < -self.__len__():
                raise IndexError("list index out of range")
            self.redis.lset(self.key, i, v)
    
    def __delitem__(self, index):
        current_list = self.redis.lrange(self.key, 0, -1)
        if index >= len(current_list) or index < -len(current_list):
            raise IndexError("list index out of range")
        self.redis.lset(self.key, index, "__DELETED__")
        self.redis.lrem(self.key, 1, "__DELETED__")
    
    def insert(self, index, value):
        if index >= self.__len__():
            self.redis.rpush(self.key, value)
        elif index == 0:
            self.redis.lpush(self.key, value)
        else:
            pivot = self.redis.lindex(self.key, index)
            self.redis.linsert(self.key, 'BEFORE', pivot, value)
    
    def __iter__(self):
        return (item for item in self.redis.lrange(self.key, 0, -1))

    def __repr__(self):
        return f"RedisList({self.redis.lrange(self.key, 0, -1)})"

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
        r: Redis = current_app.config["redis"]
        step = r.get(f"room:{room}:step")
        frame_json = r.lindex(f"room:{room}:frames", int(step))
        if frame_json is None:
            # probably default room
            frame_json = r.hget("room:default:frames", f"{step}")
        frame = znframe.Frame.from_json(frame_json)
        return frame.to_atoms()


class MethodsCollection(BaseModel):
    """Base class for collections of methods for modification, analysis, etc."""

    method: t.Type[Extension] = Field(..., description="Select a method.")

    def run(self, vis, **kwargs) -> None:
        self.method.run(vis, **kwargs)

    @classmethod
    def updated_schema(cls, extensions: list[t.Type[Extension]] | None = None) -> dict:
        methods = cls.__annotations__["method"]
        if extensions is not None:
            extensions_types = t.Union[tuple(extensions)]
            extended_methods = t.Union[methods, extensions_types]
        else:
            extended_methods = methods

        # get the description of the cls.method field
        method_description = cls.model_fields["method"].description

        extended_cls = create_model(
            cls.__name__,
            __base__=cls,
            method=(
                extended_methods,
                Field(
                    ..., description=method_description, discriminator="discriminator"
                ),
            ),
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