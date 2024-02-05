import logging
import typing as t

import ase
import numpy as np
import socketio
import splines
from znframe import Frame as ZnFrame

from .data import RoomSetData
from .utils import (
    estimate_max_batch_size_for_socket,
    split_list_into_chunks,
)

log = logging.getLogger(__name__)


class FrozenZnDraw:
    # TODO: take in _original instead and change the RoomSetData to accept token instead. Will remove a lot of the boilerplate for pushing data
    def __init__(self, token, url, cached_data: dict):
        self.socket = socketio.Client()
        self.socket.connect(url, wait_timeout=1)
        self.socket.emit("join", str(token))
        self.url = url
        self.token = token
        self._cached_data = self.cache_from_dict(cached_data)

    def __len__(self):
        return self._cached_data["length"]

    def __setitem__(self, index, value):
        assert isinstance(index, int), "Index must be an integer"
        if isinstance(value, ase.Atoms):
            value = ZnFrame.from_atoms(value)
        self.set_data(
            frames={index: value.to_dict(built_in_types=False)},
            step=index,
            update_database=True,
        )

    def __delitem__(self, index: int | slice | list[int]):
        if (
            isinstance(index, int)
            or isinstance(index, slice)
            or isinstance(index, list)
        ):
            length = len(self)
            index = self.wrap_and_check_index(index, length)
            self.set_data(frames={i: None for i in index}, update_database=True)
            self._cached_data["length"] -= len(index)
        else:
            raise TypeError("Index must be an integer, slice or list[int]")

    def append(self, value: t.Union[ase.Atoms, ZnFrame]) -> None:
        """Append atoms to the end of the list"""
        if isinstance(value, ase.Atoms):
            value = ZnFrame.from_atoms(value)
        self[len(self)] = value
        self._cached_data["length"] += 1

    def extend(
        self, values: t.Union[ase.Atoms, ZnFrame, list[ase.Atoms], list[ZnFrame]]
    ) -> None:
        """Extend the list by appending all the items in the given list"""
        size = len(self)
        if not isinstance(values, list):
            self.append(values)
        else:
            if isinstance(values[0], ase.Atoms):
                values = [ZnFrame.from_atoms(val) for val in values]
            batch_size = estimate_max_batch_size_for_socket(values)
            log.critical(f"Batch size: {batch_size}")
            indices = list(range(size, size + len(values)))
            all_data = [
                (i, val.to_dict(built_in_types=False))
                for i, val in zip(indices, values)
            ]
            chunks = list(split_list_into_chunks(all_data, batch_size))
            for i, chunk in enumerate(chunks):
                self._cached_data["length"] += len(chunk)
                batch = {tup[0]: tup[1] for tup in chunk}
                if i == len(chunks) - 1:
                    # Only send the step if it's the last batch, otherwise it will be set to the last index of the batch
                    self.set_data(frames=batch, update_database=True, step=len(self)-1)
                else:
                    self.set_data(frames=batch, update_database=True)
            

    def log(self, message: str) -> None:
        """Log a message to the console"""
        print(message)
        self.socket.emit(
            "message:log",
            {
                "message": message,
                "token": self.token,
            },
        )

    def set_data(self, **data: dict) -> None:
        data = RoomSetData(**data)
        log.critical("Trying to set data")
        self.socket.emit("room:set", data.to_dict())

    @property
    def atoms(self):
        return self._cached_data["atoms"]

    @property
    def points(self):
        return self._cached_data["points"]

    @points.setter
    def points(self, value):
        self.set_data(points=value, update_database=True)

    @property
    def segments(self):
        return self.calculate_segments(self._cached_data["points"])

    @property
    def step(self):
        return self._cached_data["step"]

    @step.setter
    def step(self, value):
        self.set_data(step=value, update_database=True)

    @property
    def selection(self):
        return self._cached_data["selection"]

    @selection.setter
    def selection(self, value):
        self.set_data(selection=value, update_database=True)

    @property
    def bookmarks(self):
        return self._cached_data["bookmarks"]

    @bookmarks.setter
    def bookmarks(self, value):
        self.set_data(bookmarks=value, update_database=True)

    @staticmethod
    def calculate_segments(points: np.ndarray) -> np.ndarray:
        if points.shape[0] <= 1:
            return points
        t = np.linspace(0, len(points) - 1, len(points) * 50)
        return splines.CatmullRom(points).evaluate(t)

    @staticmethod
    def cache_from_dict(data: dict):
        return dict(
            atoms=ZnFrame.from_dict(data["frames"][0]).to_atoms(),
            points=np.array(data["points"]),
            step=data["step"],
            selection=data["selection"],
            bookmarks=data["bookmarks"],
            length=data["length"],
        )

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
