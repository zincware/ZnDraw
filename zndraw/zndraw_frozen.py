import logging
import typing as t

import ase
import numpy as np
import socketio
from znframe import Frame as ZnFrame

from .base import ZnDrawBase
from .data import RoomSetData
from .utils import (
    check_selection,
    estimate_max_batch_size_for_socket,
    split_list_into_chunks,
    wrap_and_check_index,
)

log = logging.getLogger(__name__)


class ZnDrawFrozen(ZnDrawBase):
    # TODO: take in _original instead and change the RoomSetData to accept token instead. Will remove a lot of the boilerplate for pushing data
    def __init__(self, token, url, cached_data: dict):
        self.socket = socketio.Client()
        self.socket.connect(url, wait_timeout=5)
        self.socket.emit("join", {"token": str(token), "auth_token": None})
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

    def __getitem__(self, *args, **kwargs):
        raise NotImplementedError(
            "This method is not implemented on the frozen object. Please use the zndraw.ZnDraw object instead."
        )

    def insert(self, *args, **kwargs):
        raise NotImplementedError("Don't use inserts. Use append or extend instead.")

    def __delitem__(self, index: int | slice | list[int]):
        if (
            isinstance(index, int)
            or isinstance(index, slice)
            or isinstance(index, list)
        ):
            length = len(self)
            index = wrap_and_check_index(index, length)
            self.set_data(frames={i: None for i in index}, update_database=True)
            self._cached_data["length"] -= len(index)
            for idx in index:
                self.bookmarks.pop(idx, None)
                self.bookmarks.pop(str(idx), None)
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
            log.debug(f"Batch size: {batch_size}")
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
                    self.set_data(
                        frames=batch, update_database=True, step=len(self) - 1
                    )
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
        self.socket.emit("room:set", data.to_dict())

    @property
    def atoms(self):
        return self._cached_data["atoms"]

    @property
    def points(self):
        return self._cached_data["points"]

    @points.setter
    def points(self, value):
        if isinstance(value, np.ndarray):
            value = value.tolist()
        # TODO: add extra checks for lengths of sublists
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
        check_selection(value, len(self.atoms))
        self.set_data(selection=value, update_database=True)

    @property
    def bookmarks(self) -> dict[int | str, str]:
        return self._cached_data["bookmarks"]

    @bookmarks.setter
    def bookmarks(self, value):
        self.set_data(bookmarks=value, update_database=True)

    @staticmethod
    def cache_from_dict(data: dict):
        return dict(
            atoms=ZnFrame.from_dict(data["frames"][0]).to_atoms(),
            points=np.array(data["points"]),
            step=data["step"],
            selection=data["selection"],
            bookmarks=data["bookmarks"],
            length=data["length"],
            camera=data["camera"],
        )

    @property
    def camera(self):
        return self._cached_data["camera"]

    @camera.setter
    def camera(self, value):
        if set(value) != {"position", "target"}:
            raise ValueError("camera must have keys 'position' and 'target'")
        self.set_data(camera=value, update_database=True)
