from typing import List, Union

import ase
import numpy as np
from znframe.frame import Frame as ZnFrame

from .base import ZnDrawBase
from .db import Session
from .db.schema import Frame, Room
from .server.utils import add_frames_to_room, get_room_by_token


class ZnDrawWorker(ZnDrawBase):
    def __len__(self) -> int:
        with Session() as session:
            room = session.query(Room).get(self.token)
            return len(room.frames)

    def __setitem__(self, index: int, frame: ZnFrame):
        add_frames_to_room(self.token, (index, frame))

    def __getitem__(self, index: int) -> ase.Atoms:
        with Session() as session:
            room = session.query(Room).get(self.token)
            frame = session.query(Frame).filter_by(index=index, room=room).first()
            if frame is None:
                raise IndexError(f"Index {index} not found")
            return ZnFrame.from_dict(frame.data).to_atoms()

    def __delitem__(self, index: int):
        with Session() as session:
            room = session.query(Room).get(self.token)
            frame = room.frames.filter_by(index=index).first()
            if frame is None:
                raise IndexError(f"Index {index} not found")
            session.delete(frame)
            session.commit()

    def append(self, atoms: ase.Atoms):
        frame = ZnFrame.from_atoms(atoms)
        add_frames_to_room(self.token, (len(self), frame))

    def extend(self, atoms_list: list[ase.Atoms]):
        frames = [ZnFrame.from_atoms(atoms) for atoms in atoms_list]
        add_frames_to_room(
            self.token, [(len(self) + i, frame) for i, frame in enumerate(frames)]
        )

    @property
    def points(self) -> np.ndarray:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return np.array(room.points)

    @points.setter
    def points(self, value: Union[np.ndarray, list]):
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if len(value) > 0:
            try:
                assert len(value[0]) == 3
            except (TypeError, AssertionError):
                raise ValueError("Points must be a list of 3D coordinates")

        with Session() as session:
            room = get_room_by_token(session, self.token)
            room.points = value  # type: ignore
            session.commit()

        self.socket.emit("points:set", value)

    @property
    def segments(self) -> np.ndarray:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return np.array(room.segments)

    @property
    def step(self) -> int:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return room.currentStep

    @step.setter
    def step(self, idx: int):
        if idx < 0 or idx >= len(self):
            raise IndexError(f"Index {idx} out of range")
        with Session() as session:
            room = get_room_by_token(session, self.token)
            room.step = idx
            session.commit()

        self.socket.emit("room:set", idx)

    @property
    def selection(self) -> Union[List[int], List[None]]:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return room.selection

    @selection.setter
    def selection(self, value: Union[List[int], List[None]]):
        with Session() as session:
            room = get_room_by_token(session, self.token)
            room.selection = value
            session.commit()

        self.socket.emit("selection:set", value)

    @property
    def bookmarks(self) -> dict:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return room.bookmarks

    @bookmarks.setter
    def bookmarks(self, value: dict):
        with Session() as session:
            room = get_room_by_token(session, self.token)
            room.bookmarks = value
            session.commit()

        self.socket.emit("bookmarks:set", value)

    def insert(self, index: int, atoms: ase.Atoms):
        raise NotImplementedError

    def log(self, message: str):
        raise NotImplementedError
