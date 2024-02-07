import dataclasses
import logging
from typing import List, Union

import ase
import numpy as np
from znframe.frame import Frame as ZnFrame

from zndraw.data import CeleryTaskData, RoomSetData
from zndraw.settings import GlobalConfig

log = logging.getLogger(__name__)

from .base import ZnDrawBase
from .db import Session
from .db.schema import Bookmark, Frame, Room
from .server.utils import get_room_by_token
from .utils import check_selection, wrap_and_check_index


def _any_to_list(
    value: ZnFrame | ase.Atoms | list[ase.Atoms] | list[ZnFrame],
) -> list[ZnFrame]:
    if isinstance(value, ase.Atoms):
        return [ZnFrame.from_atoms(value)]
    elif isinstance(value, ZnFrame):
        return [value]
    elif isinstance(value, list):
        data = []
        for v in value:
            if isinstance(v, ase.Atoms):
                data.append(ZnFrame.from_atoms(v))
            elif isinstance(v, ZnFrame):
                data.append(v)
            else:
                raise ValueError(f"Invalid type for value {v}")
        return data
    raise ValueError("Invalid type for value")


@dataclasses.dataclass
class ZnDrawWorker(ZnDrawBase):
    emit: bool = True

    def __len__(self) -> int:
        with Session() as session:
            room = session.query(Room).get(self.token)
            return len(room.frames)

    def __setitem__(
        self,
        index: int | list[int] | slice,
        value: ZnFrame | ase.Atoms | list[ase.Atoms] | list[ZnFrame],
    ):
        value = _any_to_list(value)
        if isinstance(index, str):
            index = int(index)
        if isinstance(index, int):
            index = [index]
        if isinstance(index, slice):
            index = list(range(len(self)))[index]
        if len(index) != len(value):
            # TODO: support all the ways python lists can be slice set
            raise ValueError(
                f"Length of index ({len(index)}) and value ({len(value)}) must match"
            )

        # we emit first, because sending the data takes longer, but emit is faster
        if self.emit:
            self.socket.emit(
                "room:set",
                RoomSetData(
                    frames={
                        idx: frame.to_dict(built_in_types=False)
                        for idx, frame in zip(index, value)
                    },
                    step=len(self) - 1 + len(index),
                ).to_dict(),
            )
        with Session() as session:
            room = session.query(Room).get(self.token)
            for _index, _value in zip(index, value):
                frame = session.query(Frame).filter_by(index=_index, room=room).first()
                if frame is None:
                    # add new frame
                    frame = Frame(index=_index, room=room)
                    session.add(frame)
                frame.data = _value.to_dict(built_in_types=False)
            session.commit()

    def __getitem__(self, index: int | list[int] | slice) -> ase.Atoms:
        single_index = False
        if isinstance(index, int):
            index = [index]
            single_index = True
        if isinstance(index, slice):
            index = list(range(len(self)))[index]

        with Session() as session:
            room = session.query(Room).get(self.token)
            frames = (
                session.query(Frame)
                .filter_by(room=room)
                .filter(Frame.index.in_(index))
                .all()
            )
            if frames is None:
                raise IndexError(f"Index {index} not found")
            frames = sorted(frames, key=lambda f: f.index)
            data = [ZnFrame.from_dict(frame.data).to_atoms() for frame in frames]
            if single_index:
                return data[0]
            return data

    def __delitem__(self, index: int | list[int] | slice):
        if isinstance(index, int):
            index = [index]
        if isinstance(index, slice):
            index = list(range(len(self)))[index]

        with Session() as session:
            room = session.query(Room).get(self.token)
            for _index in index:
                session.query(Frame).filter_by(index=_index, room=room).delete()
            # ensure indices are contiguous
            frames = session.query(Frame).filter_by(room=room).all()
            for idx, frame in enumerate(frames):
                frame.index = idx
            session.commit()
        if self.emit:
            self.socket.emit(
                "room:set", RoomSetData(frames={idx: None for idx in index}).to_dict()
            )

    def append(self, data: ase.Atoms | ZnFrame):
        if isinstance(data, ase.Atoms):
            data = ZnFrame.from_atoms(data)
        self[len(self)] = data

    def extend(self, atoms_list: list[ase.Atoms] | list[ZnFrame]):
        indices = list(range(len(self), len(self) + len(atoms_list)))
        frames = _any_to_list(atoms_list)
        self[indices] = frames

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
        if self.emit:
            self.socket.emit("room:set", RoomSetData(points=value).to_dict())

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
            room.currentStep = idx
            session.commit()
        if self.emit:
            self.socket.emit("room:set", RoomSetData(step=idx).to_dict())

    @property
    def selection(self) -> Union[List[int], List[None]]:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return room.selection

    @selection.setter
    def selection(self, value: list[int]):
        check_selection(value)
        with Session() as session:
            room = get_room_by_token(session, self.token)
            room.selection = value
            session.commit()
        if self.emit:
            self.socket.emit("room:set", RoomSetData(selection=value).to_dict())

    @property
    def bookmarks(self) -> dict:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return {bm.step: bm.text for bm in room.bookmarks}

    @bookmarks.setter
    def bookmarks(self, value: dict) -> None:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            # delete all bookmarks
            session.query(Bookmark).filter_by(room=room).delete()
            for step, text in value.items():
                bookmark = Bookmark(step=step, text=text, room=room)
                session.add(bookmark)
            session.commit()
        if self.emit:
            self.socket.emit("room:set", RoomSetData(bookmarks=value).to_dict())

    def insert(self, index: int, atoms: ase.Atoms | ZnFrame):
        if index < 0 or index > len(self):
            raise IndexError(f"Index {index} out of range")
        data_after = self[index:]
        self[index] = atoms
        del self[index + 1 :]
        self.extend(data_after)

    def get_properties(self, **kwargs):
        with Session() as session:
            room = get_room_by_token(session, self.token)

            answer = {}
            for key, collect in kwargs.items():
                if collect:
                    if key == "frames":
                        indices = [
                            x if not x == "current" else room.currentStep
                            for x in kwargs["frames"]
                        ]
                        indices = wrap_and_check_index(indices, len(room.frames))
                        log.critical(f"Indices: {indices}")
                        collected_frames = (
                            session.query(Frame)
                            .filter(Frame.index.in_(indices), Frame.room == room)
                            .all()
                        )
                        answer["frames"] = [
                            ZnFrame.from_dict(frame.data).to_dict(built_in_types=False)
                            for frame in collected_frames
                        ]
                    if key == "length":
                        answer["length"] = len(room.frames)
                    elif key == "points":
                        answer["points"] = room.points
                    elif key == "step":
                        answer["step"] = room.currentStep
                    elif key == "selection":
                        answer["selection"] = room.selection
                    elif key == "bookmarks":
                        answer["bookmarks"] = {
                            bm.step: bm.text for bm in room.bookmarks
                        }
        return answer

    def log(self, message: str):
        self.socket.emit("message:log", {"message": message, "token": self.token})

    def upload(self, target: str):
        """Emit all frames to the target (webclient)."""
        if not self.emit:
            raise ValueError("Emit is disabled")
        config = GlobalConfig.load()
        frame_list = []
        for idx in range(len(self)):
            frame_list.append(self[idx])
            if len(frame_list) == config.read_batch_size:
                msg = CeleryTaskData(
                    target=target,
                    event="room:set",
                    data=RoomSetData(
                        frames={
                            idx
                            - jdx: ZnFrame.from_atoms(atoms).to_dict(
                                built_in_types=False
                            )
                            for jdx, atoms in enumerate(reversed(frame_list))
                        },
                        step=idx if idx < self.step else self.step,
                    ).to_dict(),
                )

                self.socket.emit("celery:task:emit", msg.to_dict())
                frame_list = []

        msg = CeleryTaskData(
            target=target,
            event="room:set",
            data=RoomSetData(
                frames={
                    idx - jdx: ZnFrame.from_atoms(atoms).to_dict(built_in_types=False)
                    for jdx, atoms in enumerate(reversed(frame_list))
                },
                selection=self.selection,
                points=self.points.tolist(),
                bookmarks=self.bookmarks,
                step=self.step,
            ).to_dict(),
        )

        self.socket.emit("celery:task:emit", msg.to_dict())
