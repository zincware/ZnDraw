import dataclasses
import logging
from typing import List, Union

import ase
import numpy as np
from sqlalchemy import func as sql_func
from znframe.frame import Frame as ZnFrame

from zndraw.data import RoomSetData

log = logging.getLogger(__name__)

from .base import ZnDrawBase
from .db import Session
from .db.schema import Bookmark, Frame, Room
from .server.utils import get_room_by_token
from .utils import (
    check_selection,
    estimate_max_batch_size_for_socket,
    split_list_into_chunks,
    wrap_and_check_index,
)


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
            return self._get_len(session, self.token)

    @staticmethod
    def _get_len(session, token) -> int:
        max_idx = (
            session.query(sql_func.max(Frame.index))
            .filter(Frame.room_token == token)
            .scalar()
        )
        return max_idx + 1 if max_idx is not None else 0

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
        # TODO: we need the batch size here!
        if self.emit:
            frames = {idx: frame for idx, frame in zip(index, value)}
            chunk_size = estimate_max_batch_size_for_socket(list(frames.values()))
            for frame_ids in split_list_into_chunks(list(frames), chunk_size):
                msg = RoomSetData(
                    frames={
                        idx: frames[idx].to_dict(built_in_types=False)
                        for idx in frame_ids
                    },
                )
                self.socket.emit(
                    "room:set",
                    msg.to_dict(),
                )
            self.socket.emit(
                "room:set",
                RoomSetData(
                    step=index[-1],
                ).to_dict(),
            )
        with Session() as session:
            room = session.query(Room).get(self.token)
            room.currentStep = index[-1]
            self.write_frames_to_session(session, index, value, room)
            session.commit()

    @staticmethod
    def write_frames_to_session(session, index, values, room):
        for _index, _value in zip(index, values):
            frame = session.query(Frame).filter_by(index=_index, room=room).first()
            if frame is None:
                # add new frame
                frame = Frame(index=_index, room=room)
                session.add(frame)
            frame.data = _value.to_dict(built_in_types=False)

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
            room = get_room_by_token(session, self.token)
            self.delete_indices_from_session(session, index, room)
            session.commit()
        if self.emit:
            self.socket.emit(
                "room:set", RoomSetData(frames={idx: None for idx in index}).to_dict()
            )

    @staticmethod
    def delete_indices_from_session(session, indices, room):
        session.query(Frame).filter(
            Frame.index.in_(indices), Frame.room == room
        ).delete(synchronize_session=False)
        # ensure indices are contiguous
        frames = session.query(Frame).filter_by(room=room).all()
        for idx, frame in enumerate(frames):
            frame.index = idx

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
        return self.calculate_segments(self.points)

    @property
    def step(self) -> int:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return room.currentStep

    @step.setter
    def step(self, idx: int):
        length = len(self)
        if idx < 0 or idx >= length:
            raise IndexError(f"Index {idx} out of range for {length} frames")
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
        check_selection(value, len(self.atoms))
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

    @property
    def camera(self) -> dict:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            return room.camera

    @camera.setter
    def camera(self, payload: dict) -> None:
        expected_keys = {"position", "target"}
        if set(payload) != expected_keys:
            raise ValueError("camera must have keys 'position' and 'target'")
        for k in expected_keys:
            v = payload[k]
            if isinstance(v, np.ndarray):
                payload[k] = v.tolist()

        self.socket.emit("camera:update", payload)

        with Session() as session:
            room = get_room_by_token(session, self.token)
            room.camera = payload
            session.commit()

    @bookmarks.setter
    def bookmarks(self, value: dict) -> None:
        with Session() as session:
            room = get_room_by_token(session, self.token)
            self.write_bookmark_dictionary_to_db(session, room, value)
            session.commit()
        if self.emit:
            self.socket.emit("room:set", RoomSetData(bookmarks=value).to_dict())

    @staticmethod
    def write_bookmark_dictionary_to_db(session, room, value):
        # delete all bookmarks
        session.query(Bookmark).filter_by(room=room).delete()
        for step, text in value.items():
            bookmark = Bookmark(step=step, text=text, room=room)
            session.add(bookmark)

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
                        indices = wrap_and_check_index(
                            indices, self._get_len(session, self.token)
                        )
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
                        answer["length"] = self._get_len(session, self.token)
                    elif key == "points":
                        answer["points"] = room.points
                    elif key == "step":
                        answer["step"] = room.currentStep
                    elif key == "selection":
                        answer["selection"] = room.selection
                    elif key == "camera":
                        answer["camera"] = room.camera
                    elif key == "bookmarks":
                        answer["bookmarks"] = {
                            bm.step: bm.text for bm in room.bookmarks
                        }
        return answer

    def set_properties(self, **kwargs):
        with Session() as session:
            room = get_room_by_token(session, self.token)
            for key, payload in kwargs.items():
                if payload is not None:
                    if key == "frames":
                        # values: dict{int:Frame}
                        if len({type(frame) for frame in payload.values()}) != 1:
                            raise ValueError("All frames must be None or not None")
                        is_removing = all(frame is None for frame in payload.values())
                        indices = list(payload.keys())
                        if is_removing:
                            self.delete_indices_from_session(session, indices, room)
                        else:
                            frames = [
                                ZnFrame.from_dict(frame) for frame in payload.values()
                            ]
                            self.write_frames_to_session(session, indices, frames, room)
                    elif key == "points":
                        room.points = payload
                    elif key == "step":
                        room.currentStep = payload
                    elif key == "selection":
                        room.selection = payload
                    elif key == "camera":
                        room.camera = payload
                    elif key == "bookmarks":
                        self.write_bookmark_dictionary_to_db(session, room, payload)
            session.commit()

    def log(self, message: str):
        self.socket.emit("message:log", {"message": message, "token": self.token})

    def upload(self, target: str):
        """Emit all frames to the target (webclient)."""
        if not self.emit:
            raise ValueError("Emit is disabled")

        self.extend(list(self))
