import unittest.mock as mock

import numpy.testing as npt
import pytest
import znframe
from ase.collections import s22
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from zndraw.db import schema
from zndraw.db.schema import Base
from zndraw.zndraw_worker import ZnDrawWorker

s22 = list(s22)


@pytest.fixture
def session() -> sessionmaker:
    """pytest fixture to setup the database"""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)


@pytest.fixture
def room_session(session) -> sessionmaker:
    with session() as s:
        room = schema.Room(
            token="test_token",
            currentStep=5,
            selection=[1, 2],
            points=[[0, 0, 0], [1, 1, 1]],
        )
        s.add(room)
        for idx, atoms in enumerate(s22):
            frame = schema.Frame(
                data=znframe.Frame.from_atoms(atoms).to_dict(built_in_types=False),
                index=idx,
                room=room,
            )
            s.add(frame)

        # add two bookmarks
        bookmark1 = schema.Bookmark(step=1, text="bm-1", room=room)
        s.add(bookmark1)
        bookmark2 = schema.Bookmark(step=2, text="bm-2", room=room)
        s.add(bookmark2)
        s.commit()

    return session


def test_room(room_session):
    with room_session() as s:
        room = s.query(schema.Room).filter_by(token="test_token").first()
        assert room.token == "test_token"
        assert room.currentStep == 5
        assert room.selection == [1, 2]


def test_frame(room_session):
    with room_session() as s:
        room = s.query(schema.Room).filter_by(token="test_token").first()
        assert len(room.frames) == 22
        for idx, (frame, atoms) in enumerate(zip(room.frames, s22)):
            assert frame.room == room
            assert frame.to_frame() == znframe.Frame.from_atoms(atoms)
            assert frame.index == idx


def test_zndraw_worker_get(room_session):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", socket=None)
        assert len(worker) == 22
        for idx, atoms in enumerate(s22):
            assert worker[idx] == atoms
        for idx, atoms in enumerate(worker):
            assert atoms == s22[idx]
        assert worker.step == 5
        assert worker.selection == [1, 2]
        npt.assert_array_equal(worker.points, [[0, 0, 0], [1, 1, 1]])
        assert worker.bookmarks == {1: "bm-1", 2: "bm-2"}

        with pytest.raises(IndexError):
            worker[22]

        # TODO: test slicing


def test_zndraw_worker_set_atoms(room_session):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", socket=None)
        del worker[:]
        assert len(worker) == 0
        worker[0] = s22[0]
        assert len(worker) == 1
        assert worker[0] == s22[0]
        worker.extend(s22[1:3])
        assert len(worker) == 3
        assert worker[1] == s22[1]
        assert worker[2] == s22[2]
        assert worker[:] == s22[0:3]
        worker.append(s22[3])
        assert len(worker) == 4
        assert worker[3] == s22[3]
        del worker[1]
        assert len(worker) == 3
        assert worker[:] == [s22[0], s22[2], s22[3]]
        worker[0] = s22[1]
        assert len(worker) == 3
        assert worker[:] == [s22[1], s22[2], s22[3]]

        worker.insert(0, s22[0])
        assert len(worker) == 4
        assert worker[:] == s22[0:4]

        with pytest.raises(IndexError):
            worker.insert(5, s22[0])


def test_zndraw_worker_set(room_session):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", socket=None)
        assert worker.step == 5
        worker.step = 6
        assert worker.step == 6
        assert worker.selection == [1, 2]
        worker.selection = [2, 3]
        assert worker.selection == [2, 3]
        npt.assert_array_equal(worker.points, [[0, 0, 0], [1, 1, 1]])
        worker.points = [[1, 1, 1], [2, 2, 2]]
        npt.assert_array_equal(worker.points, [[1, 1, 1], [2, 2, 2]])
        assert worker.bookmarks == {1: "bm-1", 2: "bm-2"}
        worker.bookmarks = {2: "bm-3", 3: "bm-4"}
        assert worker.bookmarks == {2: "bm-3", 3: "bm-4"}
