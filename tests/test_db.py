import unittest.mock as mock

import numpy.testing as npt
import pytest
import znframe
from ase.collections import s22
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from zndraw.data import RoomGetData, RoomSetData
from zndraw.db import schema
from zndraw.db.schema import Base
from zndraw.utils import typecast
from zndraw.zndraw_worker import ZnDrawWorker
from zndraw.app import ZnDrawServer

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

        # add two modifier clients
        client1 = schema.ModifierClient(sid="sid1", timeout=55, available=True)
        client2 = schema.ModifierClient(sid="sid2", timeout=66, available=False)
        s.add(client1)
        s.add(client2)

        # add 3 queue items one running, one failed, one queued
        queue = schema.Queue(name="test_queue")
        s.add(queue)
        item1 = schema.QueueItem(status="running", queue=queue)
        item2 = schema.QueueItem(status="failed", queue=queue)
        item3 = schema.QueueItem(status="queued", queue=queue)
        s.add(item1)
        s.add(item2)
        s.add(item3)

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


def test_zndraw_worker_get(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
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


def test_zndraw_worker_get_multiple(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
        request = RoomGetData(points=True, bookmarks=True, selection=True, step=True)
        answer = worker.get_properties(**request.to_dict())
        npt.assert_array_equal(answer["points"], [[0, 0, 0], [1, 1, 1]])
        assert answer["bookmarks"] == {1: "bm-1", 2: "bm-2"}
        assert answer["selection"] == [1, 2]
        assert answer["step"] == 5

        request_idx = [0, 1, 5]
        request = RoomGetData(
            points=True, bookmarks=True, selection=True, step=True, frames=request_idx
        )
        answer = worker.get_properties(**request.to_dict())
        assert len(answer["frames"]) == len(request_idx)
        for i, frame in enumerate(answer["frames"]):
            atoms = znframe.Frame.from_dict(frame).to_atoms()
            assert atoms == s22[request_idx[i]]


def test_zndraw_worker_set_atoms(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
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


def test_zndraw_worker_set(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
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


def test_set_bookmarks(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
        global answer
        answer = None

        @typecast
        def on_answer(data: RoomSetData):
            global answer
            answer = data

        worker.socket.on("room:set", on_answer)

        worker.bookmarks = {2: "bm-3", 3: "bm-4"}
        worker.socket.sleep(0.1)
        # TODO: keys being converted to strings somewhere
        assert answer.bookmarks == {"2": "bm-3", "3": "bm-4"}
        assert worker.bookmarks == {2: "bm-3", 3: "bm-4"}


def test_set_step(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
        global answer
        answer = None

        @typecast
        def on_answer(data: RoomSetData):
            global answer
            answer = data

        worker.socket.on("room:set", on_answer)

        worker.step = 6
        worker.socket.sleep(0.1)
        assert answer.step == 6
        assert worker.step == 6


def test_set_selection(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
        global answer
        answer = None

        @typecast
        def on_answer(data: RoomSetData):
            global answer
            answer = data

        worker.socket.on("room:set", on_answer)

        worker.selection = [2, 3]
        worker.socket.sleep(0.1)
        assert answer.selection == [2, 3]
        assert worker.selection == [2, 3]


def test_set_points(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
        global answer
        answer = None

        @typecast
        def on_answer(data: RoomSetData):
            global answer
            answer = data

        worker.socket.on("room:set", on_answer)

        worker.points = [[1, 1, 1], [2, 2, 2]]
        worker.socket.sleep(0.1)
        npt.assert_array_equal(answer.points, [[1, 1, 1], [2, 2, 2]])
        npt.assert_array_equal(worker.points, [[1, 1, 1], [2, 2, 2]])


def test_set_atoms(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
        global answer
        answer = []

        @typecast
        def on_answer(data: RoomSetData):
            global answer
            answer.append(data)

        worker.socket.on("room:set", on_answer)

        worker[0] = s22[10]
        worker.socket.sleep(0.1)
        data = answer[0]
        assert answer[1].step == 0  # The step is sent after the data
        answer = []
        # TODO: keys being converted to strings somewhere
        assert data.frames["0"] == znframe.Frame.from_atoms(s22[10]).to_dict(
            built_in_types=False
        )
        assert len(data.frames) == 1
        assert worker[0] == s22[10]
        assert len(worker) == 22

        worker.extend(s22[11:13])
        worker.socket.sleep(0.1)
        data = answer[0]
        assert answer[1].step == 23
        answer = []
        assert len(data.frames) == 2
        assert data.frames.keys() == {str(len(worker) - 1), str(len(worker) - 2)}
        assert len(worker) == 24
        assert data.frames["22"] == znframe.Frame.from_atoms(s22[11]).to_dict(
            built_in_types=False
        )
        assert data.frames["23"] == znframe.Frame.from_atoms(s22[12]).to_dict(
            built_in_types=False
        )
        assert worker[22] == s22[11]
        assert worker[23] == s22[12]


def test_del_atoms(room_session, sio_server):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        worker = ZnDrawWorker(token="test_token", url=sio_server)
        global answer
        answer = None

        @typecast
        def on_answer(data: RoomSetData):
            global answer
            answer = data

        worker.socket.on("room:set", on_answer)

        del worker[0]
        worker.socket.sleep(0.1)
        assert len(answer.frames) == 1
        assert answer.frames.keys() == {"0"}
        assert answer.frames["0"] == None
        assert len(worker) == 21

        del worker[4:7]
        worker.socket.sleep(0.1)
        assert len(answer.frames) == 3
        assert answer.frames.keys() == {"4", "5", "6"}
        assert answer.frames["4"] == None
        assert answer.frames["5"] == None
        assert answer.frames["6"] == None
        assert len(worker) == 18


def test_app_startup_removes_modifier_clients(room_session):
    with room_session() as s:
        assert s.query(schema.ModifierClient).count() == 2
    ZnDrawServer._purge_old_modifier_clients(room_session)
    with room_session() as s:
        assert s.query(schema.ModifierClient).count() == 0


def test_app_startup_removes_queue_items(room_session):
    with room_session() as s:
        assert s.query(schema.QueueItem).count() == 3
        assert s.query(schema.QueueItem).filter_by(status="running").count() == 1
        assert s.query(schema.QueueItem).filter_by(status="failed").count() == 1
        assert s.query(schema.QueueItem).filter_by(status="queued").count() == 1
    ZnDrawServer._mark_old_queue_items_as_failed(room_session)
    with room_session() as s:
        assert s.query(schema.QueueItem).count() == 3
        assert s.query(schema.QueueItem).filter_by(status="running").count() == 1
        assert s.query(schema.QueueItem).filter_by(status="queued").count() == 0
        assert (
            s.query(schema.QueueItem)
            .filter(schema.QueueItem.status.like("%failed%"))
            .count()
            == 2
        )
