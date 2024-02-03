from sqlalchemy.orm import Session
from typing import Tuple

from ..db import Session as ses
from ..db.schema import Queue, QueueItem, Room, Frame
from znframe import Frame as ZnFrame

def get_queue(session: Session, queue_name: str) -> Queue:
    queue = session.query(Queue).filter_by(name=queue_name).first()
    if queue is None:
        queue = Queue(name=queue_name)
        session.add(queue)
        session.commit()
    return queue


def insert_into_queue(queue_name: str, job_id: str) -> None:
    with ses() as session:
        queue = get_queue(session, queue_name)
        job = QueueItem(job_id=job_id)
        queue.jobs.append(job)
        session.commit()


def remove_job_from_queue(queue_name: str, job_id: str) -> None:
    with ses() as session:
        queue = get_queue(session, queue_name)
        for job in queue.jobs:
            if job.job_id == job_id:
                queue.jobs.remove(job)

        session.commit()


def get_queue_position(queue_name: str, job_id: str) -> int:
    with ses() as session:
        queue = get_queue(session, queue_name)
        for i, job in enumerate(queue.jobs):
            if job.job_id == job_id:
                return i
        return -1

def get_room_by_token(session:Session, token:str):
    return session.query(Room).filter_by(token=token).one()

# write custom type that corresponds to a tuple of (int, ZnFrame)
frame_data = Tuple[int, ZnFrame]

def add_frames_to_room(room_token:str, data: frame_data | list[frame_data]):
    if isinstance(data, (int, ZnFrame)):
        list_data:list[frame_data] = [data] # type: ignore
    else:
        list_data = data # type: ignore
    with ses() as session:
        room = get_room_by_token(session, room_token)
        for idx, frame in list_data:
            old_frame = room.frames.filter_by(index=idx).first()
            if old_frame is not None:
                old_frame.data = frame.to_dict(built_in_types=False)
            else:
                new_frame = Frame(index=idx, data=frame.to_dict(built_in_types=False))
                room.frames.append(new_frame)