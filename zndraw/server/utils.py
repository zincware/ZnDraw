
from sqlalchemy.orm import Session

from ..db import Session as ses
from ..db.schema import Queue, QueueItem, Room


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


def get_room_by_token(session: Session, token: str):
    return session.query(Room).filter_by(token=token).one()
