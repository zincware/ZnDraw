import uuid
from datetime import datetime

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


def get_room_by_token(session: Session, token: str) -> Room:
    return session.query(Room).filter_by(token=token).one()


def insert_into_queue(queue_name: str, job_name: str, room_token: str, parameters: dict) -> str:
    job_id = uuid.uuid4().hex
    with ses() as session:
        queue = get_queue(session, queue_name)
        room = get_room_by_token(session, room_token)
        job = QueueItem(job_name=job_name, job_id=job_id, datetime=datetime.utcnow(), parameters=parameters, queue=queue, room=room)
        session.add(job)
        session.commit()
    return job_id


def update_job_status(job_id: str | int, status: str) -> None:
    with ses() as session:
        session.query(QueueItem).filter_by(job_id=job_id).update(dict(status=status))
        session.commit()


def get_queue_position(queue_name: str) -> list[tuple[int, str]]:
    results = []
    with ses() as session:
        queueing_jobs = (
            session.query(QueueItem)
            .filter_by(status="queued", queue_name=queue_name)
            .all()
        )
        rooms = [job.room_token for job in queueing_jobs]
        for room in rooms:
            largest_idx = max(
                [i + 1 for i, job in enumerate(queueing_jobs) if job.room_token == room]
            )
            results.append((largest_idx, room))
    return results
