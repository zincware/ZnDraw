from sqlalchemy.orm import Session
from sqlalchemy import Engine
from ..db.schema import Queue

def get_queue(session: Session, queue_name: str) -> Queue:
    queue = session.query(Queue).filter_by(name=queue_name).first()
    if queue is None:
        queue = Queue(name=queue_name)
        session.add(queue)
        session.commit()
    return queue

def insert_into_queue(engine:Engine, queue_name:str, job_id:str) -> None:
    with Session(engine) as session:
        queue = get_queue(session, queue_name)
        queue.jobs.append(job_id)
        session.commit()
        
def remove_job_from_queue(engine:Engine, queue_name:str, job_id:str) -> None:
    with Session(engine) as session:
        queue = get_queue(session, queue_name)
        queue.jobs.remove(job_id)
        session.commit()