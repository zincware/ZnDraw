from celery import shared_task
import logging
import typing as t

from zndraw.base import FileIO
from redis import Redis
import ase.io
import znframe
import tqdm


log = logging.getLogger(__name__)


@shared_task
def read_file(data: dict) -> None:
    file_io = FileIO(**data)
    r = Redis(host='localhost', port=6379, db=0, decode_responses=True)

    r.delete("room:default:frames")

    for i, atoms in tqdm.tqdm(enumerate(ase.io.iread(file_io.name))):
        frame = znframe.Frame.from_atoms(atoms)
        r.hset("room:default:frames", f"{i}", frame.to_json())
        if i > file_io.stop:
            break
