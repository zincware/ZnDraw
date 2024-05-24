import logging

import ase.io
import tqdm
import znframe
from celery import shared_task
from redis import Redis

from zndraw.base import FileIO

log = logging.getLogger(__name__)


@shared_task
def read_file(data: dict) -> None:
    file_io = FileIO(**data)
    r = Redis(host="localhost", port=6379, db=0, decode_responses=True)
    # r = znsocket.Client("http://127.0.0.1:5000")
    r.delete("room:default:frames")

    for i, atoms in tqdm.tqdm(enumerate(ase.io.iread(file_io.name))):
        if file_io.stop is not None and i >= file_io.stop:
            break
        frame = znframe.Frame.from_atoms(atoms)
        r.hset("room:default:frames", f"{i}", frame.to_json())
