from dataclasses import asdict

import ase.io
import znframe
from celery import shared_task
from flask import current_app
from socketio import Client

from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig

from ..app import cache
from .data import CeleryTaskData


def get_client(url) -> Client:
    client = Client()
    client.connect(url)
    return client


@shared_task
def get_selection_schema(url: str, target: str):
    print(f"emitting selection_schema to {target}")

    config = GlobalConfig.load()
    cls = get_selection_class(config.get_selection_methods())
    schema = cls.model_json_schema()
    msg = CeleryTaskData(
        target=target,
        event="selection:schema",
        data=schema,
    )
    con = get_client(url)
    con.emit("celery:task:results", asdict(msg))


@shared_task
def read_file(url: str, target: str):
    fileio = cache.get("FILEIO")
    print(f"emitting {fileio} to {target}")
    FILENAME = "/home/rokas/Programming/data/zinc_fragments_difflinker.xyz"
    con = get_client(url)
    cache2 = next(iter(current_app.extensions["cache"]))
    print(cache2)
    print(f"file io from cache2: {cache2.get('FILEIO')}")
    for idx, atoms in enumerate(ase.io.iread(FILENAME)):
        frame = znframe.Frame.from_atoms(atoms)

        data = {
            idx: frame.to_dict(built_in_types=False),
            "display_new": True,
        }

        msg = CeleryTaskData(
            target=target,
            event="atoms:upload",
            data=data,
        )

        con.emit("celery:task:results", asdict(msg))
        if idx > 10:
            break
    print("done")
