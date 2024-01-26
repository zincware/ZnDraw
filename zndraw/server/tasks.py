from dataclasses import asdict

import ase.io
import znframe
from celery import shared_task
from socketio import Client
from flask_caching import Cache
from zndraw.app import create_app, cache
import tqdm

from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig

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
    FILENAME = "/Users/fzills/tools/ZnDraw/tmp/BMIM_BF4_303_15K.extxyz"
    con = get_client(url)
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
