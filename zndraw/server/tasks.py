from dataclasses import asdict

import ase.io
import znframe
from celery import shared_task
from flask import current_app
from socketio import Client
import pathlib
import tqdm
import znh5md

from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig

from ..app import cache
from .data import CeleryTaskData, FrameData


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
    con = get_client(url)
    fileio = cache.get("FILEIO")

    msg = CeleryTaskData(
        target=target,
        event="atoms:upload",
        data=FrameData(
            index=0,
            data=znframe.Frame.from_atoms(ase.Atoms()).to_dict(built_in_types=False),
            update=True,
        ),
    )
    con.emit("celery:task:results", asdict(msg))
    if fileio.name is None:
        return
    
    if fileio.remote is not None:
        node_name, attribute = fileio.name.split(".", 1)
        try:
            import zntrack

            node = zntrack.from_rev(
                node_name, remote=fileio.remote, rev=fileio.rev
            )
            generator = getattr(node, attribute)
        except ImportError as err:
            raise ImportError(
                "You need to install ZnTrack to use the remote feature"
            ) from err
    elif pathlib.Path(fileio.name).suffix == ".h5":
        reader = znh5md.ASEH5MD(fileio.name)
        generator = reader.get_atoms_list()
    else:
        generator = ase.io.iread(fileio.name)

    frame = 0

    for idx, atoms in tqdm.tqdm(enumerate(generator), ncols=100):
        if fileio.start and idx < fileio.start:
            continue
        if fileio.stop and idx >= fileio.stop:
            break
        if fileio.step and idx % fileio.step != 0:
            continue
        
        msg = CeleryTaskData(
            target=target,
            event="atoms:upload",
            data=FrameData(
                index=idx,
                data=znframe.Frame.from_atoms(atoms).to_dict(built_in_types=False),
                update=True,
            ),
        )
        con.emit("celery:task:results", asdict(msg))

        frame += 1
