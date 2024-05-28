import logging

import ase.io
import tqdm
import znframe
from celery import shared_task
from redis import Redis
from socketio import SimpleClient

from zndraw.base import FileIO, RedisList

log = logging.getLogger(__name__)


@shared_task
def read_file(fileio: dict, io_port: int, storage: str) -> None:
    file_io = FileIO(**fileio)
    # r = Redis(host="localhost", port=6379, db=0, decode_responses=True)

    if storage.startswith("redis"):
        r = Redis.from_url(storage, decode_responses=True)
    elif storage.startswith("znsocket"):
        raise NotImplementedError

    io = SimpleClient()

    # r = znsocket.Client("http://127.0.0.1:5000")

    # TODO: make everyone join room main
    # send update here to everyone in room, because this is only called once in the beginnig
    # chain this with compute_bonds. So this will load much faster
    r.delete("room:default:frames")

    lst = RedisList(r, "room:default:frames")

    for i, atoms in tqdm.tqdm(enumerate(ase.io.iread(file_io.name))):
        if file_io.stop is not None and i >= file_io.stop:
            break
        frame = znframe.Frame.from_atoms(atoms)
        # r.hset("room:default:frames", f"{i}", frame.to_json())
        # r.rpush("room:default:frames", frame.to_json())
        lst.append(frame.to_json())
        if i == 0:
            io.connect(f"http://127.0.0.1:{io_port}")
            io.emit("room:all:frames:refresh", [0])


@shared_task
def run_modifier(url, room, data: dict) -> None:
    from zndraw import ZnDraw
    from zndraw.modify import Modifier

    # cls = get_cls_from_json_schema(modifier["schema"], modifier["name"])
    vis = ZnDraw(url=url, token=room)
    vis.socket.emit("modifier:run:running")
    try:
        modifier = Modifier(**data)
        modifier.run(vis)
    finally:
        # vis[list(range(10))] = [ase.build.molecule("H2O") for _ in range(10)]
        # TODO: why is everything after 10 configuration removed?
        vis.socket.emit("modifier:run:finished")

    # 1. run modifier and update redis
    # 2. use to update frames in real time "room:frames:refresh"


@shared_task
def run_selection(url, room, data: dict) -> None:
    from zndraw import ZnDraw
    from zndraw.select import Selection

    vis = ZnDraw(url=url, token=room)
    vis.socket.emit("selection:run:running")
    try:
        selection = Selection(**data)
        selection.run(vis)
    finally:
        vis.socket.emit("selection:run:finished")


@shared_task
def run_analysis(url, room, data: dict) -> None:
    from zndraw import ZnDraw
    from zndraw.analyse import Analysis

    vis = ZnDraw(url=url, token=room)
    vis.socket.emit("analysis:run:running")
    try:
        analysis = Analysis(**data)
        analysis.run(vis)
    finally:
        vis.socket.emit("analysis:run:finished")
