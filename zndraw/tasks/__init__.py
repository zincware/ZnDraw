import logging
import typing as t

import ase.io
import tqdm
import znframe
import znsocket
from celery import shared_task
from redis import Redis
from socketio import SimpleClient

from zndraw.base import FileIO

log = logging.getLogger(__name__)


@shared_task
def run_znsocket_server(port: int) -> None:
    znsocket.Server(port=port).run()
    log.critical("ZnSocket server closed.")


@shared_task
def read_file(fileio: dict, io_port: int, storage: str) -> None:
    file_io = FileIO(**fileio)
    # r = Redis(host="localhost", port=6379, db=0, decode_responses=True)

    if storage.startswith("redis"):
        r = Redis.from_url(storage, decode_responses=True)
    elif storage.startswith("znsocket"):
        r = znsocket.Client.from_url(storage)

    io = SimpleClient()

    # r = znsocket.Client("http://127.0.0.1:5000")

    # TODO: make everyone join room main
    # send update here to everyone in room, because this is only called once in the beginning
    # chain this with compute_bonds. So this will load much faster
    r.delete("room:default:frames")

    lst = znsocket.List(r, "room:default:frames")

    if file_io.name is None:

        def _generator():
            yield ase.Atoms()

        generator = _generator()
    elif file_io.remote is not None:
        node_name, attribute = file_io.name.split(".", 1)
        try:
            import zntrack

            node = zntrack.from_rev(node_name, remote=file_io.remote, rev=file_io.rev)
            generator = getattr(node, attribute)
        except ImportError as err:
            raise ImportError(
                "You need to install ZnTrack to use the remote feature (or `pip install zndraw[all]`)."
            ) from err
    elif file_io.name.endswith((".h5", ".hdf5", ".h5md")):
        try:
            import znh5md

            reader = znh5md.ASEH5MD(file_io.name)
            generator = reader.get_atoms_list()
        except ImportError as err:
            raise ImportError(
                "You need to install ZnH5MD to use the remote feature (or `pip install zndraw[all]`)."
            ) from err
    else:
        generator = ase.io.iread(file_io.name)

    generator: t.Iterable[ase.Atoms]

    for idx, atoms in tqdm.tqdm(enumerate(generator)):
        if file_io.start and idx < file_io.start:
            continue
        if file_io.stop and idx >= file_io.stop:
            break
        if file_io.step and idx % file_io.step != 0:
            continue
        frame = znframe.Frame.from_atoms(atoms)
        lst.append(frame.to_json())
        if idx == 0:
            io.connect(f"http://127.0.0.1:{io_port}")
            io.emit("room:all:frames:refresh", [0])
        # TODO: trigger length refresh


@shared_task
def run_modifier(url, room, data: dict) -> None:
    from zndraw import ZnDraw
    from zndraw.modify import Modifier

    # cls = get_cls_from_json_schema(modifier["schema"], modifier["name"])
    vis = ZnDraw(url=url, token=room)
    vis.socket.emit("room:modifier:queue", 0)
    try:
        modifier = Modifier(**data)
        modifier.run(vis)
    except Exception as e:
        vis.log(str(e))
    finally:
        vis.socket.emit("room:modifier:queue", -1)

    # wait and then disconnect
    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_selection(url, room, data: dict) -> None:
    from zndraw import ZnDraw
    from zndraw.select import Selection

    vis = ZnDraw(url=url, token=room)
    vis.socket.emit("room:selection:queue", 0)
    try:
        selection = Selection(**data)
        selection.run(vis)
    finally:
        vis.socket.emit("room:selection:queue", -1)

    # wait and then disconnect
    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_analysis(url, room, data: dict) -> None:
    from zndraw import ZnDraw
    from zndraw.analyse import Analysis

    vis = ZnDraw(url=url, token=room)
    vis.socket.emit("room:analysis:queue", 0)
    try:
        analysis = Analysis(**data)
        analysis.run(vis)
    except Exception as e:
        vis.log(str(e))
    finally:
        vis.socket.emit("room:analysis:queue", -1)

    # wait and then disconnect
    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_geometry(url, room, data: dict) -> None:
    from zndraw import ZnDraw
    from zndraw.draw import Geometry

    vis = ZnDraw(url=url, token=room)
    vis.socket.emit("room:geometry:queue", 0)
    try:
        geom = Geometry(**data)
        # TODO: set the position / rotation / scale
        geom.run(vis)
    finally:
        vis.socket.emit("room:geometry:queue", -1)

    # wait and then disconnect
    vis.socket.sleep(1)
    vis.socket.disconnect()
