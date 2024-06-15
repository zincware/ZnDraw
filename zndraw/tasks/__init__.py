import logging
import typing as t

import ase.io
import socketio.exceptions
import tqdm
import znjson
import znsocket
from celery import shared_task
from flask import current_app

from zndraw.base import FileIO
from zndraw.bonds import ASEComputeBonds
from zndraw.utils import ASEConverter

log = logging.getLogger(__name__)


@shared_task
def run_znsocket_server(port: int) -> None:
    # Does not work with eventlet enabled!
    znsocket.Server(port=port).run()
    log.critical("ZnSocket server closed.")


@shared_task
def read_file(fileio: dict) -> None:
    file_io = FileIO(**fileio)
    # r = Redis(host="localhost", port=6379, db=0, decode_responses=True)
    r = current_app.extensions["redis"]

    io = socketio.Client()

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
        lst.append(
            znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter]))
        )
        if idx == 0:
            try:
                io.connect(current_app.config["SERVER_URL"], wait_timeout=10)
                io.emit("room:all:frames:refresh", [0])
            except socketio.exceptions.ConnectionError:
                pass

    while True:
        try:
            if not io.connected:
                io.connect(current_app.config["SERVER_URL"], wait_timeout=10)

            # updates len after all frames are loaded
            io.emit("room:all:frames:refresh", [idx])
            break
        except socketio.exceptions.ConnectionError:
            pass

    io.sleep(1)
    io.disconnect()


@shared_task
def compute_bonds(room=None) -> None:
    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token="default" if room is None else room,
    )

    bonds_calculator = ASEComputeBonds()
    for idx, atoms in enumerate(vis):
        try:
            atoms.connectivity = bonds_calculator.get_bonds(atoms)
            vis[idx] = atoms
        except Exception as e:
            vis.log(str(e))


@shared_task
def run_modifier(room, data: dict) -> None:
    from zndraw.modify import Modifier
    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )
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
def run_selection(room, data: dict) -> None:
    from zndraw.selection import Selection
    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )
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
def run_analysis(room, data: dict) -> None:
    from zndraw.analyse import Analysis
    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )
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
def run_geometry(room, data: dict) -> None:
    from zndraw.draw import Geometry
    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )
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


@shared_task
def run_upload_file(room, data: dict):
    from io import StringIO

    import ase.io

    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    vis.log(f"Uploading {data['filename']}")

    format = data["filename"].split(".")[-1]
    format = format if format != "xyz" else "extxyz"

    if format == "h5":
        raise ValueError("H5MD format not supported for uploading yet")

    stream = StringIO(bytes(data["content"]).decode("utf-8"))

    atoms_list = list(ase.io.iread(stream, format=format))
    if len(atoms_list) == 1 and len(vis.points) != 0:
        scene = vis.atoms
        atoms = atoms_list[0]
        if hasattr(scene, "connectivity"):
            del scene.connectivity
        for point in vis.points:
            atoms.positions -= atoms.get_center_of_mass() - point
            scene += atoms
        vis.append(scene)
    else:
        vis.extend(atoms_list)

    vis.step = len(vis) - 1

    vis.socket.sleep(1)
    vis.socket.disconnect()
