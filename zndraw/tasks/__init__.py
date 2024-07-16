import logging
import typing as t
import urllib.request
from io import StringIO

import ase.io
import socketio.exceptions
import tqdm
import znjson
import znsocket
from celery import shared_task
from flask import current_app

from zndraw.base import FileIO
from zndraw.bonds import ASEComputeBonds
from zndraw.exceptions import RoomLockedError
from zndraw.utils import ASEConverter

log = logging.getLogger(__name__)


def _get_default_generator(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    return [ase.Atoms()]


def _get_zntrack_generator(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    node_name, attribute = file_io.name.split(".", 1)
    try:
        import zntrack

        node = zntrack.from_rev(node_name, remote=file_io.remote, rev=file_io.rev)
        images = getattr(node, attribute)
    except ImportError as err:
        raise ImportError(
            "You need to install ZnTrack to use the remote feature."
        ) from err

    return images[file_io.start : file_io.stop : file_io.step]


def _get_znh5md_generator(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    try:
        import znh5md

        io = znh5md.IO(file_io.name)
    except ImportError as err:
        raise ImportError(
            "You need to install ZnH5MD>=0.3 to use the remote feature."
        ) from err
    log.critical(file_io)
    return io[file_io.start : file_io.stop : file_io.step]


def _get_http_generator(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    format = file_io.name.split(".")[-1]
    format = format if format != "xyz" else "extxyz"
    content = urllib.request.urlopen(file_io.name).read().decode("utf-8")
    stream = StringIO(content)

    generator = ase.io.iread(
        stream, format=format, index=slice(file_io.start, file_io.stop, file_io.step)
    )

    return generator


def _get_ase_generator(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    generator = ase.io.iread(
        file_io.name, index=slice(file_io.start, file_io.stop, file_io.step)
    )

    return generator


def get_generator_from_filename(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    if file_io.name is None:
        return _get_default_generator(file_io)
    elif file_io.remote is not None:
        return _get_zntrack_generator(file_io)
    elif file_io.name.endswith((".h5", ".hdf5", ".h5md")):
        return _get_znh5md_generator(file_io)
    elif file_io.name.startswith(("http", "https")):
        return _get_http_generator(file_io)
    else:
        return _get_ase_generator(file_io)


@shared_task
def read_file(fileio: dict) -> None:
    file_io = FileIO(**fileio)
    r = current_app.extensions["redis"]

    io = socketio.Client()
    r.delete("room:default:frames")

    lst = znsocket.List(r, "room:default:frames")
    bonds_calculator = ASEComputeBonds()

    generator = get_generator_from_filename(file_io)

    for idx, atoms in tqdm.tqdm(enumerate(generator)):
        if current_app.config.get("COMPUTE_BONDS", False):
            if not hasattr(atoms, "connectivity"):
                atoms.connectivity = bonds_calculator.get_bonds(atoms)
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
def run_modifier(room, data: dict) -> None:
    from zndraw.modify import Modifier
    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )
    if not current_app.config.get("COMPUTE_BONDS", False):
        vis.bond_calculator = None
    vis.socket.emit("room:modifier:queue", 0)
    try:
        if vis.locked:
            raise RoomLockedError("The room you are trying to modify is locked.")
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
