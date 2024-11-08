import logging
import typing as t
import urllib.request
from io import StringIO

import ase.io
import znsocket
from celery import shared_task
from flask import current_app

from zndraw.analyse import analyses
from zndraw.base import FileIO
from zndraw.bonds import ASEComputeBonds
from zndraw.draw import geometries
from zndraw.exceptions import RoomLockedError
from zndraw.modify import modifier
from zndraw.scene import Scene
from zndraw.selection import selections
from zndraw.utils import load_plots_to_dict

log = logging.getLogger(__name__)


def _get_default_generator(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    return [ase.Atoms()]


def _get_zntrack_generator(file_io: FileIO) -> t.Iterable[ase.Atoms]:
    node_name, attribute = file_io.name.split(".", 1)
    try:
        import os

        ## FIX for zntrack bug https://github.com/zincware/ZnTrack/issues/806
        import sys

        import zntrack

        sys.path.insert(0, os.getcwd())
        ##

        node = zntrack.from_rev(node_name, remote=file_io.remote, rev=file_io.rev)
        images = node
        for key in attribute.split("."):
            if isinstance(images, list):
                images = images[int(key)]
                if not isinstance(images, list):
                    images = [images]
            else:
                images = getattr(images, key)

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


def get_generator_from_filename(
    file_io: FileIO, bond_calculator: ASEComputeBonds | None = None
) -> t.Iterable[ase.Atoms]:
    if file_io.name is None:
        gen = _get_default_generator(file_io)
    elif file_io.remote is not None:
        gen = _get_zntrack_generator(file_io)
    elif file_io.name.endswith((".h5", ".hdf5", ".h5md")):
        gen = _get_znh5md_generator(file_io)
    elif file_io.name.startswith(("http", "https")):
        gen = _get_http_generator(file_io)
    else:
        gen = _get_ase_generator(file_io)

    if bond_calculator is not None:
        for atoms in gen:
            if not hasattr(atoms, "connectivity"):
                atoms.connectivity = bond_calculator.get_bonds(atoms)
            yield atoms
    else:
        yield from gen


@shared_task
def read_file(fileio: dict) -> None:
    file_io = FileIO(**fileio)

    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token="default",
    )
    bonds_calculator = ASEComputeBonds()
    if current_app.config.get("COMPUTE_BONDS", False):
        generator = get_generator_from_filename(file_io, bonds_calculator)
    else:
        generator = get_generator_from_filename(file_io)

    # TODO: vis.extend(generator) # vis does not yet support consuming a generator
    atoms_buffer = []
    for atoms in generator:
        atoms_buffer.append(atoms)
        if len(atoms_buffer) > 10:
            vis.extend(atoms_buffer)
            atoms_buffer = []
    vis.extend(atoms_buffer)

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


@shared_task
def read_plots(paths: list[str], remote, rev) -> None:
    from zndraw.zndraw import ZnDrawLocal

    vis = ZnDrawLocal(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token="default",
    )

    vis.figures = load_plots_to_dict(paths, remote, rev)

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def load_zntrack_frames(room: str, remote: str, rev: str, attribute: str, name: str):
    if len(rev) == 0:
        rev = None
    import zntrack

    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )
    try:
        vis.extend(getattr(zntrack.from_rev(name, remote=remote, rev=rev), attribute))
    except Exception as e:
        vis.log(str(e))

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def load_zntrack_figures(room: str, remote: str, rev: str, attribute: str, name: str):
    if len(rev) == 0:
        rev = None
    import zntrack

    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )
    try:
        vis.figures.update(
            getattr(zntrack.from_rev(name, remote=remote, rev=rev), attribute)
        )
    except Exception as e:
        vis.log(str(e))

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def inspect_zntrack_node(name, rev, remote):
    import dataclasses

    import zntrack
    from zntrack.config import ZNTRACK_OPTION, ZnTrackOptionEnum

    node = zntrack.from_rev(name=name, remote=remote, rev=rev)

    # find all @property decorated methods
    def find_properties(cls):
        exclude = set(
            [
                "_init_descriptors_",
                "_init_subclass_basecls_",
                "_use_repr_",
                "nwd",
                "state",
                "uuid",
            ]
        )
        names = []
        for name in dir(cls):
            if name not in exclude:
                attr = getattr(cls, name)
                if isinstance(attr, property):
                    # get the type hint of the property
                    hint = attr.fget.__annotations__.get("return", None)
                    names.append((name, str(hint)))
        return names

    def find_deps_outs(cls):
        names = []
        for field in dataclasses.fields(cls):
            if field.metadata.get(ZNTRACK_OPTION) in [
                ZnTrackOptionEnum.OUTS,
                ZnTrackOptionEnum.DEPS,
            ]:
                # names.append(field.name)
                # append a tuple with the name and the type
                names.append((field.name, str(field.type)))
        return names

    return find_deps_outs(node) + find_properties(node.__class__)


@shared_task
def run_geometry_schema(room) -> None:
    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    dct = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"schema:{room}:geometry",
    )
    for key, val in geometries.items():
        dct[key] = val.model_json_schema()
    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_scene_schema(room) -> None:
    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    dct = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"schema:{room}:scene",
    )

    scene_schema = Scene.model_json_schema_from_atoms(vis.atoms)

    # we also want to initialize the `vis.config`
    # calling the config will set the default values
    # but not overwrite the existing ones, if they are set
    orig_config = vis.config

    for key, val in scene_schema["properties"].items():
        try:
            scene_schema["properties"][key]["default"] = orig_config["scene"][key]
        except KeyError:
            vis.log(f"KeyError: {key}")

    dct["Scene"] = scene_schema

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_analysis_schema(room) -> None:
    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    dct = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"schema:{room}:analysis",
    )
    for key, val in analyses.items():
        dct[key] = val.model_json_schema_from_atoms(vis.atoms)
    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_modify_schema(room) -> None:
    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    dct = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key="schema:default:modifier",
    )
    for key, val in modifier.items():
        dct[key] = val.model_json_schema()

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def setup_public_modifier(room):
    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    room_modifier_queue = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"queue:{room}:modifier",
    )

    # We link the room queue to the default queue
    # so that the public modifier can be used in the room
    # public modifiers iterate the default_modifier_queue
    # and check for applicable tasks to pick up.
    # TODO: remove rooms from this list if they are closed
    # TODO: when to close a room?
    default_modifier_queue = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key="queue:default:modifier",
    )

    default_modifier_queue[room] = room_modifier_queue

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_selection_schema(room) -> None:
    from zndraw import ZnDraw
    from zndraw.selection import selections

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    dct = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"schema:{room}:selection",
    )
    for key, val in selections.items():
        dct[key] = val.model_json_schema()
    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_room_worker(room):
    from zndraw import ZnDraw

    vis = ZnDraw(
        r=current_app.extensions["redis"],
        url=current_app.config["SERVER_URL"],
        token=room,
    )

    geometry_queue = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"queue:{room}:geometry",
    )
    for key in geometry_queue:
        if key in geometries:
            try:
                task = geometry_queue.pop(key)
                vis.geometries.append(geometries[key](**task))
            except IndexError:
                pass

    selection_queue = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"queue:{room}:selection",
    )

    for key in selection_queue:
        if key in selections:
            try:
                task = selection_queue.pop(key)
                selections[key](**task).run(vis)
            except IndexError:
                pass

    analysis_queue = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"queue:{room}:analysis",
    )

    for key in analysis_queue:
        if key in analyses:
            try:
                task = analysis_queue.pop(key)
                analyses[key](**task).run(vis)
            except IndexError:
                pass

    scene_queue = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"queue:{room}:scene",
    )

    for key, val in scene_queue.items():
        if key == "Scene":  # hotfix
            vis.config["scene"] = val
        else:
            vis.config[key] = val

    modifier_queue = znsocket.Dict(
        r=current_app.extensions["redis"],
        socket=vis._refresh_client,
        key=f"queue:{room}:modifier",
    )

    for key in modifier_queue:
        if key in modifier:
            try:
                task = modifier_queue.pop(key)
                modifier[key](**task).run(vis)
            except IndexError:
                pass

    # wait and then disconnect
    vis.socket.sleep(1)
    vis.socket.disconnect()
