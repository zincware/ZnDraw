import datetime
import logging
import pathlib
from dataclasses import asdict

import ase.io
import tqdm
import znframe
import znh5md
from celery import shared_task
from socketio import Client

from zndraw.analyse import get_analysis_class
from zndraw.db import Session
from zndraw.db import schema as db_schema
from zndraw.draw import Geometry
from zndraw.modify import get_modify_class
from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig
from zndraw.utils import get_cls_from_json_schema, hide_discriminator_field
from zndraw.zndraw import ZnDraw

from ..app import cache
from ..data import CeleryTaskData, FrameData

log = logging.getLogger(__name__)


def get_client(url) -> Client:
    client = Client()
    client.connect(url, wait_timeout=10)
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
        disconnect=True,
    )
    con = get_client(url)
    con.emit("celery:task:emit", asdict(msg))


@shared_task
def scene_schema(url: str, target: str):
    import enum

    from pydantic import BaseModel, Field

    class Material(str, enum.Enum):
        MeshBasicMaterial = "MeshBasicMaterial"
        MeshLambertMaterial = "MeshLambertMaterial"
        MeshMatcapMaterial = "MeshMatcapMaterial"
        MeshPhongMaterial = "MeshPhongMaterial"
        MeshPhysicalMaterial = "MeshPhysicalMaterial"
        MeshStandardMaterial = "MeshStandardMaterial"
        MeshToonMaterial = "MeshToonMaterial"

    # create a class for the material, resolution, etc.
    class Scene(BaseModel):
        fps: int = Field(60, ge=1, le=120, description="Maxiumm frames per second")
        material: Material = Field(Material.MeshPhongMaterial, description="Material")
        resolution: int = Field(10, ge=1, le=50, description="Resolution")
        particle_size: float = Field(1.0, ge=0.1, le=5, description="Particle Size")
        bonds_size: float = Field(1.0, ge=0.1, le=5, description="Bonds Size")
        wireframe: bool = Field(False, description="Wireframe")
        loop: bool = Field(
            False,
            alias="Animation Loop",
            description="Automatically restart animation when finished.",
        )
        simulation_box: bool = Field(
            False,
            description="Show the simulation box.",
        )
        bonds: bool = Field(
            True,
            description="Show bonds.",
        )
        line_label: bool = Field(
            True,
            description="Show the length of the line.",
        )
        label_offset: int = Field(
            0,
            ge=-7,
            le=7,
            description="Move the label to the left or right (keypress i).",
        )

    schema = Scene.model_json_schema()

    schema["properties"]["wireframe"]["format"] = "checkbox"
    schema["properties"]["Animation Loop"]["format"] = "checkbox"
    schema["properties"]["simulation_box"]["format"] = "checkbox"
    schema["properties"]["resolution"]["format"] = "range"
    schema["properties"]["label_offset"]["format"] = "range"
    schema["properties"]["particle_size"]["format"] = "range"
    schema["properties"]["fps"]["format"] = "range"
    schema["properties"]["particle_size"]["step"] = 0.1
    schema["properties"]["bonds_size"]["format"] = "range"
    schema["properties"]["bonds_size"]["step"] = 0.1
    schema["properties"]["bonds"]["format"] = "checkbox"
    schema["properties"]["line_label"]["format"] = "checkbox"

    msg = CeleryTaskData(
        target=target,
        event="scene:schema",
        data=schema,
        disconnect=True,
    )
    con = get_client(url)
    print(f"emitting scene_schema to {target}")
    con.emit("celery:task:emit", asdict(msg))


@shared_task
def geometries_schema(url: str, target: str):
    msg = CeleryTaskData(
        target=target,
        event="draw:schema",
        data=Geometry.updated_schema(),
        disconnect=True,
    )
    con = get_client(url)
    print(f"emitting scene_schema to {target}")
    con.emit("celery:task:emit", asdict(msg))


@shared_task
def analysis_schema(url: str, token: str):
    vis = ZnDraw(url=url, token=token)

    config = GlobalConfig.load()

    cls = get_analysis_class(config.get_analysis_methods())

    schema = cls.model_json_schema_from_atoms(vis[0])
    hide_discriminator_field(schema)

    msg = CeleryTaskData(
        target=f"webclients_{vis.token}",
        event="analysis:schema",
        data=schema,
        disconnect=True,
    )
    con = get_client(url)
    print(f"emitting analysis_schema to {vis.token}")
    con.emit("celery:task:emit", asdict(msg))
    vis.socket.disconnect()


@shared_task
def modifier_schema(url: str, token: str):
    config = GlobalConfig.load()
    include = []

    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        room_modifiers = ses.query(db_schema.RoomModifier).filter_by(room=room).all()
        modifiers = ses.query(db_schema.GlobalModifier).all()

        for modifier in modifiers:
            include.append(get_cls_from_json_schema(modifier.schema, modifier.name))

        for room_modifier in room_modifiers:
            include.append(
                get_cls_from_json_schema(room_modifier.schema, room_modifier.name)
            )

    cls = get_modify_class(
        config.get_modify_methods(include=include)
    )  # todo: include=include)
    schema = cls.model_json_schema()

    hide_discriminator_field(schema)
    msg = CeleryTaskData(
        target=f"webclients_{token}",
        event="modifier:schema",
        data=schema,
        disconnect=True,
    )
    con = get_client(url)
    print(f"emitting modifier_schema to {token}")
    con.emit("celery:task:emit", asdict(msg))


@shared_task
def scene_trash(url: str, token: str):
    vis = ZnDraw(url=url, token=token)
    del vis[vis.step + 1 :]
    if len(vis.selection) == 0:
        vis.append(ase.Atoms())
    else:
        # remove the selected atoms
        atoms = vis.atoms
        del atoms[vis.selection]
        if hasattr(atoms, "connectivity"):
            del atoms.connectivity
        vis.append(atoms)
    vis.selection = []
    vis.points = []

    vis.socket.sleep(10)
    vis.socket.disconnect()


@shared_task
def read_file(url: str, target: str, token: str):
    vis = ZnDraw(url=url, token=token)
    con = vis.socket

    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        # get all frames and iterate over them
        frames = ses.query(db_schema.Frame).filter_by(room=room)
        if frames.count() == 0:
            pass
        else:
            for frame in frames.all():
                msg = CeleryTaskData(
                    target=target,
                    event="atoms:upload",
                    data=FrameData(
                        index=frame.index,
                        data=frame.data,
                        update=True,
                    ),
                )
                con.emit("celery:task:emit", asdict(msg))
                con.sleep(0.01)
                # it is too fast for the socket to handle otherwise
            con.disconnect()
            return

    fileio = cache.get("FILEIO")

    if fileio.name is None:
        msg = CeleryTaskData(
            target=target,
            event="atoms:upload",
            data=FrameData(
                index=0,
                data=znframe.Frame.from_atoms(ase.Atoms()).to_dict(
                    built_in_types=False
                ),
                update=True,
            ),
            disconnect=True,
        )
        con.emit("celery:task:emit", asdict(msg))
        return

    if fileio.remote is not None:
        node_name, attribute = fileio.name.split(".", 1)
        try:
            import zntrack

            node = zntrack.from_rev(node_name, remote=fileio.remote, rev=fileio.rev)
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
                index=frame,
                data=znframe.Frame.from_atoms(atoms).to_dict(built_in_types=False),
                update=True,
            ),
        )
        con.emit("celery:task:emit", asdict(msg))

        frame += 1

    con.sleep(1)
    con.disconnect()


@shared_task
def run_selection(url: str, token: str, data: dict):
    print(datetime.datetime.now().isoformat())
    vis = ZnDraw(url=url, token=token)

    config = GlobalConfig.load()
    cls = get_selection_class(config.get_selection_methods())

    try:
        selection = cls(**data)
        selection.run(vis)
    except ValueError as err:
        vis.log.critical(err)

    print(datetime.datetime.now().isoformat())

    vis.socket.sleep(10)
    vis.socket.disconnect()


@shared_task
def run_analysis(url: str, token: str, data: dict):
    vis = ZnDraw(url=url, token=token)

    msg = CeleryTaskData(
        target=f"webclients_{vis.token}", event="analysis:run:running", data=None
    )

    vis.socket.emit("celery:task:emit", asdict(msg))

    config = GlobalConfig.load()
    cls = get_analysis_class(config.get_analysis_methods())

    analysis = cls(**data)
    try:
        analysis.run(vis)
    except Exception as err:
        vis.log(f"Error: {err}")

    msg = CeleryTaskData(
        target=f"webclients_{vis.token}",
        event="analysis:run:finished",
        data=None,
        disconnect=True,
    )

    vis.socket.emit("celery:task:emit", asdict(msg))


def run_global_modifier(vis, name, data):
    while True:
        with Session() as ses:
            # get the available hosts for the modifier
            modifier = ses.query(db_schema.GlobalModifier).filter_by(name=name).first()
            host = (
                ses.query(db_schema.GlobalModifierClient)
                .filter_by(global_modifier=modifier, available=True)
                .first()
            )
            assigned_hosts = (
                ses.query(db_schema.GlobalModifierClient)
                .filter_by(global_modifier=modifier)
                .count()
            )
        if assigned_hosts == 0:
            msg = CeleryTaskData(
                target=f"webclients_{vis.token}",
                event="modifier:run:finished",
                data=None,
            )
            vis.socket.emit("celery:task:emit", asdict(msg))
            msg = CeleryTaskData(
                target=f"webclients_{vis.token}",
                event="message:alert",
                data=f"Could not find any available modifier for {name}.",
                disconnect=True,
            )
            vis.socket.emit("celery:task:emit", asdict(msg))
            return

        if host is None:
            vis.socket.sleep(1)
            log.critical("No modifier available")
            continue

        # run the modifier
        msg = CeleryTaskData(
            target=host.sid,
            event="modifier:run",
            data={"params": data, "token": vis.token},
        )
        vis.socket.emit("celery:task:emit", asdict(msg))

        # add additional 5 seconds for communication overhead
        for _ in range(int(host.timeout + 5)):
            if vis.socket.connected:
                # We have seen that socket.connected is not always reliable
                # if there are timeout issues, check here
                vis.socket.sleep(1)
            else:
                log.critical("Modifier finished")
                return

        print("modifier timed out")
        msg = CeleryTaskData(
            target=f"webclients_{vis.token}",
            event="message:alert",
            data=f"Modifier {name} did not finish in time.",
        )
        vis.socket.emit("celery:task:emit", asdict(msg))

        msg = CeleryTaskData(
            target=f"webclients_{vis.token}",
            event="modifier:run:finished",
            data=None,
            disconnect=True,
        )
        vis.socket.emit("celery:task:emit", asdict(msg))
        return


def _run_room_modifier(vis, name, data):
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=vis.token).first()
        room_modifier = (
            ses.query(db_schema.RoomModifier).filter_by(room=room, name=name).first()
        )
        hosts = (
            ses.query(db_schema.RoomModifierClient)
            .filter_by(room_modifier=room_modifier)
            .first()
        )
        if hosts is None:
            msg = CeleryTaskData(
                target=f"webclients_{vis.token}",
                event="modifier:run:finished",
                data=None,
            )
            vis.socket.emit("celery:task:emit", asdict(msg))
            msg = CeleryTaskData(
                target=f"webclients_{vis.token}",
                event="message:alert",
                data=f"Could not find any available modifier for {name}.",
                disconnect=True,
            )
            vis.socket.emit("celery:task:emit", asdict(msg))
            return

        # run the modifier
        msg = CeleryTaskData(
            target=hosts.sid,
            event="modifier:run",
            data={"params": data, "token": vis.token},
            disconnect=True,
        )
        vis.socket.emit("celery:task:emit", asdict(msg))


@shared_task
def run_modifier(url: str, token: str, data: dict):
    vis = ZnDraw(url=url, token=token)

    def on_finished():
        # TODO: disconnect the modifier, release the worker
        vis.socket.disconnect()

    vis.socket.on("modifier:run:finished", on_finished)

    NAME = data["method"]["discriminator"]

    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        room_modifiers = ses.query(db_schema.RoomModifier).filter_by(room=room).all()
        modifiers = ses.query(db_schema.GlobalModifier).all()
        names = [modifier.name for modifier in modifiers]
        room_modifiers = [modifier.name for modifier in room_modifiers]

    if NAME in names:
        run_global_modifier(vis, NAME, data)
    elif NAME in room_modifiers:
        _run_room_modifier(vis, NAME, data)
    else:
        msg = CeleryTaskData(
            target=f"webclients_{vis.token}", event="modifier:run:running", data=None
        )

        vis.socket.emit("celery:task:emit", asdict(msg))

        config = GlobalConfig.load()
        cls = get_modify_class(config.get_modify_methods())

        modifier = cls(**data)

        try:
            modifier.run(vis)
        except Exception as err:
            vis.log(f"Error: {err}")

        msg = CeleryTaskData(
            target=f"webclients_{vis.token}",
            event="modifier:run:finished",
            data=None,
            disconnect=True,
        )

        vis.socket.emit("celery:task:emit", asdict(msg))
