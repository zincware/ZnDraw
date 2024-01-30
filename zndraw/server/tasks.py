import datetime
import pathlib
from dataclasses import asdict

import ase.io
import tqdm
import znframe
import znh5md
from celery import shared_task
from socketio import Client

from zndraw.analyse import get_analysis_class
from zndraw.draw import Geometry
from zndraw.modify import get_modify_class
from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig
from zndraw.utils import get_cls_from_json_schema, hide_discriminator_field
from zndraw.zndraw import ZnDraw

from ..app import cache
from ..data import CeleryTaskData, FrameData


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

    MODIFIER_SCHEMA = cache.get("MODIFIER_SCHEMA")
    ROOM_MODIFIER_SCHEMA = cache.get("ROOM_MODIFIER_SCHEMA").get(token, {})

    include = []
    for name, schema in MODIFIER_SCHEMA.items():
        include.append(get_cls_from_json_schema(schema, name))
    for name, schema in ROOM_MODIFIER_SCHEMA.items():
        include.append(get_cls_from_json_schema(schema, name))

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
def read_file(url: str, target: str):
    con = get_client(url)
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


@shared_task
def run_modifier(url: str, token: str, data: dict):
    vis = ZnDraw(url=url, token=token)
    vis.available = False

    def on_finished():
        vis.available = True

    vis.socket.on("modifier:run:finished", on_finished)

    NAME = data["method"]["discriminator"]

    MODIFIER_HOSTS = cache.get("MODIFIER_HOSTS")
    ROOM_MODIFIER_HOSTS = cache.get("ROOM_MODIFIER_HOSTS")

    if NAME in MODIFIER_HOSTS:
        while len(_hosts := cache.get("MODIFIER_HOSTS")[NAME]) > 0:
            _available = cache.get("MODIFIER_AVAILABLE")
            for pyclient in _hosts:
                if _available.get(pyclient, False):
                    # run the modifier
                    msg = CeleryTaskData(
                        target=pyclient,
                        event="modifier:run",
                        data={"params": data, "token": token},
                    )
                    vis.socket.emit("celery:task:emit", asdict(msg))
                    _timeout = cache.get("MODIFIER_TIMEOUT").get(pyclient, 0)
                    vis.socket.sleep(_timeout)
                    _available = cache.get("MODIFIER_AVAILABLE")
                    if not _available.get(pyclient, False):
                        print("modifier timed out")
                        msg = CeleryTaskData(
                            target=f"webclients_{vis.token}",
                            event="message:alert",
                            data=f"Modifier {NAME} did not finish in time.",
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
            else:
                vis.socket.sleep(1)
                print("No modifier available")

        msg = CeleryTaskData(
            target=f"webclients_{vis.token}",
            event="modifier:run:finished",
            data=None,
        )
        vis.socket.emit("celery:task:emit", asdict(msg))
        msg = CeleryTaskData(
            target=f"webclients_{vis.token}",
            event="message:alert",
            data=f"Could not find any available modifier for {NAME}.",
            disconnect=True,
        )
        vis.socket.emit("celery:task:emit", asdict(msg))
    elif NAME in ROOM_MODIFIER_HOSTS.get(vis.token, []):
        _available = cache.get("MODIFIER_AVAILABLE")
        for pyclient in ROOM_MODIFIER_HOSTS[vis.token][NAME]:
            if _available.get(pyclient, False):
                # run the modifier
                msg = CeleryTaskData(
                    target=pyclient,
                    event="modifier:run",
                    data={"params": data, "token": token},
                    disconnect=True,
                )
                vis.socket.emit("celery:task:emit", asdict(msg))
                # TODO: Do we want to support the timeout here?
        else:
            msg = CeleryTaskData(
                target=f"webclients_{vis.token}",
                event="modifier:run:finished",
                data=None,
            )
            vis.socket.emit("celery:task:emit", asdict(msg))

            msg = CeleryTaskData(
                target=f"webclients_{vis.token}",
                event="message:alert",
                data=f"All modifiers registerd for {NAME} are busy - Please try again.",
                disconnect=True,
            )
            vis.socket.emit("celery:task:emit", asdict(msg))
            return

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
