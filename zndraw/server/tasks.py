import datetime
import pathlib
from dataclasses import asdict

import ase.io
import tqdm
import znframe
import znh5md
from celery import shared_task
from socketio import Client

from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig
from zndraw.zndraw import ZnDraw

from ..app import cache
from .data import CeleryTaskData, FrameData, SelectionRunData


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
        event="selection:schema",
        data=schema,
    )
    con = get_client(url)
    con.emit("celery:task:results", asdict(msg))


@shared_task
def scene_trash(url: str, token: str, sid: str):
    vis = ZnDraw(url=url, token=token, _target_sid=sid)
    del vis[vis.step + 1 :]
    print("AHOY!")
    if len(vis.selection) == 0:
        print("here")
        vis.append(ase.Atoms())
    else:
        # remove the selected atoms
        atoms = vis.atoms
        del atoms[vis.selection]
        vis.append(atoms)
    vis.selection = []
    vis.points = []


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
        )
        con.emit("celery:task:results", asdict(msg))
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
                index=idx,
                data=znframe.Frame.from_atoms(atoms).to_dict(built_in_types=False),
                update=True,
            ),
        )
        con.emit("celery:task:results", asdict(msg))

        frame += 1


@shared_task
def run_selection(url: str, token: str, data: SelectionRunData):
    data = SelectionRunData(**data)
    print(datetime.datetime.now().isoformat())
    vis = ZnDraw(url=url, token=token)

    config = GlobalConfig.load()
    cls = get_selection_class(config.get_selection_methods())

    try:
        selection = cls(**data.params)
        selection.run(vis)
    except ValueError as err:
        vis.log.critical(err)

    print(datetime.datetime.now().isoformat())
