import datetime
import logging
import pathlib
from dataclasses import asdict

import ase.io
import znh5md
from celery import chain, shared_task
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
from ..data import CeleryTaskData, RoomGetData, RoomSetData
from ..utils import typecast
from .utils import get_queue_position, insert_into_queue, update_job_status

log = logging.getLogger(__name__)


def get_client(url) -> Client:
    client = Client()
    client.connect(url, wait_timeout=10)
    return client


@shared_task
def update_atoms(token: str, data: list) -> None:
    """Update the atoms in the database.

    Attributes
    ----------
    token : str
        The token of the room.
    index : int
        The index of the frame.
    data : dict
        The data of the frame.
    """
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        # check if the index is already in the db
        for new_frame in data:
            index = new_frame["index"]
            frame = ses.query(db_schema.Frame).filter_by(index=index, room=room).first()
            if frame is not None:
                # if so, update the data
                frame.data = new_frame["data"]
            else:
                # create a db_schema.Frame
                frame = db_schema.Frame(index=index, data=new_frame["data"], room=room)
                ses.add(frame)
        ses.commit()


@shared_task
def get_selection_schema(url: str, target: str):
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
    con.emit("celery:task:emit", asdict(msg))


@shared_task
def modifier_schema(url: str, token: str):
    config = GlobalConfig.load()
    include = []

    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        room_modifiers = ses.query(db_schema.Modifier).filter_by(room=room).all()
        modifiers = ses.query(db_schema.Modifier).filter_by(room_token=None).all()

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
    con.emit("celery:task:emit", asdict(msg))


@shared_task
def read_file(url: str, target: str, token: str):
    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDrawWorker(token=token, url=url)

    if len(vis) == 0:
        fileio = cache.get("FILEIO")
        if fileio is None or fileio.name is None:
            vis.append(ase.Atoms())
        elif fileio.remote is not None:
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

        atoms_list = []

        for idx, atoms in enumerate(generator):
            if fileio.start and idx < fileio.start:
                continue
            if fileio.stop and idx >= fileio.stop:
                break
            if fileio.step and idx % fileio.step != 0:
                continue
            atoms_list.append(atoms)

        vis.extend(atoms_list)
        vis.step = len(vis) - 1
    else:
        vis.upload(target)

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_selection(url: str, token: str, data: dict):
    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDrawWorker(token=str(token), url=url)
    print(datetime.datetime.now().isoformat())

    config = GlobalConfig.load()
    cls = get_selection_class(config.get_selection_methods())

    try:
        selection = cls(**data)
        selection.run(vis)
    except ValueError as err:
        vis.log(str(err))

    print(datetime.datetime.now().isoformat())

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def run_analysis(url: str, token: str, data: dict):
    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDraw(url=url, token=token)
    vis = ZnDrawWorker(token=str(token), url=url)

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
def update_queue_positions(queue_name, url):
    from zndraw.zndraw_worker import ZnDrawWorker

    if queue_name == "slow":
        queue_positions = get_queue_position(queue_name)
        worker = ZnDrawWorker(token="None", url=url)
        for position, room_token in queue_positions:
            msg = CeleryTaskData(
                target=f"webclients_{room_token}",
                event="modifier:queue:update",
                data=position,
                disconnect=False,
            )
            worker.socket.emit("celery:task:emit", msg.to_dict())
        worker.socket.sleep(0.5)
        worker.socket.disconnect()
    else:
        return None


@shared_task(bind=True)
def _run_global_modifier(self, url: str, token: str, data, queue_job_id: str):
    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDrawWorker(token=str(token), url=url)
    vis.socket.on("modifier:run:finished", lambda: vis.socket.disconnect())

    name = data["method"]["discriminator"]
    while True:
        with Session() as ses:
            # get the available hosts for the modifier
            modifier = ses.query(db_schema.Modifier).filter_by(name=name, room_token=None).first()
            host = (
                ses.query(db_schema.ModifierClient)
                .filter_by(modifier=modifier, available=True)
                .first()
            )
            assigned_hosts = (
                ses.query(db_schema.ModifierClient)
                .filter_by(modifier=modifier)
                .count()
            )

        if assigned_hosts == 0:
            update_job_status(job_id=queue_job_id, status="failed:no_host")
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
        to_get = RoomGetData.get_current_state()
        cache = vis.get_properties(**to_get.to_dict())
        msg = CeleryTaskData(
            target=host.sid,
            event="modifier:run",
            data={"params": data, "token": vis.token, "cache": cache},
        )
        vis.socket.emit("celery:task:emit", asdict(msg))
        update_job_status(job_id=queue_job_id, status="running")
        # add additional 5 seconds for communication overhead
        for _ in range(int(host.timeout + 5)):
            if vis.socket.connected:
                # We have seen that socket.connected is not always reliable
                # if there are timeout issues, check here
                vis.socket.sleep(1)
            else:
                log.critical("Modifier finished")
                status = "finished"
                update_job_status(job_id=queue_job_id, status=status)
                update_queue_positions(queue_name="slow", url=url)
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
        update_job_status(job_id=queue_job_id, status="failed:timeout")
        update_queue_positions(queue_name="slow", url=url)
        return


@shared_task(bind=True)
def _run_room_modifier(self, url: str, token: str, data, queue_job_id: str):
    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDrawWorker(token=str(token), url=url)

    name = data["method"]["discriminator"]
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=vis.token).first()
        room_modifier = (
            ses.query(db_schema.Modifier).filter_by(room=room, name=name).first()
        )
        hosts = (
            ses.query(db_schema.ModifierClient)
            .filter_by(modifier=room_modifier)
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


@shared_task(bind=True)
def _run_default_modifier(self, url: str, token: str, data: dict, queue_job_id: str):
    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDrawWorker(token=str(token), url=url)
    config = GlobalConfig.load()
    cls = get_modify_class(config.get_modify_methods(include_private=True))
    modifier = cls(**data)

    msg = CeleryTaskData(
        target=f"webclients_{vis.token}", event="modifier:run:running", data=None
    )

    vis.socket.emit("celery:task:emit", asdict(msg))

    try:
        modifier.run(vis)
        status = "finished"
    except Exception as err:
        # exception type
        exception_type = type(err).__name__
        vis.log(f"Error: {exception_type}:{err}")
        status = f"failed: {exception_type}:{err}"

    msg = CeleryTaskData(
        target=f"webclients_{vis.token}",
        event="modifier:run:finished",
        data=None,
        disconnect=True,
    )

    vis.socket.emit("celery:task:emit", asdict(msg))
    update_job_status(job_id=queue_job_id, status=status)


def route_modifier_to_queue(name: str, token: str) -> str:
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        room_modifiers = ses.query(db_schema.Modifier).filter_by(room=room).all()
        modifiers = ses.query(db_schema.Modifier).filter_by(room_token=None).all()
        custom_global_modifiers = [modifier.name for modifier in modifiers]
        custom_room_modifiers = [modifier.name for modifier in room_modifiers]
    if name in custom_global_modifiers:
        queue_name = "slow"
    elif name in custom_room_modifiers:
        queue_name = "custom"
    else:
        queue_name = "default"
    return queue_name


def run_modifier(url: str, token: str, data: dict):
    name = data["method"]["discriminator"]
    queue_name = route_modifier_to_queue(name, token)
    queue_job_id = insert_into_queue(
        queue_name=queue_name, job_name=name, room_token=str(token)
    )
    if queue_name == "slow":
        task_chain = chain(
            update_queue_positions.si(queue_name, url),
            _run_global_modifier.si(url, token, data, queue_job_id),
        )
    elif queue_name == "custom":
        task_chain = chain(
            update_queue_positions.si(queue_name, url),
            _run_room_modifier.si(url, token, data, queue_job_id),
        )
    else:
        task_chain = chain(
            update_queue_positions.si(queue_name, url),
            _run_default_modifier.si(url, token, data, queue_job_id),
        )
    task_chain.delay()


@shared_task
@typecast
def handle_room_get(data: RoomGetData, token: str, url: str, target: str):
    from zndraw.zndraw_worker import ZnDrawWorker

    worker = ZnDrawWorker(token=token, url=url)
    answer = worker.get_properties(**data.to_dict())
    msg = CeleryTaskData(target=target, event="room:get", data=answer, disconnect=True)
    worker.socket.emit("celery:task:emit", msg.to_dict())


@shared_task
@typecast
def handle_room_set(data: RoomSetData, token: str, url: str, source: str):
    from zndraw.zndraw_worker import ZnDrawWorker

    worker = ZnDrawWorker(token=token, url=url, emit=False)
    worker.set_properties(**data.to_dict())

    msg = CeleryTaskData(
        target=source,
        event="room:set:finished",
        data=None,
        disconnect=True,
    )
    worker.socket.emit("celery:task:emit", msg.to_dict())
    # worker.commit() and a mode, that waits for all updates before committing


@shared_task
def activate_modifier(sid: str, available: bool):
    with Session() as ses:
        global_modifier_client = (
            ses.query(db_schema.ModifierClient).filter_by(sid=sid).all()
        )
        room_modifier_client = (
            ses.query(db_schema.ModifierClient).filter_by(sid=sid).all()
        )
        for gmc in global_modifier_client:
            gmc.available = available
        for rmc in room_modifier_client:
            rmc.available = available
        ses.commit()


@shared_task
def on_disconnect(token: str, sid: str, url: str):
    # from zndraw.zndraw_worker import ZnDrawWorker

    # worker = ZnDrawWorker(token=token, url=url)
    with Session() as ses:
        room = ses.query(db_schema.Room).filter_by(token=token).first()
        client = ses.query(db_schema.WebClient).filter_by(sid=sid).first()
        if client is not None:
            ses.delete(client)
            if client.host:
                new_host = ses.query(db_schema.WebClient).filter_by(room=room).first()
                if new_host is not None:
                    new_host.host = True
            ses.commit()

    with Session() as ses:
        global_modifier_client = (
            ses.query(db_schema.ModifierClient).filter_by(sid=sid).all()
        )
        for gmc in global_modifier_client:
            ses.delete(gmc)

        room_modifier_client = (
            ses.query(db_schema.ModifierClient).filter_by(sid=sid).all()
        )
        for rmc in room_modifier_client:
            ses.delete(rmc)
        ses.commit()

    # with Session() as ses:
    # room = ses.query(db_schema.Room).filter_by(token=token).first()
    # clients = ses.query(db_schema.Client).filter_by(room=room).all()
    # connected_users = [{"name": client.name} for client in clients]

    # msg = CeleryTaskData(
    #     target=f"webclients_{token}",
    #     event="connectedUsers",
    #     data=list(reversed(connected_users)),
    #     disconnect=True,
    # )
    # # TODO this can not work, because it triggers itself

    # worker.socket.emit("celery:task:emit", msg.to_dict())


@shared_task
def upload_file(url: str, token: str, filename: str, content: str):
    from io import StringIO

    import ase.io

    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDrawWorker(token=token, url=url)
    vis.log(f"Uploading {filename}")

    format = filename.split(".")[-1]
    format = format if format != "xyz" else "extxyz"

    if format == "h5":
        raise ValueError("H5MD format not supported for uploading yet")

    stream = StringIO(content)

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

    vis.socket.sleep(1)
    vis.socket.disconnect()


@shared_task
def download_file(url: str, token: str, sid: str):
    from io import StringIO

    import ase.io

    from zndraw.zndraw_worker import ZnDrawWorker

    vis = ZnDrawWorker(token=token, url=url)
    atoms = vis.atoms
    if len(vis.selection) != 0:
        atoms = atoms[vis.selection]
    file = StringIO()
    ase.io.write(file, atoms, format="extxyz")
    file.seek(0)
    msg = CeleryTaskData(
        target=sid,
        event="file:download",
        data=file.read(),
        disconnect=True,
    )
    vis.socket.emit("celery:task:emit", asdict(msg))
