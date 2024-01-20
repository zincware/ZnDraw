import contextlib
import dataclasses
import datetime
import logging
import pathlib
import threading
import typing as t
import uuid
from io import StringIO

import ase
import ase.io
import numpy as np
import socketio
import tqdm
import znh5md
from znframe.frame import Frame

from zndraw.analyse import get_analysis_class
from zndraw.draw import Geometry
from zndraw.modify import UpdateScene, get_modify_class
from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig
from zndraw.utils import (
    ZnDrawLoggingHandler,
    get_cls_from_json_schema,
    hide_discriminator_field,
)

log = logging.getLogger(__name__)


@dataclasses.dataclass
class Config:
    """Configuration for ZnDraw client.

    Attributes
    ----------
    call_timeout : int
        Timeout for socket calls in seconds.
        Set to a smaller value to fail faster.
    """

    call_timeout: int = 60


@dataclasses.dataclass
class FileIO:
    name: str = None
    start: int = 0
    stop: int = None
    step: int = 1
    remote: str = None
    rev: str = None


@dataclasses.dataclass
class ZnDrawBase:  # collections.abc.MutableSequence
    """

    Attributes
    ----------
    display_new : bool
        Display new atoms in the webclient, when they are added.

    token : str
        Identifies the session this instances is being connected to.
        Tokens can be shared.
    _uuid : uuid.UUID
        Unique identifier for this instance. Can be set for reconnecting
        but only ONE instance with the same uuid can be connected at the same time.
    auth_token : str
        Authentication token, used e.g. for registering modifiers to all users and
        not just the current session.
    """

    url: str
    token: str = "notoken"
    display_new: bool = True
    _uuid: uuid.UUID = dataclasses.field(default_factory=uuid.uuid4)
    auth_token: str = None
    config: Config = dataclasses.field(default_factory=Config)

    _target_sid: str = None

    def __post_init__(self):
        self._uuid = str(self._uuid)
        self.socket = socketio.Client()
        self.socket.on(
            "connect",
            lambda: self.socket.emit(
                "join",
                {
                    "token": self.token,
                    "uuid": self._uuid,
                    "auth_token": self.auth_token,
                },
            ),
        )
        self.socket.on("disconnect", lambda: self.socket.disconnect())
        self.socket.on("modifier:run", self._pre_modifier_run)
        self.socket.on("message:log", lambda data: print(data))

    def _connect(self):
        for _ in range(100):
            try:
                self.socket.connect(self.url)
                break
            except socketio.exceptions.ConnectionError:
                self.socket.sleep(0.1)
        else:
            raise socketio.exceptions.ConnectionError

    def reconnect(self) -> None:
        """Reconnect to the server."""
        self.socket.disconnect()
        self._connect()

    @contextlib.contextmanager
    def _set_sid(self, sid):
        self._target_sid = sid
        yield
        self._target_sid = None

    @contextlib.contextmanager
    def _set_token(self, token):
        old_token = self.token
        self.token = token
        yield
        self.token = old_token

    def __len__(self) -> int:
        return int(
            self.socket.call(
                "atoms:length", {"token": self.token}, timeout=self.config.call_timeout
            )
        )

    def __setitem__(self, index, value):
        if not isinstance(value, ase.Atoms) and not isinstance(value, Frame):
            raise ValueError("Must be an ase.Atoms or Frame object")

        assert isinstance(index, int), "Index must be an integer"

        if isinstance(value, ase.Atoms):
            value = Frame.from_atoms(value)
        data = {
            index: value.to_dict(built_in_types=False),
            "display_new": self.display_new,
            "token": self.token,
        }
        if self._target_sid is not None:
            # this only affects initial loading of data if a new webclient connects
            data["sid"] = self._target_sid
        self.socket.emit("atoms:upload", data)

    def __delitem__(self, index):
        if (
            isinstance(index, int)
            or isinstance(index, slice)
            or isinstance(index, list)
        ):
            length = len(self)
            is_slice = isinstance(index, slice)
            if is_slice:
                index = range(*index.indices(length))

            index = [index] if isinstance(index, int) else index
            index = [i if i >= 0 else length + i for i in index]

            data = {"index": index, "token": self.token}

            self.socket.emit("atoms:delete", data)
            if not is_slice and (index[0] >= length or index[-1] >= length):
                raise IndexError("Index out of range")
        else:
            raise TypeError("Index must be an integer, slice or list[int]")

    def insert(self, index, value: t.Union[ase.Atoms, Frame]):
        """Insert atoms before index"""
        if isinstance(value, ase.Atoms):
            value = Frame.from_atoms(value)

        data = list(self)
        data.insert(index, value)
        for idx, val in enumerate(data):
            self[idx] = val

        # TODO: why is this not working at the moment?
        # self.socket.emit("atoms:insert", {index: value.to_dict()})

    def append(self, value: t.Union[ase.Atoms, Frame]) -> None:
        """Append atoms to the end of the list"""
        if isinstance(value, ase.Atoms):
            value = Frame.from_atoms(value)
        self[len(self)] = value

    def extend(
        self, values: t.Union[ase.Atoms, Frame, list[ase.Atoms], list[Frame]]
    ) -> None:
        """Extend the list by appending all the items in the given list"""
        size = len(self)
        for idx, value in enumerate(values):
            if isinstance(value, ase.Atoms):
                value = Frame.from_atoms(value)
            self[size + idx] = value

    def __getitem__(self, index) -> t.Union[ase.Atoms, list[ase.Atoms]]:
        length = len(self)
        is_scalar = isinstance(index, int)
        is_sclice = isinstance(index, slice)
        if is_sclice:
            index = range(*index.indices(length))

        index = [index] if isinstance(index, int) else index
        index = [i if i >= 0 else length + i for i in index]

        data = {"indices": index, "token": self.token}
        downloaded_data = self.socket.call(
            "atoms:download", data, timeout=self.config.call_timeout
        )

        atoms_list = []

        for val in downloaded_data.values():
            atoms_list.append(Frame.from_dict(val).to_atoms())

        data = atoms_list[0] if is_scalar else atoms_list
        if data == [] and not is_sclice:
            raise IndexError("Index out of range")
        return data

    def log(self, message: str) -> None:
        """Log a message to the console"""
        print(message)
        self.socket.emit(
            "message:log",
            {
                "message": message,
                "token": self.token,
            },
        )

    @property
    def atoms(self) -> ase.Atoms:
        """Return the atoms at the current step."""
        return self[self.step]

    @property
    def points(self) -> np.ndarray:
        data = self.socket.call(
            "points:get", {"token": self.token}, timeout=self.config.call_timeout
        )
        return np.array([[val["x"], val["y"], val["z"]] for val in data])

    @points.setter
    def points(self, value: t.Union[np.ndarray, list]) -> None:
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if len(value) > 0:
            try:
                assert len(value[0]) == 3
            except (TypeError, AssertionError):
                raise ValueError("Points must be a list of 3D coordinates")
        data = {"value": value, "token": self.token}
        self.socket.emit("points:set", data)

    @property
    def segments(self) -> np.ndarray:
        data = self.socket.call(
            "scene:segments", {"token": self.token}, timeout=self.config.call_timeout
        )
        return np.array(data)

    @property
    def step(self) -> int:
        step = int(
            self.socket.call(
                "scene:step",
                {"token": self.token},
                timeout=self.config.call_timeout,
            )
        )
        return step

    @step.setter
    def step(self, index):
        if index > len(self) - 1:
            raise IndexError(f"Index {index} out of range for length {len(self)}")
        data = {"index": index, "token": self.token}
        self.socket.emit("scene:set", data)

    @property
    def selection(self) -> list[int]:
        return self.socket.call(
            "selection:get", {"token": self.token}, timeout=self.config.call_timeout
        )

    @selection.setter
    def selection(self, value: list[int]):
        data = {"selection": value, "token": self.token}
        self.socket.emit("selection:set", data)

    def play(self):
        self.socket.emit("scene:play", {"token": self.token})

    def pause(self):
        self.socket.emit("scene:pause", {"token": self.token})

    @property
    def figure(self):
        raise NotImplementedError("Gathering figure from webclient not implemented yet")

    @figure.setter
    def figure(self, fig: str):
        data = {"figure": fig, "token": self.token}
        self.socket.emit("analysis:figure", data)

    @property
    def bookmarks(self) -> dict:
        return self.socket.call(
            "bookmarks:get", {"token": self.token}, timeout=self.config.call_timeout
        )

    @bookmarks.setter
    def bookmarks(self, value: dict):
        data = {"bookmarks": value, "token": self.token}
        self.socket.emit("bookmarks:set", data)

    def _pre_modifier_run(self, data) -> None:
        # TODO: this should be inside the _set_sid ?
        self.socket.emit(
            "modifier:run:running",
            {
                "token": data["target"],
                "name": data["params"]["method"]["discriminator"],
            },
        )
        try:
            self._modifier_run(data)
        except Exception as err:
            self.log(f"Modifier failed with error: {repr(err)}")
        self.socket.emit(
            "modifier:run:finished",
            {
                "token": data["target"],
            },
        )

    def _modifier_run(self, data) -> None:
        raise NotImplementedError


@dataclasses.dataclass
class ZnDrawDefault(ZnDrawBase):
    file_io: FileIO = None

    def __post_init__(self):
        super().__post_init__()

        self.socket.on("webclient:available", self.initialize_webclient)
        self.socket.on("selection:run", self.selection_run)
        self.socket.on("analysis:run", self.analysis_run)
        self.socket.on("upload", self.upload_file)
        self.socket.on("download:request", self.download_file)
        self.socket.on("modifier:register", self.register_modifier)
        self.socket.on("scene:trash", self.trash_scene)
        self._connect()
        self.socket.wait()

    def initialize_webclient(self, sid):
        start_time = datetime.datetime.now()
        with self._set_sid(sid):
            for idx, atoms in enumerate(self.read_data()):
                if idx == 0:
                    self.analysis_schema(atoms)
                    self.selection_schema()
                    self.draw_schema()
                self[idx] = atoms
                # self.step = idx # double the message count ..., replace with part of the setitem message, benchmark
        log.warning(f"{datetime.datetime.now() - start_time} Finished sending data.")

    def read_data(self) -> t.Generator[ase.Atoms, None, None]:
        if self.file_io is None or self.file_io.name is None:
            yield ase.Atoms()
            return

        if self.file_io.remote is not None:
            node_name, attribute = self.file_io.name.split(".", 1)
            try:
                import zntrack

                node = zntrack.from_rev(
                    node_name, remote=self.file_io.remote, rev=self.file_io.rev
                )
                generator = getattr(node, attribute)
            except ImportError as err:
                raise ImportError(
                    "You need to install ZnTrack to use the remote feature"
                ) from err
        elif pathlib.Path(self.file_io.name).suffix == ".h5":
            reader = znh5md.ASEH5MD(self.file_io.name)
            generator = reader.get_atoms_list()
        else:
            generator = ase.io.iread(self.file_io.name)

        frame = 0

        for idx, atoms in tqdm.tqdm(enumerate(generator), ncols=100):
            if self.file_io.start and idx < self.file_io.start:
                continue
            if self.file_io.stop and idx >= self.file_io.stop:
                break
            if self.file_io.step and idx % self.file_io.step != 0:
                continue
            yield atoms
            frame += 1

    def analysis_schema(self, atoms: ase.Atoms):
        config = GlobalConfig.load()

        cls = get_analysis_class(config.get_analysis_methods())

        schema = cls.model_json_schema_from_atoms(atoms)
        hide_discriminator_field(schema)

        self.socket.emit(
            "analysis:schema",
            {
                "schema": schema,
                "sid": self._target_sid,
            },
        )

    def selection_schema(self):
        config = GlobalConfig.load()
        cls = get_selection_class(config.get_selection_methods())

        self.socket.emit(
            "selection:schema",
            {
                "schema": cls.model_json_schema(),
                "sid": self._target_sid,
            },
        )

    def draw_schema(self):
        self.socket.emit(
            "draw:schema",
            {"schema": Geometry.updated_schema(), "sid": self._target_sid},
        )

    def _modifier_run(self, data):
        with self._set_token(data["target"]):
            config = GlobalConfig.load()
            cls = get_modify_class(config.get_modify_methods())
            modifier = cls(**data["params"])
            modifier.run(self)

    def selection_run(self, data):
        with self._set_token(data["target"]):
            config = GlobalConfig.load()
            cls = get_selection_class(config.get_selection_methods())

            try:
                selection = cls(**data["params"])
                selection.run(self)
            except ValueError as err:
                log.critical(err)

    def analysis_run(self, data):
        with self._set_token(data["target"]):
            config = GlobalConfig.load()
            cls = get_analysis_class(config.get_analysis_methods())

            try:
                instance = cls(**data["params"])
                instance.run(self)
            except ValueError as err:
                log.critical(err)

    def upload_file(self, data):
        with self._set_token(data["target"]):
            data = data["data"]

            format = data["filename"].split(".")[-1]

            if format == "h5":
                raise ValueError("H5MD format not supported for uploading yet")
            else:
                stream = StringIO(data["content"])
                del self[:]
                for idx, atoms in tqdm.tqdm(
                    enumerate(ase.io.iread(stream, format=format))
                ):
                    self.append(Frame.from_atoms(atoms))
                    self.step = idx

    def download_file(self, data):
        with self._set_sid(data["sid"]):
            atoms_list = list(self)
            if "selection" in data:
                atoms_list = [atoms_list[data["selection"]] for atoms in atoms_list]

            file = StringIO()
            ase.io.write(file, atoms_list, format="extxyz")
            file.seek(0)
            self.socket.emit(
                "download:response", {"data": file.read(), "sid": self._target_sid}
            )

    def register_modifier(self, data):
        include = [
            get_cls_from_json_schema(conf["schema"], conf["name"])
            for conf in data.get("modifiers", [])
        ]
        config = GlobalConfig.load()
        cls = get_modify_class(config.get_modify_methods(include=include))
        sid = data["token"]

        schema = cls.model_json_schema()

        hide_discriminator_field(schema)

        data = {"schema": schema, "token": sid}
        self.socket.emit("modifier:schema", data)

    def trash_scene(self, data):
        with self._set_token(data["target"]):
            # remove everything after the current step
            del self[self.step + 1 :]
            if len(self.selection) == 0:
                self.append(ase.Atoms())
            else:
                # remove the selected atoms
                atoms = self.atoms
                del atoms[self.selection]
                self.append(atoms)
            self.selection = []
            self.points = []


@dataclasses.dataclass
class ZnDraw(ZnDrawBase):
    """ZnDraw client."""

    url: str = None

    jupyter: bool = False

    _modifiers: dict = dataclasses.field(default_factory=dict)

    def __post_init__(self):
        super().__post_init__()
        if self.url is None:
            from zndraw.utils import get_port
            from zndraw.view import view

            port = get_port()
            self._view_thread = threading.Thread(
                target=view,
                kwargs={
                    "filename": None,
                    "port": port,
                    "open_browser": False,
                    "webview": False,
                    "fullscreen": False,
                    "start": 0,
                    "stop": None,
                    "step": 1,
                    "compute_bonds": True,
                },
            )
            self._view_thread.start()
            self.url = f"http://127.0.0.1:{port}"

        self._connect()

    def close(self):
        import time
        import urllib.request

        self.socket.disconnect()

        time.sleep(1)

        # open self.url/exit
        try:
            urllib.request.urlopen(f"{self.url}/exit")
        except Exception:
            pass
        if hasattr(self, "_view_thread"):
            log.debug("Waiting for ZnDraw client to close")
            # self._view_thread.terminate()
            self._view_thread.join()
            # raise ValueError("ZnDraw client closed")

    def _repr_html_(self):
        from IPython.display import IFrame

        return IFrame(src=self.url, width="100%", height="600px")._repr_html_()

    def get_logging_handler(self) -> ZnDrawLoggingHandler:
        return ZnDrawLoggingHandler(self)

    def register_modifier(
        self, cls: UpdateScene, run_kwargs: dict = None, default: bool = False
    ):
        """Register a modifier class.

        Attributes
        ----------
        cls : UpdateScene
            The modifier class to register.
        run_kwargs : dict, optional
            Keyword arguments to pass to the run method of the modifier class.
        default : bool, optional
            Whether to enable the modifier for ALL sessions of the ZnDraw client,
            or just the session for the given token.
        """
        if run_kwargs is None:
            run_kwargs = {}
        self.socket.emit(
            "modifier:register",
            {
                "uuid": self._uuid,
                "modifiers": [
                    {
                        "schema": cls.model_json_schema(),
                        "name": cls.__name__,
                        "default": default,
                    }
                ],
            },
        )
        self._modifiers[cls.__name__] = {"cls": cls, "run_kwargs": run_kwargs}

    def _modifier_run(self, data):
        # TODO: send back a response that the modifier has been received, otherwise log a warning
        # that the modifier has not been received to the user
        with self._set_token(data["target"]):
            config = GlobalConfig.load()
            cls = get_modify_class(
                config.get_modify_methods(
                    include=[x["cls"] for x in self._modifiers.values()]
                )
            )
            modifier = cls(**data["params"])
            modifier.run(
                self,
                **self._modifiers[modifier.method.__class__.__name__]["run_kwargs"],
            )
