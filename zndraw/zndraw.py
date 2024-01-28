import dataclasses
import logging
import threading
import typing as t
import uuid

import ase
import ase.io
import numpy as np
import socketio
from znframe.frame import Frame

from zndraw.data import CeleryTaskData, FrameData, ModifierRegisterData
from zndraw.modify import UpdateScene, get_modify_class
from zndraw.settings import GlobalConfig
from zndraw.utils import (
    ZnDrawLoggingHandler,
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
    modifier_timeout : int
        Maximum runtime for modifiers in seconds.
    """

    call_timeout: int = 60
    retries: int = 100
    modifier_timeout: int = 60


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
    token: str = None
    display_new: bool = True
    _uuid: uuid.UUID = dataclasses.field(default_factory=uuid.uuid4)
    auth_token: str = None
    config: Config = dataclasses.field(default_factory=Config)
    available: bool = True

    _target_sid: str = None

    def __post_init__(self):
        self._uuid = str(self._uuid)
        self.socket = socketio.Client()
        if isinstance(self.config, dict):
            self.config = Config(**self.config)
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
        self.socket.on("available", self._on_available)

    def _on_available(self):
        return {"available": self.available, "timeout": self.config.modifier_timeout}

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
        log.critical("Reconnecting to server")
        self.socket.disconnect()
        self._connect()
        self.socket.sleep(0.1)

    def __len__(self) -> int:
        return int(self.socket.call("atoms:length", timeout=self.config.call_timeout))

    def __setitem__(self, index, value):
        if not isinstance(value, ase.Atoms) and not isinstance(value, Frame):
            raise ValueError("Must be an ase.Atoms or Frame object")

        assert isinstance(index, int), "Index must be an integer"
        if isinstance(value, ase.Atoms):
            value = Frame.from_atoms(value)
        self.socket.emit(
            "atoms:upload",
            dataclasses.asdict(
                FrameData(
                    index=index, data=value.to_dict(built_in_types=False), update=True
                )
            ),
        )

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

        downloaded_data = self.socket.call(
            "atoms:download", index, timeout=self.config.call_timeout
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
        data = self.socket.call("points:get", timeout=self.config.call_timeout)
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
        self.socket.emit("points:set", value)

    @property
    def segments(self) -> np.ndarray:
        data = self.socket.call("scene:segments", timeout=self.config.call_timeout)
        return np.array(data)

    @property
    def step(self) -> int:
        step = int(
            self.socket.call(
                "scene:step",
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
        return self.socket.call("selection:get", timeout=self.config.call_timeout)

    @selection.setter
    def selection(self, value: list[int]):
        self.socket.emit("selection:set", value)

    def play(self):
        self.socket.emit("scene:play")

    def pause(self):
        self.socket.emit("scene:pause")

    @property
    def figure(self):
        raise NotImplementedError("Gathering figure from webclient not implemented yet")

    @figure.setter
    def figure(self, fig: str):
        data = {"figure": fig, "token": self.token}
        self.socket.emit("analysis:figure", data)

    @property
    def bookmarks(self) -> dict:
        return self.socket.call("bookmarks:get", timeout=self.config.call_timeout)

    @bookmarks.setter
    def bookmarks(self, value: dict):
        self.socket.emit("bookmarks:set", value)

    def _pre_modifier_run(self, data) -> None:
        self.available = False
        log.critical(f"Modifier running {self.available = }")


        vis = type(self)(self.url, data["token"])

        msg = CeleryTaskData(
            target=f"webclients_{vis.token}", event="modifier:run:running", data=None
        )

        vis.socket.emit("celery:task:results", dataclasses.asdict(msg))
        try:
            config = GlobalConfig.load()
            cls = get_modify_class(
                config.get_modify_methods(
                    include=[x["cls"] for x in self._modifiers.values()]
                )
            )
            modifier = cls(**data["params"])
            modifier.run(
                vis,
                **self._modifiers[modifier.method.__class__.__name__]["run_kwargs"],
            )
        except Exception as err:
            vis.log(f"Modifier failed with error: {repr(err)}")
        msg = CeleryTaskData(
            target=f"{vis.token}", event="modifier:run:finished", data=None
        )
        vis.socket.emit("celery:task:results", dataclasses.asdict(msg))
        self.available = True
        log.critical(f"Modifier finished {self.available = }")

    def _modifier_run(self, data: dict) -> None:
        config = GlobalConfig.load()
        cls = get_modify_class(
            config.get_modify_methods(
                include=[x["cls"] for x in self._modifiers.values()]
            )
        )
        modifier = cls(**data)
        modifier.run(
            self,
            **self._modifiers[modifier.method.__class__.__name__]["run_kwargs"],
        )


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

        msg = ModifierRegisterData(
            schema=cls.model_json_schema(),
            name=cls.__name__,
            default=default,
        )

        self.socket.emit(
            "modifier:register",
            dataclasses.asdict(msg),
        )
        self._modifiers[cls.__name__] = {"cls": cls, "run_kwargs": run_kwargs}
