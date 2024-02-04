import dataclasses
import logging
import threading
import typing as t
from threading import Lock

import ase
import ase.io
import numpy as np
import socketio
from znframe.frame import Frame

from zndraw.data import CeleryTaskData, ModifierRegisterData
from zndraw.modify import UpdateScene, get_modify_class
from zndraw.settings import GlobalConfig
from zndraw.utils import (
    ZnDrawLoggingHandler,
)

from .base import ZnDrawBase
from .data import RoomGetData, RoomSetData
from .utils import (
    estimate_max_batch_size_for_socket,
    split_list_into_chunks,
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

    call_timeout: int = 3


class ZnDraw(ZnDrawBase):
    """

    Attributes
    ----------
    display_new : bool
        Display new atoms in the webclient, when they are added.

    token : str
        Identifies the session this instances is being connected to.
        Tokens can be shared.
    auth_token : str
        Authentication token, used e.g. for registering modifiers to all users and
        not just the current session.
    """

    display_new: bool = True
    auth_token: str | None = None
    config: Config = dataclasses.field(default_factory=Config)
    _modifiers: dict = dataclasses.field(default_factory=dict)

    _target_sid: str | None = None

    _lock = Lock()
    _data = None

    def __post_init__(self):
        if isinstance(self.config, dict):
            self.config = Config(**self.config)
        self.socket.on("disconnect", self._on_disconnect)
        self.socket.on("connect", self._on_connect)
        self.socket.on("modifier:run", self._pre_modifier_run)
        self.socket.on("message:log", lambda data: print(data))
        self.socket.on("room:get", lambda data: setattr(self, "_data", data))
        self.socket.on("room:set:finished", lambda *args: self._lock.release())
        super().__post_init__()

    def _on_disconnect(self):
        log.critical(f"Disconnected from server: {self._modifiers}")

    def _on_connect(self):
        log.critical(f"Connected to server: {self._modifiers} with token {self.token}")
        self.socket.emit(
            "join",
            {
                "token": self.token,
                "auth_token": self.auth_token,
            },
        )

    def _connect(self):
        for _ in range(100):
            try:
                self.socket.connect(self.url)
                break
            except socketio.exceptions.ConnectionError:
                self.socket.sleep(0.1)  # this can't work?
        else:
            raise socketio.exceptions.ConnectionError

    def reconnect(self) -> None:
        """Reconnect to the server."""
        log.critical("Reconnecting to server")
        self.socket.disconnect()
        self._connect()
        self._on_connect()
        self.socket.sleep(0.1)
        log.critical("Reconnected to server")

    def get_data(self, **data: dict) -> RoomGetData:
        data = RoomGetData(**data)
        with self._lock:
            self._data = None
            self.socket.emit("room:get", data.to_dict())
            while self._data is None:
                self.socket.sleep(seconds=1)
                # generous timeout
            # self._data.pop("update_database", None) # TODO: this should not happen
            data = RoomGetData(**self._data)
            self._data = None
            return data

    def set_data(self, **data: dict) -> None:
        data = RoomSetData(**data)
        self._lock.acquire(blocking=True)  # might need to have a while loop here
        self.socket.emit("room:set", data.to_dict())

    def __len__(self) -> int:
        return self.get_data(length=True).length

    def __setitem__(self, index, value):
        assert isinstance(index, int), "Index must be an integer"
        if isinstance(value, ase.Atoms):
            value = Frame.from_atoms(value)
        self.set_data(
            frames={index: value.to_dict(built_in_types=False)},
            step=index,
            update_database=True,
        )

    def __delitem__(self, index: int | slice | list[int]):
        if (
            isinstance(index, int)
            or isinstance(index, slice)
            or isinstance(index, list)
        ):
            length = len(self)
            index = self.wrap_and_check_index(index, length)
            self.set_data(frames={i: None for i in index}, update_database=True)
        else:
            raise TypeError("Index must be an integer, slice or list[int]")

    def insert(self, index, value: t.Union[ase.Atoms, Frame]):
        """Insert atoms before index"""
        if isinstance(value, ase.Atoms):
            value = Frame.from_atoms(value)
        log.warning(
            "Currently `insert` is very taxing on the server, use with caution!"
        )
        index = self.wrap_and_check_index(index, len(self))[0]
        data_after = self[index:]
        self[index] = value
        del self[index + 1 :]
        self.extend(data_after)

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
        if not isinstance(values, list):
            self[size] = values
        else:
            if isinstance(values[0], ase.Atoms):
                values = [Frame.from_atoms(val) for val in values]
            batch_size = estimate_max_batch_size_for_socket(values)
            indices = list(range(size, size + len(values)))
            all_data = [
                (i, val.to_dict(built_in_types=False))
                for i, val in zip(indices, values)
            ]

            for chunk in split_list_into_chunks(all_data, batch_size):
                batch = {tup[0]: tup[1] for tup in chunk}
                self.set_data(frames=batch, update_database=True)

    def __getitem__(self, index) -> t.Union[ase.Atoms, list[ase.Atoms]]:
        length = len(self)
        index = self.wrap_and_check_index(index, length)
        data = self.get_data(frames=index).frames

        atoms_list = []
        for idx, val in zip(index, data):
            if val is None:
                raise IndexError(f"Index {idx} out of range")
            atoms_list.append(Frame.from_dict(val).to_atoms())

        return_data = atoms_list[0] if len(index) == 1 else atoms_list
        return return_data

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
        return self[self.step]  # type: ignore

    @property
    def points(self) -> np.ndarray:
        data = self.get_data(points=True).points
        return np.array(data)

    @points.setter
    def points(self, value: t.Union[np.ndarray, list]) -> None:
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if len(value) > 0:
            try:
                assert len(value[0]) == 3
            except (TypeError, AssertionError):
                raise ValueError("Points must be a list of 3D coordinates")
        self.set_data(points=value, update_database=True)

    @property
    def segments(self) -> np.ndarray:
        data = self.get_data(segments=True).segments
        return np.array(data)

    @property
    def step(self) -> int:
        step = self.get_data(step=True).step
        return step

    @step.setter
    def step(self, index):
        index = self.wrap_and_check_index(index, len(self))[0]
        self.set_data(step=index, update_database=True)

    @property
    def selection(self) -> list[int]:
        return self.get_data(selection=True).selection

    @selection.setter
    def selection(self, value: list[int]):
        self.set_data(selection=value, update_database=True)

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
        return self.get_data(bookmarks=True).bookmarks

    @bookmarks.setter
    def bookmarks(self, value: dict):
        self.set_data(bookmarks=value, update_database=True)

    def _pre_modifier_run(self, data) -> None:
        self.socket.emit("modifier:available", False)
        vis = type(self)(self.url, data["token"])
        msg = CeleryTaskData(
            target=f"webclients_{vis.token}", event="modifier:run:running", data=None
        )

        vis.socket.emit("celery:task:emit", dataclasses.asdict(msg))
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
            target=f"{vis.token}",
            event="modifier:run:finished",
            data=None,
            disconnect=True,
        )
        vis.socket.emit("celery:task:emit", dataclasses.asdict(msg))
        self.socket.emit("modifier:available", True)
        print("Modifier finished!!!!!!!!")

    @staticmethod
    def wrap_and_check_index(index: int | slice | list[int], length: int) -> list[int]:
        is_slice = isinstance(index, slice)
        if is_slice:
            index = list(range(*index.indices(length)))
        index = [index] if isinstance(index, int) else index
        index = [i if i >= 0 else length + i for i in index]
        # check if index is out of range
        for i in index:
            if i >= length:
                raise IndexError(f"Index {i} out of range for length {length}")
        return index


@dataclasses.dataclass
class ZnDrawOld(ZnDrawBase):
    """ZnDraw client."""

    url: str = None

    jupyter: bool = False

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

    def reconnect(self) -> None:
        super().reconnect()

        for k, v in self._modifiers.items():
            log.critical(f"Re-registering modifier {k}")
            msg = ModifierRegisterData(
                schema=v["cls"].model_json_schema(),
                name=k,
                default=v["default"],
                timeout=v["run_kwargs"]["timeout"],
            )
            self.socket.emit(
                "modifier:register",
                dataclasses.asdict(msg),
            )

            self.socket.emit(
                "modifier:available",
                True,
            )

    def register_modifier(
        self,
        cls: UpdateScene,
        run_kwargs: dict = None,
        default: bool = False,
        timeout: float = 60,
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
        timeout : float, optional
            Timeout for the modifier to run in seconds. The Webclient
            will alert the user if the modifier takes longer than this time and
            release the modify button (no further changes are expected, but they
            can happen).
        """
        if run_kwargs is None:
            run_kwargs = {}
        run_kwargs["timeout"] = timeout
        if len(self._modifiers):
            raise ValueError(
                "Only one modifier can be registered at the moment. "
                "This is a limitation of the current implementation."
            )

        msg = ModifierRegisterData(
            schema=cls.model_json_schema(),
            name=cls.__name__,
            default=default,
            timeout=timeout,
        )

        self.socket.emit(
            "modifier:register",
            dataclasses.asdict(msg),
        )
        self._modifiers[cls.__name__] = {
            "cls": cls,
            "run_kwargs": run_kwargs,
            "default": default,
        }
        self.socket.emit(
            "modifier:available",
            True,
        )
