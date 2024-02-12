import dataclasses
import logging
import typing as t
from threading import Lock

import ase
import ase.io
import numpy as np
import socketio
from socketio import exceptions as socketio_exceptions
from znframe.frame import Frame

from zndraw.data import CeleryTaskData, ModifierRegisterData
from zndraw.modify import UpdateScene, get_modify_class
from zndraw.settings import GlobalConfig

from .base import ZnDrawBase
from .data import RoomGetData, RoomSetData
from .utils import (
    check_selection,
    estimate_max_batch_size_for_socket,
    split_list_into_chunks,
    wrap_and_check_index,
)
from .zndraw_frozen import ZnDrawFrozen

log = logging.getLogger(__name__)


@dataclasses.dataclass
class Config:
    """Configuration for ZnDraw client.

    Attributes
    ----------
    timeout : int
        Timeout for socket calls in seconds.
        Set to a smaller value to fail faster.
    """

    timeout: int = 3


@dataclasses.dataclass
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

    token: str = None
    display_new: bool = True
    auth_token: str | None = None
    config: Config = dataclasses.field(default_factory=Config)
    _modifiers: dict = dataclasses.field(default_factory=dict)

    _target_sid: str | None = None

    _lock: Lock = dataclasses.field(default_factory=Lock)
    _data = None
    _available = True

    def __post_init__(self):
        if isinstance(self.config, dict):
            self.config = Config(**self.config)
        self.socket.on("disconnect", self._on_disconnect)
        self.socket.on("modifier:run", self._pre_modifier_run)
        self.socket.on("message:log", lambda data: print(data))
        self.socket.on("room:get", lambda data: setattr(self, "_data", data))
        self.socket.on("room:set:finished", lambda *args: self._lock.release())
        super().__post_init__()

    def _on_disconnect(self):
        log.info(f"Disconnected from server")

    def _on_connect(self):
        log.info(f"Connected to server")
        self.socket.emit(
            "join",
            {
                "token": self.token,
                "auth_token": self.auth_token,
            },
        )
        for k, v in self._modifiers.items():
            log.debug(f"Re-registering modifier {k}")
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
                self._available,
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

    def get_data(self, **data: dict) -> RoomGetData:
        data = RoomGetData(**data)
        with self._lock:
            self._data = None
            self.socket.emit("room:get", data.to_dict())
            for _ in range(self.config.timeout):
                if self._data is not None:
                    break
                self.socket.sleep(seconds=1)
            else:
                raise TimeoutError("Timeout while waiting for data")
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
            index = wrap_and_check_index(index, length)
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
        index = wrap_and_check_index(index, len(self))[0]
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
        index = wrap_and_check_index(index, length)
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
        return self.calculate_segments(self.points)

    @property
    def step(self) -> int:
        step = self.get_data(step=True).step
        return step

    @step.setter
    def step(self, index):
        if not isinstance(index, int):
            raise TypeError("Index must be an integer")
        if index < 0:
            raise IndexError(f"Index {index} out of range")
        if index >= len(self):
            raise IndexError(f"Index {index} out of range")
        self.set_data(step=index, update_database=True)

    @property
    def selection(self) -> list[int]:
        return self.get_data(selection=True).selection

    @selection.setter
    def selection(self, value: list[int]):
        check_selection(value, len(self.atoms))
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

    @property
    def camera(self) -> dict:
        return self.get_data(camera=True).camera

    @camera.setter
    def camera(self, camera: dict):
        """Set the camera position and orientation

        camera: dict
            A dictionary with the following
            - position: list[float]
                The position of the camera
            - target: list[float]
                The target of the camera
        """
        if set(camera) != {"position", "target"}:
            raise ValueError("camera must have keys 'position' and 'target'")
        msg = CeleryTaskData(
            target=str(self.token),
            event="camera:update",
            data=camera,
        )
        self.socket.emit("celery:task:emit", msg.to_dict())

    @bookmarks.setter
    def bookmarks(self, value: dict):
        self.set_data(bookmarks=value, update_database=True)

    def _pre_modifier_run(self, data) -> None:
        self._available = False
        self.socket.emit("modifier:available", self._available)
        msg = CeleryTaskData(
            target=f"webclients_{data['token']}",
            event="modifier:run:running",
            data=None,
        )

        self.socket.emit("celery:task:emit", dataclasses.asdict(msg))
        try:
            vis = ZnDrawFrozen(
                url=self.url, token=data["token"], cached_data=data["cache"]
            )
            config = GlobalConfig.load()
            cls = get_modify_class(
                config.get_modify_methods(
                    include=[x["cls"] for x in self._modifiers.values()]
                )
            )
            modifier = cls(**data["params"])
            try:
                modifier.run(
                    vis,
                    **self._modifiers[modifier.method.__class__.__name__]["run_kwargs"],
                )
            except Exception as err:
                vis.log(f"Modifier failed with error: {repr(err)}")

            vis.socket.sleep(1)
            vis.socket.disconnect()
        except (
            socketio_exceptions.ConnectionError,
            socketio_exceptions.BadNamespaceError,
        ) as err:
            msg = CeleryTaskData(
                target=f"webclients_{data['token']}",
                event="message:log",
                data="Could not establish connection " + str(err),
                disconnect=False,
            )
            self.socket.emit("celery:task:emit", dataclasses.asdict(msg))

        msg = CeleryTaskData(
            target=f"{data['token']}",
            event="modifier:run:finished",
            data=None,
            disconnect=False,
        )
        self.socket.emit("celery:task:emit", dataclasses.asdict(msg))
        self._available = True
        self.socket.emit("modifier:available", self._available)
        print("Modifier finished!!!!!!!!")

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
        if cls.__name__ in self._modifiers:
            raise ValueError(f"Modifier {cls.__name__} already registered")

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
            self._available,
        )
