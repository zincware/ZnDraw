import dataclasses
import json
import logging
import typing as t

import ase
import numpy as np
import socketio
import znframe

from zndraw.modify import UpdateScene
from zndraw.base import ZnDrawBase

log = logging.getLogger(__name__)


@dataclasses.dataclass
class ZnDraw(ZnDrawBase):
    url: str
    token: str | None = None
    auth_token: str | None = None

    socket: socketio.Client = dataclasses.field(
        default_factory=socketio.Client, repr=False
    )

    _modifiers: dict = dataclasses.field(default_factory=dict)
    _available: bool = True

    def __post_init__(self):
        def on_wakeup():
            if self._available:
                self.socket.emit("modifier:available", list(self._modifiers))

        self.url = self.url.replace("http", "ws")
        self.socket.on("connect", self._on_connect)
        self.socket.connect(self.url, wait_timeout=1)
        self.socket.on("modifier:run", self._run_modifier)
        self.socket.on("modifier:wakeup", on_wakeup)

    def _on_connect(self):
        self.socket.emit(
            "join",
            {
                "token": str(self.token),
                "auth_token": self.auth_token,
            },
        )

    def __getitem__(self, index) -> ase.Atoms | list[ase.Atoms]:
        single_item = isinstance(index, int)
        if single_item:
            index = [index]
        
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        
        if any(x < 0 for x in index):
            raise IndexError("Index must be positive")
        if any(x >= len(self) for x in index):
            raise IndexError("Index out of range")

        data: dict = self.socket.call("room:frames:get", index)
        structures = [znframe.Frame(**x).to_atoms() for x in data.values()]
        return structures[0] if single_item else structures

    def __setitem__(self, index: list[int], value: list[ase.Atoms]):
        # TODO: send in chunks
        data = {i: znframe.Frame.from_atoms(x).to_json() for i, x in zip(index, value)}

        self.socket.emit("room:frames:set", data)

    def __len__(self) -> int:
        return int(self.socket.call("room:length:get"))
    
    def __delitem__(self, index: int | slice | list[int]):
        if isinstance(index, int):
            index = [index]
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        self.socket.emit("room:frames:set", {i: None for i in index})
    
    def insert(self, index: int, value: ase.Atoms):
        self.socket.emit("room:frames:insert", {index: znframe.Frame.from_atoms(value).to_json()})

    @property
    def selection(self) -> list[int]:
        data: dict = self.socket.call("room:selection:get")
        return json.loads(data["0"])

    @selection.setter
    def selection(self, value: list[int]):
        self.socket.emit("room:selection:set", {"0": value})

    @property
    def step(self) -> int:
        return int(self.socket.call("room:step:get"))

    def log(self, message: str) -> None:
        self.socket.emit("message:log", message)

    @step.setter
    def step(self, value: int):
        # Have only one step per room!
        # shared rooms are rare anyhow and making per-client steps
        # and room hosts is anoying
        # what about the camera?
        # or collect the steps of all clients in a dict
        # and save the host and go from there, also fine and not too much worker.
        self.socket.emit("room:step:set", value)

    @property
    def figure(self) -> str:
        return self.socket.call("analysis:figure:get")

    @figure.setter
    def figure(self, fig: str) -> None:
        self.socket.emit("analysis:figure:set", fig)

    @property
    def atoms(self) -> ase.Atoms:
        return self[self.step]

    @atoms.setter
    def atoms(self, atoms: ase.Atoms) -> None:
        self[[self.step]] = [atoms]

    @property
    def points(self) -> list[list[float]]:
        return np.array(json.loads(self.socket.call("room:points:get")["0"]))

    @points.setter
    def points(self, points: np.ndarray) -> None:
        self.socket.emit("room:points:set", {"0": np.array(points).tolist()})

    @property
    def bookmarks(self) -> dict:
        return json.loads(self.socket.call("room:bookmarks:get"))
    
    @bookmarks.setter
    def bookmarks(self, value: dict):
        self.socket.emit("room:bookmarks:set", value)

    @property
    def camera(self) -> dict:
        return json.loads(self.socket.call("room:camera:get"))
    
    @camera.setter
    def camera(self, value: dict):
        self.socket.emit("room:camera:set", value)
    

    def register_modifier(
        self,
        cls: t.Type[UpdateScene],
        run_kwargs: dict = None,
        public: bool = False,
        timeout: float = 60,
        use_frozen: bool = True,
    ):
        """Register a modifier class.

        Attributes
        ----------
        cls : UpdateScene
            The modifier class to register.
        run_kwargs : dict, optional
            Keyword arguments to pass to the run method of the modifier class.
        public : bool, optional
            Whether to enable the modifier for ALL sessions of the ZnDraw client,
            or just the session for the given token.
        timeout : float, optional
            Timeout for the modifier to run in seconds. The Webclient
            will alert the user if the modifier takes longer than this time and
            release the modify button (no further changes are expected, but they
            can happen).
        use_frozen : bool, default=True
            Whether to use the ZnDrawFrozen class to run the modifier.
            The frozen class only allows provides cached data and
            e.g. access to other steps than the current one is not possible.
            If set to false, a full ZnDraw instance will be created for the modifier.
            This can have a performance impact and may lead to timeouts.

        """
        if timeout < 1:
            raise ValueError("Timeout must be at least 1 second")
        if timeout > 300:
            log.critical(
                "Timeout is set to more than 300 seconds. Modifiers might be killed automatically."
            )
        if run_kwargs is None:
            run_kwargs = {}
        run_kwargs["timeout"] = timeout
        if cls.__name__ in self._modifiers:
            raise ValueError(f"Modifier {cls.__name__} already registered")

        self.socket.emit(
            "modifier:register",
            {
                "schema": cls.model_json_schema(),
                "name": cls.__name__,
                "public": public,
                "timeout": timeout,
            },
        )
        self._modifiers[cls.__name__] = {
            "cls": cls,
            "run_kwargs": run_kwargs,
            "public": public,
            "frozen": use_frozen,
        }
        if self._available:
            self.socket.emit("modifier:available", list(self._modifiers))

    def _run_modifier(self, data: dict):
        self._available = False
        room = data.pop("ZNDRAW_CLIENT_ROOM")
        vis = type(self)(url=self.url, token=room)

        try:
            # TODO: for public modifiers the vis object must not be in the same room, create a new one!!!!!
            vis.socket.emit("modifier:run:running")
            name = data["method"]["discriminator"]

            instance = self._modifiers[name]["cls"](**data["method"])
            instance.run(vis, **self._modifiers[name]["run_kwargs"])

            vis.socket.emit("modifier:run:finished")
        finally:
            self._available = True
            self.socket.emit("modifier:available", list(self._modifiers))


# import logging
# import typing as t
# from threading import Lock

# import ase
# import ase.io
# import numpy as np
# import socketio
# from socketio import exceptions as socketio_exceptions
# from znframe.frame import Frame

# from zndraw.data import CeleryTaskData, ModifierRegisterData
# from zndraw.modify import UpdateScene, get_modify_class
# from zndraw.settings import GlobalConfig

# from .base import ZnDrawBase
# from .data import RoomGetData, RoomSetData
# from .utils import (
#     check_selection,
#     estimate_max_batch_size_for_socket,
#     split_list_into_chunks,
#     wrap_and_check_index,
# )
# from .zndraw_frozen import ZnDrawFrozen

# log = logging.getLogger(__name__)


# @dataclasses.dataclass
# class Config:
#     """Configuration for ZnDraw client.

#     Attributes
#     ----------
#     timeout : int
#         Timeout for socket calls in seconds.
#         Set to a smaller value to fail faster.
#     """

#     timeout: int = 3


# @dataclasses.dataclass
# class ZnDraw(ZnDrawBase):
#     """

#     Attributes
#     ----------
#     display_new : bool
#         Display new atoms in the webclient, when they are added.

#     token : str
#         Identifies the session this instances is being connected to.
#         Tokens can be shared.
#     auth_token : str
#         Authentication token, used e.g. for registering modifiers to all users and
#         not just the current session.
#     """

#     token: str = None
#     display_new: bool = True
#     auth_token: str | None = None
#     config: Config = dataclasses.field(default_factory=Config)
#     _modifiers: dict = dataclasses.field(default_factory=dict)

#     _target_sid: str | None = None

#     _lock: Lock = dataclasses.field(default_factory=Lock)
#     _data = None
#     _available = True

#     def __post_init__(self):
#         if isinstance(self.config, dict):
#             self.config = Config(**self.config)
#         self.socket.on("disconnect", self._on_disconnect)
#         self.socket.on("modifier:run", self._pre_modifier_run)
#         self.socket.on("message:log", lambda data: print(data))
#         self.socket.on("room:get", lambda data: setattr(self, "_data", data))
#         self.socket.on("room:set:finished", lambda *args: self._lock.release())
#         super().__post_init__()

#     def _on_disconnect(self):
#         log.info("Disconnected from server")

#     def _on_connect(self):
#         log.info("Connected to server")
#         self.socket.emit(
#             "join",
#             {
#                 "token": self.token,
#                 "auth_token": self.auth_token,
#             },
#         )
#         for k, v in self._modifiers.items():
#             log.debug(f"Re-registering modifier {k}")
#             msg = ModifierRegisterData(
#                 schema=v["cls"].model_json_schema(),
#                 name=k,
#                 default=v["default"],
#                 timeout=v["run_kwargs"]["timeout"],
#             )
#             self.socket.emit(
#                 "modifier:register",
#                 dataclasses.asdict(msg),
#             )
#             self.socket.emit(
#                 "modifier:available",
#                 self._available,
#             )

#     def _connect(self):
#         for _ in range(100):
#             try:
#                 self.socket.connect(self.url)
#                 break
#             except socketio.exceptions.ConnectionError:
#                 self.socket.sleep(0.1)  # this can't work?
#         else:
#             raise socketio.exceptions.ConnectionError

#     def get_data(self, **data: dict) -> RoomGetData:
#         data = RoomGetData(**data)
#         with self._lock:
#             self._data = None
#             self.socket.emit("room:get", data.to_dict())
#             for _ in range(self.config.timeout):
#                 if self._data is not None:
#                     break
#                 self.socket.sleep(seconds=1)
#             else:
#                 raise TimeoutError("Timeout while waiting for data")
#             # self._data.pop("update_database", None) # TODO: this should not happen
#             data = RoomGetData(**self._data)
#             self._data = None
#             return data

#     def set_data(self, **data: dict) -> None:
#         data = RoomSetData(**data)
#         self._lock.acquire(blocking=True)  # might need to have a while loop here
#         self.socket.emit("room:set", data.to_dict())

#     def __len__(self) -> int:
#         return self.get_data(length=True).length

#     def __setitem__(self, index, value):
#         assert isinstance(index, int), "Index must be an integer"
#         if isinstance(value, ase.Atoms):
#             value = Frame.from_atoms(value)
#         self.set_data(
#             frames={index: value.to_dict(built_in_types=False)},
#             step=index,
#             update_database=True,
#         )

#     def __delitem__(self, index: int | slice | list[int]):
#         if (
#             isinstance(index, int)
#             or isinstance(index, slice)
#             or isinstance(index, list)
#         ):
#             length = len(self)
#             index = wrap_and_check_index(index, length)
#             self.set_data(frames={i: None for i in index}, update_database=True)
#         else:
#             raise TypeError("Index must be an integer, slice or list[int]")

#     def insert(self, index, value: t.Union[ase.Atoms, Frame]):
#         """Insert atoms before index"""
#         if isinstance(value, ase.Atoms):
#             value = Frame.from_atoms(value)
#         log.warning(
#             "Currently `insert` is very taxing on the server, use with caution!"
#         )
#         index = wrap_and_check_index(index, len(self))[0]
#         data_after = self[index:]
#         self[index] = value
#         del self[index + 1 :]
#         self.extend(data_after)

#     def append(self, value: t.Union[ase.Atoms, Frame]) -> None:
#         """Append atoms to the end of the list"""
#         if isinstance(value, ase.Atoms):
#             value = Frame.from_atoms(value)
#         self[len(self)] = value

#     def extend(
#         self, values: t.Union[ase.Atoms, Frame, list[ase.Atoms], list[Frame]]
#     ) -> None:
#         """Extend the list by appending all the items in the given list"""
#         size = len(self)
#         if not isinstance(values, list):
#             self[size] = values
#         else:
#             if isinstance(values[0], ase.Atoms):
#                 values = [Frame.from_atoms(val) for val in values]
#             batch_size = estimate_max_batch_size_for_socket(values)
#             indices = list(range(size, size + len(values)))
#             all_data = [
#                 (i, val.to_dict(built_in_types=False))
#                 for i, val in zip(indices, values)
#             ]

#             for chunk in split_list_into_chunks(all_data, batch_size):
#                 batch = {tup[0]: tup[1] for tup in chunk}
#                 self.set_data(frames=batch, update_database=True)
#             self.set_data(step=max(indices), update_database=True)

#     def __getitem__(self, index) -> t.Union[ase.Atoms, list[ase.Atoms]]:
#         length = len(self)
#         index = wrap_and_check_index(index, length)
#         data = self.get_data(frames=index).frames

#         atoms_list = []
#         for idx, val in zip(index, data):
#             if val is None:
#                 raise IndexError(f"Index {idx} out of range")
#             atoms_list.append(Frame.from_dict(val).to_atoms())

#         return_data = atoms_list[0] if len(index) == 1 else atoms_list
#         return return_data

#     def log(self, message: str) -> None:
#         """Log a message to the console"""
#         print(message)
#         self.socket.emit(
#             "message:log",
#             {
#                 "message": message,
#                 "token": self.token,
#             },
#         )

#     @property
#     def atoms(self) -> ase.Atoms:
#         """Return the atoms at the current step."""
#         return self[self.step]  # type: ignore

#     @property
#     def points(self) -> np.ndarray:
#         data = self.get_data(points=True).points
#         return np.array(data)

#     @points.setter
#     def points(self, value: t.Union[np.ndarray, list]) -> None:
#         if isinstance(value, np.ndarray):
#             value = value.tolist()
#         if len(value) > 0:
#             try:
#                 assert len(value[0]) == 3
#             except (TypeError, AssertionError):
#                 raise ValueError("Points must be a list of 3D coordinates")
#         self.set_data(points=value, update_database=True)

#     @property
#     def segments(self) -> np.ndarray:
#         return self.calculate_segments(self.points)

#     @property
#     def step(self) -> int:
#         step = self.get_data(step=True).step
#         return step

#     @step.setter
#     def step(self, index):
#         if not isinstance(index, int):
#             raise TypeError("Index must be an integer")
#         if index < 0:
#             raise IndexError(f"Index {index} out of range")
#         if index >= len(self):
#             raise IndexError(f"Index {index} out of range")
#         self.set_data(step=index, update_database=True)

#     @property
#     def selection(self) -> list[int]:
#         return self.get_data(selection=True).selection

#     @selection.setter
#     def selection(self, value: list[int]):
#         check_selection(value, len(self.atoms))
#         self.set_data(selection=value, update_database=True)

#     def play(self):
#         self.socket.emit("scene:play")

#     def pause(self):
#         self.socket.emit("scene:pause")

#     @property
#     def figure(self):
#         raise NotImplementedError("Gathering figure from webclient not implemented yet")

#     @figure.setter
#     def figure(self, fig: str):
#         data = {"figure": fig, "token": self.token}
#         self.socket.emit("analysis:figure", data)

#     @property
#     def bookmarks(self) -> dict:
#         return self.get_data(bookmarks=True).bookmarks

#     @property
#     def camera(self) -> dict:
#         return self.get_data(camera=True).camera

#     @camera.setter
#     def camera(self, camera: dict):
#         """Set the camera position and orientation

#         camera: dict
#             A dictionary with the following
#             - position: list[float]
#                 The position of the camera
#             - target: list[float]
#                 The target of the camera
#         """
#         if set(camera) != {"position", "target"}:
#             raise ValueError("camera must have keys 'position' and 'target'")
#         msg = CeleryTaskData(
#             target=str(self.token),
#             event="camera:update",
#             data=camera,
#         )
#         self.socket.emit("celery:task:emit", msg.to_dict())

#     @bookmarks.setter
#     def bookmarks(self, value: dict):
#         self.set_data(bookmarks=value, update_database=True)

#     def _pre_modifier_run(self, data) -> None:
#         self._available = False
#         try:
#             self._modifier_run(data)
#         except Exception as err:
#             self.log(f"Modifier failed with error: {repr(err)}")
#         finally:
#             self._available = (
#                 True  # always execute this, even if an exception is raised
#             )
#         self.socket.emit("modifier:available", self._available)

#     def _modifier_run(self, data: dict) -> None:
#         self.socket.emit("modifier:available", self._available)
#         msg = CeleryTaskData(
#             target=f"webclients_{data['token']}",
#             event="modifier:run:running",
#             data=None,
#         )

#         self.socket.emit("celery:task:emit", dataclasses.asdict(msg))
#         try:
#             config = GlobalConfig.load()
#             cls = get_modify_class(
#                 config.get_modify_methods(
#                     include=[x["cls"] for x in self._modifiers.values()]
#                 )
#             )
#             modifier = cls(**data["params"])
#             use_frozen = self._modifiers[modifier.method.__class__.__name__]["frozen"]
#             if use_frozen:
#                 vis = ZnDrawFrozen(
#                     url=self.url, token=data["token"], cached_data=data["cache"]
#                 )
#             else:
#                 vis = ZnDraw(url=self.url, token=data["token"])
#             try:
#                 modifier.run(
#                     vis,
#                     **self._modifiers[modifier.method.__class__.__name__]["run_kwargs"],
#                 )
#             except Exception as err:
#                 vis.log(f"Modifier failed with error: {repr(err)}")

#             vis.socket.sleep(1)
#             vis.socket.disconnect()
#         except (
#             socketio_exceptions.ConnectionError,
#             socketio_exceptions.BadNamespaceError,
#         ) as err:
#             msg = CeleryTaskData(
#                 target=f"webclients_{data['token']}",
#                 event="message:log",
#                 data="Could not establish connection " + str(err),
#                 disconnect=False,
#             )
#             self.socket.emit("celery:task:emit", dataclasses.asdict(msg))
#         msg = CeleryTaskData(
#             target=f"{data['token']}",
#             event="modifier:run:finished",
#             data=None,
#             disconnect=False,
#         )
#         self.socket.emit("celery:task:emit", dataclasses.asdict(msg))
