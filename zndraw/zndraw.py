import dataclasses
import json
import logging
import typing as t

import ase
import numpy as np
import socketio
import tqdm
import znframe

from zndraw.base import Extension, ZnDrawBase
from zndraw.draw import Geometry, Object3D

log = logging.getLogger(__name__)


class RegisterModifier(t.TypedDict):
    cls: t.Type[Extension]
    run_kwargs: dict
    public: bool
    frozen: bool
    timeout: float


class TimeoutConfig(t.TypedDict):
    """Timeout configuration for the ZnDraw client."""

    connection: int
    modifier: float


def _register_modifier(vis: "ZnDraw", data: RegisterModifier) -> None:
    log.debug(f"Registering modifier `{data['cls'].__name__}`")
    vis.socket.emit(
        "modifier:register",
        {
            "schema": data["cls"].model_json_schema(),
            "name": data["cls"].__name__,
            "public": data["public"],
            "timeout": data["timeout"],
        },
    )
    if vis._available:
        vis.socket.emit("modifier:available", list(vis._modifiers))


@dataclasses.dataclass
class ZnDraw(ZnDrawBase):
    url: str
    token: str | None = None
    auth_token: str | None = None

    socket: socketio.Client = dataclasses.field(
        default_factory=socketio.Client, repr=False
    )
    timeout: TimeoutConfig = dataclasses.field(
        default_factory=lambda: {"connection": 10, "modifier": 0.25}
    )
    maximum_message_size: int = dataclasses.field(default=1_000_000, repr=False)

    _modifiers: dict[str, RegisterModifier] = dataclasses.field(default_factory=dict)
    _available: bool = True

    def __post_init__(self):
        def on_wakeup():
            if self._available:
                self.socket.emit("modifier:available", list(self._modifiers))

        self.url = self.url.replace("http", "ws")
        self.socket.on("connect", self._on_connect)
        self.socket.on("modifier:run", self._run_modifier)
        self.socket.on("modifier:wakeup", on_wakeup)
        self.socket.on("room:log", lambda x: print(x))

        self.socket.connect(self.url, wait_timeout=self.timeout["connection"])

    def _on_connect(self):
        log.debug("Connected to ZnDraw server")
        self.socket.emit(
            "join",
            {
                "token": str(self.token),
                "auth_token": self.auth_token,
            },
        )
        for data in self._modifiers.values():
            _register_modifier(self, data)

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

    def __setitem__(self, index: int | list[int], value: ase.Atoms | list[ase.Atoms]):
        if isinstance(index, int):
            data = {index: znframe.Frame.from_atoms(value).to_json()}
        else:
            data = {
                i: znframe.Frame.from_atoms(val).to_json() for i, val in zip(index, value)
            }

        self.socket.emit("room:frames:set", data)

    def __len__(self) -> int:
        return int(self.socket.call("room:length:get"))

    def __delitem__(self, index: int | slice | list[int]):
        if isinstance(index, int):
            index = [index]
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        self.socket.emit("room:frames:delete", index)

    def insert(self, index: int, value: ase.Atoms):
        self.socket.emit(
            "room:frames:insert",
            {"index": index, "value": znframe.Frame.from_atoms(value).to_json()},
        )

    def extend(self, values: list[ase.Atoms]):
        msg = {}

        # enable tbar if more than 10 messages are sent
        # approximated by the size of the first frame
        show_tbar = (
            len(values)
            * len(
                json.dumps(znframe.Frame.from_atoms(values[0]).to_json()).encode("utf-8")
            )
        ) > (10 * self.maximum_message_size)
        tbar = tqdm.tqdm(
            values, desc="Sending frames", unit=" frame", disable=not show_tbar
        )

        for i, val in enumerate(tbar, start=len(self)):
            msg[i] = znframe.Frame.from_atoms(val).to_json()
            if len(json.dumps(msg).encode("utf-8")) > self.maximum_message_size:
                self.socket.emit("room:frames:set", msg)
                msg = {}
                # after each large message, wait a bit
                self.socket.sleep(self.timeout["modifier"])
        if msg:  # Only send the message if it's not empty
            self.socket.emit("room:frames:set", msg)

    @property
    def selection(self) -> list[int]:
        return self.socket.call("room:selection:get")["0"]
        data: dict = self.socket.call("room:selection:get")
        return data["0"]

    @selection.setter
    def selection(self, value: list[int]):
        self.socket.emit("room:selection:set", {"0": value})

    @property
    def step(self) -> int:
        return int(self.socket.call("room:step:get"))

    def log(self, message: str) -> None:
        self.socket.emit("room:log", message)

    @step.setter
    def step(self, value: int):
        if value < 0:
            raise ValueError("Step must be positive")
        if value >= len(self):
            raise ValueError("Step out of range")
        # Have only one step per room!
        # shared rooms are rare anyhow and making per-client steps
        # and room hosts is annoying
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
        return np.array(self.socket.call("room:points:get")["0"])

    @points.setter
    def points(self, points: np.ndarray) -> None:
        self.socket.emit("room:points:set", {"0": np.array(points).tolist()})

    @property
    def bookmarks(self) -> dict:
        return {int(k): v for k, v in self.socket.call("room:bookmarks:get").items()}

    @bookmarks.setter
    def bookmarks(self, value: dict):
        self.socket.emit("room:bookmarks:set", value)

    @property
    def camera(self) -> dict:
        return self.socket.call("room:camera:get")

    @camera.setter
    def camera(self, value: dict):
        raise NotImplementedError("Changing camera position is currently not possible")
        self.socket.emit("room:camera:set", value)

    @property
    def geometries(self) -> list[Object3D]:
        # return self.socket.call("room:geometry:get")

        return [Geometry(method=x).method for x in self.socket.call("room:geometry:get")]

    @geometries.setter
    def geometries(self, value: list[Object3D]):
        self.socket.emit("room:geometry:set", [x.model_dump() for x in value])

    def register_modifier(
        self,
        cls: t.Type[Extension],
        run_kwargs: dict | None = None,
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
        if cls.__name__ in self._modifiers:
            raise ValueError(f"Modifier {cls.__name__} already registered")

        self._modifiers[cls.__name__] = {
            "cls": cls,
            "run_kwargs": run_kwargs,
            "public": public,
            "frozen": use_frozen,
            "timeout": timeout,
        }

        _register_modifier(self, self._modifiers[cls.__name__])

    def _run_modifier(self, data: dict):
        self._available = False
        room = data.pop("ZNDRAW_CLIENT_ROOM")
        vis = type(self)(
            url=self.url, token=room, maximum_message_size=self.maximum_message_size
        )

        try:
            # TODO: for public modifiers the vis object must not be in the same room, create a new one!!!!!
            vis.socket.emit("room:modifier:queue", 0)
            name = data["method"]["discriminator"]

            instance = self._modifiers[name]["cls"](**data["method"])
            instance.run(
                vis,
                timeout=self._modifiers[name]["timeout"],
                **self._modifiers[name]["run_kwargs"],
            )
        except Exception as e:
            log.exception(e)
            vis.log(f"Error: {e}")
        finally:
            vis.socket.emit("room:modifier:queue", -1)

            # wait and then disconnect
            vis.socket.sleep(1)
            vis.socket.disconnect()

            self._available = True
            self.socket.emit("modifier:available", list(self._modifiers))
