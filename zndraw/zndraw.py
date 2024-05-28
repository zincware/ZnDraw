import dataclasses
import json
import logging
import typing as t

import ase
import numpy as np
import socketio
import znframe

from zndraw.base import Extension, ZnDrawBase

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

    def __setitem__(self, index: int | list[int], value: ase.Atoms | list[ase.Atoms]):
        if isinstance(index, int):
            data = {index: znframe.Frame.from_atoms(value).to_json()}
        else:
            data = {
                i: znframe.Frame.from_atoms(val).to_json()
                for i, val in zip(index, value)
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
        return {int(k): v for k, v in self.socket.call("room:bookmarks:get").items()}

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
