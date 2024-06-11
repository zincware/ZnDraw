import dataclasses
import datetime
import json
import logging
import typing as t

import ase
import numpy as np
import socketio.exceptions
import tqdm
import znframe

from zndraw.base import Extension, ZnDrawBase
from zndraw.draw import Geometry, Object3D
from zndraw.utils import call_with_retry, emit_with_retry

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
    between_calls: float

    emit_retries: int
    call_retries: int
    connect_retries: int


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
        default_factory=lambda: TimeoutConfig(
            connection=10,
            modifier=0.25,
            between_calls=0.1,
            emit_retries=3,
            call_retries=3,
            connect_retries=3,
        )
    )
    maximum_message_size: int = dataclasses.field(default=1_000_000, repr=False)

    _modifiers: dict[str, RegisterModifier] = dataclasses.field(default_factory=dict)
    _available: bool = True
    _last_call: datetime.datetime = dataclasses.field(
        default_factory=datetime.datetime.now
    )

    def __post_init__(self):
        def on_wakeup():
            if self._available:
                self.socket.emit("modifier:available", list(self._modifiers))

        self.url = self.url.replace("http", "ws")
        self.socket.on("connect", self._on_connect)
        self.socket.on("modifier:run", self._run_modifier)
        self.socket.on("modifier:wakeup", on_wakeup)
        self.socket.on("room:log", lambda x: print(x))

        for idx in range(self.timeout["connect_retries"] + 1):
            try:
                self.socket.connect(self.url, wait_timeout=self.timeout["connection"])
                break
            except socketio.exceptions.ConnectionError as err:
                log.warning("Connection failed. Retrying...")
                self._delay_socket()
                if idx == self.timeout["connect_retries"]:
                    raise err

    def _on_connect(self):
        log.debug("Connected to ZnDraw server")

        emit_with_retry(
            self.socket,
            "join",
            {
                "token": str(self.token),
                "auth_token": self.auth_token,
            },
            retries=self.timeout["emit_retries"],
        )

        for data in self._modifiers.values():
            _register_modifier(self, data)

    def _delay_socket(self):
        """Delay if the last call was too recent."""
        while (datetime.datetime.now() - self._last_call).total_seconds() < self.timeout[
            "between_calls"
        ]:
            self.socket.sleep(self.timeout["between_calls"] / 10)
        self._last_call = datetime.datetime.now()

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

        self._delay_socket()

        data = call_with_retry(
            self.socket,
            "room:frames:get",
            index,
            retries=self.timeout["call_retries"],
        )

        structures = [znframe.Frame(**x).to_atoms() for x in data.values()]
        return structures[0] if single_item else structures

    def __setitem__(self, index: int | list[int], value: ase.Atoms | list[ase.Atoms]):
        if isinstance(index, int):
            data = {index: znframe.Frame.from_atoms(value).to_json()}
        else:
            data = {
                i: znframe.Frame.from_atoms(val).to_json() for i, val in zip(index, value)
            }

        emit_with_retry(
            self.socket,
            "room:frames:set",
            data,
            retries=self.timeout["emit_retries"],
        )

    def __len__(self) -> int:
        self._delay_socket()

        return call_with_retry(
            self.socket,
            "room:length:get",
            retries=self.timeout["call_retries"],
        )

    def __delitem__(self, index: int | slice | list[int]):
        if isinstance(index, int):
            index = [index]
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))

        emit_with_retry(
            self.socket,
            "room:frames:delete",
            index,
            retries=self.timeout["emit_retries"],
        )

    def insert(self, index: int, value: ase.Atoms):
        emit_with_retry(
            self.socket,
            "room:frames:insert",
            {"index": index, "value": znframe.Frame.from_atoms(value).to_json()},
            retries=self.timeout["emit_retries"],
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
                emit_with_retry(
                    self.socket,
                    "room:frames:set",
                    msg,
                    retries=self.timeout["emit_retries"],
                )
                msg = {}
                # after each large message, wait a bit
                self.socket.sleep(self.timeout["modifier"])
        if msg:  # Only send the message if it's not empty
            for idx in range(self.timeout["emit_retries"] + 1):
                emit_with_retry(
                    self.socket,
                    "room:frames:set",
                    msg,
                    retries=self.timeout["emit_retries"],
                )

    @property
    def selection(self) -> list[int]:
        self._delay_socket()
        return call_with_retry(
            self.socket,
            "room:selection:get",
            retries=self.timeout["call_retries"],
        )["0"]

    @selection.setter
    def selection(self, value: list[int]):
        emit_with_retry(
            self.socket,
            "room:selection:set",
            {"0": value},
            retries=self.timeout["emit_retries"],
        )

    @property
    def step(self) -> int:
        self._delay_socket()

        return int(
            call_with_retry(
                self.socket,
                "room:step:get",
                retries=self.timeout["call_retries"],
            )
        )

    def log(self, message: str) -> None:
        emit_with_retry(
            self.socket,
            "room:log",
            message,
            retries=self.timeout["emit_retries"],
        )

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
        emit_with_retry(
            self.socket,
            "room:step:set",
            value,
            retries=self.timeout["emit_retries"],
        )

    @property
    def figure(self) -> str:
        self._delay_socket()
        return call_with_retry(
            self.socket,
            "analysis:figure:get",
            retries=self.timeout["call_retries"],
        )

    @figure.setter
    def figure(self, fig: str) -> None:
        emit_with_retry(
            self.socket,
            "analysis:figure:set",
            fig,
            retries=self.timeout["emit_retries"],
        )

    @property
    def atoms(self) -> ase.Atoms:
        return self[self.step]

    @atoms.setter
    def atoms(self, atoms: ase.Atoms) -> None:
        self[[self.step]] = [atoms]

    @property
    def points(self) -> np.ndarray:
        self._delay_socket()
        return np.array(
            call_with_retry(
                self.socket,
                "room:points:get",
                retries=self.timeout["call_retries"],
            )["0"]
        )

    @points.setter
    def points(self, points: np.ndarray) -> None:
        emit_with_retry(
            self.socket,
            "room:points:set",
            {"0": points.tolist()},
            retries=self.timeout["emit_retries"],
        )

    @property
    def bookmarks(self) -> dict[int, str]:
        self._delay_socket()

        return {
            int(k): v
            for k, v in call_with_retry(
                self.socket,
                "room:bookmarks:get",
                retries=self.timeout["call_retries"],
            ).items()
        }

    @bookmarks.setter
    def bookmarks(self, value: dict):
        emit_with_retry(
            self.socket,
            "room:bookmarks:set",
            value,
            retries=self.timeout["emit_retries"],
        )

    @property
    def camera(self) -> dict:
        self._delay_socket()

        return call_with_retry(
            self.socket,
            "room:camera:get",
            retries=self.timeout["call_retries"],
        )

    @camera.setter
    def camera(self, value: dict):
        raise NotImplementedError("Changing camera position is currently not possible")
        self.socket.emit("room:camera:set", value)

    @property
    def geometries(self) -> list[Object3D]:
        self._delay_socket()

        return [
            Geometry(method=x).method
            for x in call_with_retry(
                self.socket,
                "room:geometry:get",
                retries=self.timeout["call_retries"],
            )
        ]

    @geometries.setter
    def geometries(self, value: list[Object3D]):
        emit_with_retry(
            self.socket,
            "room:geometry:set",
            [x.model_dump() for x in value],
            retries=self.timeout["emit_retries"],
        )

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
        vis.timeout = self.timeout

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
