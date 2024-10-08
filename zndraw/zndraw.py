import dataclasses
import importlib.metadata
import json
import logging
import typing as t

import ase
import numpy as np
import plotly.graph_objects as go
import requests
import socketio.exceptions
import tqdm
import znjson
import znsocket
from plotly.io import from_json as ploty_from_json
from redis import Redis

from zndraw.base import Extension, ZnDrawBase
from zndraw.bonds import ASEComputeBonds
from zndraw.config import ArrowsConfig, ZnDrawConfig
from zndraw.draw import Geometry, Object3D
from zndraw.exceptions import RoomLockedError
from zndraw.scene import Scene
from zndraw.type_defs import (
    ATOMS_LIKE,
    CameraData,
    JupyterConfig,
    RegisterModifier,
    TimeoutConfig,
)
from zndraw.utils import (
    ASEConverter,
    call_with_retry,
    convert_url_to_http,
    emit_with_retry,
    parse_url,
)

log = logging.getLogger(__name__)
__version__ = importlib.metadata.version("zndraw")


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


def _check_version_compatibility(server_version: str) -> None:
    if server_version != __version__:
        log.warning(
            f"Server version ({server_version}) and client version ({__version__}) are not the same. "
            "This may lead to unexpected behavior."
        )


@dataclasses.dataclass
class ZnDraw(ZnDrawBase):
    url: str
    token: str | None = None
    auth_token: str | None = None

    socket: socketio.Client | None = dataclasses.field(default=None, repr=False)
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
    jupyter_config: JupyterConfig = dataclasses.field(
        default_factory=lambda: JupyterConfig(
            width="70%",
            height=600,
        )
    )
    verify: bool | str = True

    maximum_message_size: int = dataclasses.field(default=500_000, repr=False)

    _modifiers: dict[str, RegisterModifier] = dataclasses.field(default_factory=dict)
    _available: bool = True

    bond_calculator: ASEComputeBonds | None = dataclasses.field(
        default_factory=ASEComputeBonds, repr=False
    )

    def __post_init__(self):
        if self.socket is None:
            http_session = requests.Session()
            http_session.verify = self.verify
            self.socket = socketio.Client(http_session=http_session)

        def on_wakeup():
            if self._available:
                self.socket.emit("modifier:available", list(self._modifiers))

        self.url = self.url.replace("http", "ws")
        self.socket.on("connect", self._on_connect)
        self.socket.on("modifier:run", self._run_modifier)
        self.socket.on("modifier:wakeup", on_wakeup)
        self.socket.on("room:log", lambda x: print(x))
        self.socket.on("version", _check_version_compatibility)

        for idx in range(self.timeout["connect_retries"] + 1):
            try:
                _url, _path = parse_url(self.url)
                if _path:
                    self.socket.connect(
                        _url,
                        wait_timeout=self.timeout["connection"],
                        socketio_path=f"/{_path}/socket.io",
                    )
                else:
                    self.socket.connect(_url, wait_timeout=self.timeout["connection"])
                break
            except socketio.exceptions.ConnectionError as err:
                log.warning("Connection failed. Retrying...")
                self.socket.sleep(self.timeout["connection"])
                if idx == self.timeout["connect_retries"]:
                    raise socketio.exceptions.ConnectionError(
                        f"Unable to connect to ZnDraw server at '{self.url}'. Is the server running?"
                    ) from err

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

    def __getitem__(self, index) -> ase.Atoms | list[ase.Atoms]:
        single_item = isinstance(index, int)
        if single_item:
            index = [index]

        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))

        # make negative indices count from the end
        index = [i if i >= 0 else len(self) + i for i in index]

        if any(x >= len(self) for x in index):
            raise IndexError("Index out of range")

        data = call_with_retry(
            self.socket,
            "room:frames:get",
            index,
            retries=self.timeout["call_retries"],
        )

        structures = [
            znjson.loads(
                json.dumps(x), cls=znjson.ZnDecoder.from_converters([ASEConverter])
            )
            for x in data.values()
        ]
        return structures[0] if single_item else structures

    def __setitem__(
        self, index: int | list[int] | slice, value: ATOMS_LIKE | list[ATOMS_LIKE]
    ):
        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        if isinstance(index, int):
            index = [index]
            value = [value]

        data = {}
        for i, val in zip(index, value):
            if isinstance(val, ase.Atoms):
                if not hasattr(val, "connectivity") and self.bond_calculator is not None:
                    val.connectivity = self.bond_calculator.get_bonds(val)
                data[i] = znjson.dumps(
                    val, cls=znjson.ZnEncoder.from_converters([ASEConverter])
                )
            else:
                data[i] = val
            if '"_type": "ase.Atoms"' not in data[i]:
                raise ValueError("Unable to parse provided data object")
        call_with_retry(
            self.socket,
            "room:frames:set",
            data,
            retries=self.timeout["call_retries"],
        )

    def __len__(self) -> int:
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

        call_with_retry(
            self.socket,
            "room:frames:delete",
            index,
            retries=self.timeout["call_retries"],
        )

    def _repr_html_(self):
        from IPython.display import IFrame

        # TODO: save address and do not replace in post_init
        address = convert_url_to_http(f"{self.url}/token/{self.token}")
        log.info(f"Opening ZnDraw at {address}")
        return IFrame(
            address,
            width=self.jupyter_config["width"],
            height=self.jupyter_config["height"],
        )._repr_html_()

    def insert(self, index: int, value: ATOMS_LIKE):
        if isinstance(value, ase.Atoms):
            if not hasattr(value, "connectivity") and self.bond_calculator is not None:
                value.connectivity = self.bond_calculator.get_bonds(value)

            value = znjson.dumps(
                value, cls=znjson.ZnEncoder.from_converters([ASEConverter])
            )

        if '"_type": "ase.Atoms"' not in value:
            raise ValueError("Unable to parse provided data object")

        call_with_retry(
            self.socket,
            "room:frames:insert",
            {
                "index": index,
                "value": value,
            },
            retries=self.timeout["call_retries"],
        )

    def extend(self, values: list[ATOMS_LIKE]):
        msg = {}

        # enable tbar if more than 10 messages are sent
        # approximated by the size of the first frame

        if not isinstance(values, list):
            raise ValueError("Unable to parse provided data object")
        if self.locked:
            raise RoomLockedError("The room you are trying to modify is locked.")

        show_tbar = (
            len(values)
            * len(
                znjson.dumps(
                    values[0], cls=znjson.ZnEncoder.from_converters([ASEConverter])
                ).encode("utf-8")
            )
        ) > (10 * self.maximum_message_size)
        tbar = tqdm.tqdm(
            values, desc="Sending frames", unit=" frame", disable=not show_tbar
        )

        for i, val in enumerate(tbar, start=len(self)):
            if isinstance(val, ase.Atoms):
                if not hasattr(val, "connectivity") and self.bond_calculator is not None:
                    val.connectivity = self.bond_calculator.get_bonds(val)

                msg[i] = znjson.dumps(
                    val, cls=znjson.ZnEncoder.from_converters([ASEConverter])
                )
            else:
                msg[i] = val
            if '"_type": "ase.Atoms"' not in msg[i]:
                raise ValueError("Unable to parse provided data object")
            if len(json.dumps(msg).encode("utf-8")) > self.maximum_message_size:
                call_with_retry(
                    self.socket,
                    "room:frames:set",
                    msg,
                    retries=self.timeout["call_retries"],
                )
                msg = {}
                # after each large message, wait a bit
                self.socket.sleep(self.timeout["modifier"])
        if len(msg) > 0:  # Only send the message if it's not empty
            call_with_retry(
                self.socket,
                "room:frames:set",
                msg,
                retries=self.timeout["call_retries"],
            )

    @property
    def selection(self) -> list[int]:
        return call_with_retry(
            self.socket,
            "room:selection:get",
            retries=self.timeout["call_retries"],
        )["0"]

    @selection.setter
    def selection(self, value: list[int]):
        if not isinstance(value, list):
            raise ValueError("Selection must be a list")
        if not all(isinstance(x, int) for x in value):
            raise ValueError("Selection must be a list of integers")
        if len(value) != len(set(value)):
            raise ValueError("Selection must not contain duplicates")

        max_index = len(self.atoms)
        if any(x >= max_index for x in value):
            raise IndexError("Selection out of range")
        if any(x < 0 for x in value):
            raise IndexError("Selection must be positive")
        emit_with_retry(
            self.socket,
            "room:selection:set",
            {"0": value},
            retries=self.timeout["emit_retries"],
        )

    @property
    def step(self) -> int:
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
        if not isinstance(value, int):
            raise ValueError("Step must be an integer")
        if value < 0:
            raise ValueError("Step must be positive")
        if value >= len(self):
            raise IndexError("Step out of range")
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
    def figures(self) -> dict[str, go.Figure]:
        # TODO: znjson.loads
        data = call_with_retry(
            self.socket,
            "analysis:figure:get",
            retries=self.timeout["call_retries"],
        )
        return {k: ploty_from_json(v) for k, v in data.items()}

    @figures.setter
    def figures(self, data: dict[str, go.Figure]) -> None:
        """Update the figures on the remote."""
        # TODO: can you use znsocket.Dict
        # to update the data an avoid
        # sending duplicates?
        data = {k: v.to_json() for k, v in data.items()}
        emit_with_retry(
            self.socket,
            "analysis:figure:set",
            data,
            retries=self.timeout["emit_retries"],
        )

    @property
    def atoms(self) -> ase.Atoms:
        return self[self.step]

    @atoms.setter
    def atoms(self, atoms: ATOMS_LIKE) -> None:
        self[[self.step]] = [atoms]

    @property
    def points(self) -> np.ndarray:
        return np.array(
            call_with_retry(
                self.socket,
                "room:points:get",
                retries=self.timeout["call_retries"],
            )["0"]
        )

    @points.setter
    def points(self, points: np.ndarray | list) -> None:
        if isinstance(points, list):
            points = np.array(points)
        emit_with_retry(
            self.socket,
            "room:points:set",
            {"0": points.tolist()},
            retries=self.timeout["emit_retries"],
        )

    @property
    def bookmarks(self) -> dict[int, str]:
        return {
            int(k): v
            for k, v in call_with_retry(
                self.socket,
                "room:bookmarks:get",
                retries=self.timeout["call_retries"],
            ).items()
        }

    @bookmarks.setter
    def bookmarks(self, value: dict[int, str]):
        if not isinstance(value, dict):
            raise ValueError("Bookmarks must be a dictionary")
        if not all(isinstance(x, int) for x in value.keys()):
            raise ValueError("Bookmark keys must be integers")
        if not all(isinstance(x, str) for x in value.values()):
            raise ValueError("Bookmark values must be strings")

        emit_with_retry(
            self.socket,
            "room:bookmarks:set",
            value,
            retries=self.timeout["emit_retries"],
        )

    @property
    def camera(self) -> CameraData:
        return call_with_retry(
            self.socket,
            "room:camera:get",
            retries=self.timeout["call_retries"],
        )

    @camera.setter
    def camera(self, value: CameraData):
        if set(value.keys()) != {"position", "target"}:
            raise ValueError("Camera must have 'position' and 'target' keys")
        self.socket.emit("room:camera:set", {"content": value, "emit": True})

    @property
    def geometries(self) -> list[Object3D]:
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
        if not isinstance(value, list):
            raise ValueError("Geometries must be a list")
        if not all(isinstance(x, Object3D) for x in value):
            raise ValueError("Geometries must be a list of Object3D instances")
        emit_with_retry(
            self.socket,
            "room:geometry:set",
            [x.model_dump() for x in value],
            retries=self.timeout["emit_retries"],
        )

    @property
    def config(self) -> ZnDrawConfig:
        config: dict = call_with_retry(
            self.socket,
            "room:config:get",
            retries=self.timeout["call_retries"],
        )
        return ZnDrawConfig(
            vis=self,
            arrows=ArrowsConfig(**config.get("arrows", {})),
            scene=Scene(**config.get("scene", {})),
        )

    @property
    def locked(self) -> bool:
        return call_with_retry(self.socket, "room:lock:get")

    @locked.setter
    def locked(self, value: bool) -> None:
        emit_with_retry(self.socket, "room:lock:set", value)

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


@dataclasses.dataclass(kw_only=True)
class ZnDrawLocal(ZnDraw):
    """Access database directly.

    This client can / should be used if the database can be accessed directly.
    Data will not be loaded via sockets but modified directly in the database.
    """

    r: Redis

    def __getitem__(self, index: int | list | slice) -> ase.Atoms | list[ase.Atoms]:
        single_item = isinstance(index, int)
        if single_item:
            index = [index]
        if self.r.exists(f"room:{self.token}:frames"):
            data = znsocket.List(self.r, f"room:{self.token}:frames")[index]

        else:
            try:
                data = znsocket.List(self.r, "room:default:frames")[index]
            except IndexError:
                data = []

        structures = [
            znjson.loads(x, cls=znjson.ZnDecoder.from_converters([ASEConverter]))
            for x in data
        ]
        if single_item:
            return structures[0]
        return structures

    def insert(self, index: int, value: ATOMS_LIKE):
        lst = znsocket.List(self.r, f"room:{self.token}:frames")
        if not self.r.exists(f"room:{self.token}:frames"):
            default_lst = znsocket.List(self.r, "room:default:frames")
            # TODO: using a redis copy action would be faster
            lst.extend(default_lst)

        if isinstance(value, ase.Atoms):
            if not hasattr(value, "connectivity") and self.bond_calculator is not None:
                value.connectivity = self.bond_calculator.get_bonds(value)

            value = znjson.dumps(
                value, cls=znjson.ZnEncoder.from_converters([ASEConverter])
            )

        if '"_type": "ase.Atoms"' not in value:
            raise ValueError("Unable to parse provided data object")
        lst.insert(index, value)
        self.socket.emit("room:frames:refresh", [self.step])

    def extend(self, values: list[ATOMS_LIKE]):
        if not isinstance(values, list):
            raise ValueError("Unable to parse provided data object")

        # enable tbar if more than 10 messages are sent
        # approximated by the size of the first frame
        lst = znsocket.List(self.r, f"room:{self.token}:frames")
        show_tbar = (
            len(values)
            * len(
                znjson.dumps(
                    values[0], cls=znjson.ZnEncoder.from_converters([ASEConverter])
                ).encode("utf-8")
            )
        ) > (10 * self.maximum_message_size)
        tbar = tqdm.tqdm(
            values, desc="Sending frames", unit=" frame", disable=not show_tbar
        )

        msg = []

        for val in tbar:
            if isinstance(val, ase.Atoms):
                if not hasattr(val, "connectivity") and self.bond_calculator is not None:
                    val.connectivity = self.bond_calculator.get_bonds(val)

                msg.append(
                    znjson.dumps(
                        val, cls=znjson.ZnEncoder.from_converters([ASEConverter])
                    )
                )
            else:
                msg.append(val)
            if '"_type": "ase.Atoms"' not in msg[-1]:
                raise ValueError("Unable to parse provided data object")
            if len(json.dumps(msg).encode("utf-8")) > self.maximum_message_size:
                lst.extend(msg)
                msg = []
        if len(msg) > 0:  # Only send the message if it's not empty
            lst.extend(msg)

        self.socket.emit("room:frames:refresh", [self.step])

    def __setitem__(
        self,
        index: int | list | slice,
        value: ATOMS_LIKE | list[ATOMS_LIKE],
    ):
        lst = znsocket.List(self.r, f"room:{self.token}:frames")
        if not self.r.exists(f"room:{self.token}:frames"):
            default_lst = znsocket.List(self.r, "room:default:frames")
            # TODO: using a redis copy action would be faster
            lst.extend(default_lst)

        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        if isinstance(index, int):
            index = [index]
            value = [value]

        for i, val in zip(index, value):
            if isinstance(val, ase.Atoms):
                if not hasattr(val, "connectivity") and self.bond_calculator is not None:
                    val.connectivity = self.bond_calculator.get_bonds(val)
                val = znjson.dumps(
                    val, cls=znjson.ZnEncoder.from_converters([ASEConverter])
                )
            if '"_type": "ase.Atoms"' not in val:
                raise ValueError("Unable to parse provided data object")
            lst[i] = val
        self.socket.emit("room:frames:refresh", [self.step])
