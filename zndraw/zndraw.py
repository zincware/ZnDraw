import dataclasses
import datetime
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
import typing_extensions as tyex
import znjson
import znsocket
from redis import Redis

from zndraw.abc import Message
from zndraw.base import Extension, ZnDrawBase
from zndraw.bonds import ASEComputeBonds
from zndraw.config import ArrowsConfig, ZnDrawConfig
from zndraw.converter import Object3DConverter
from zndraw.draw import Object3D
from zndraw.figure import Figure, FigureConverter
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
    r: Redis | znsocket.Client | None = None

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
    name: str | None = None

    _modifiers: dict[str, RegisterModifier] = dataclasses.field(default_factory=dict)
    _available: bool = True
    _refresh_client: znsocket.Client | None = None

    bond_calculator: ASEComputeBonds | None = dataclasses.field(
        default_factory=ASEComputeBonds, repr=False
    )

    def __post_init__(self):
        if self.socket is None:
            http_session = requests.Session()
            http_session.verify = self.verify
            self.socket = socketio.Client(http_session=http_session)

        self._refresh_client = znsocket.Client.from_url(self.url)
        if self.r is None:
            self.r = self._refresh_client

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

    def __getitem__(self, index: int | list | slice) -> ase.Atoms | list[ase.Atoms]:
        single_item = isinstance(index, int)
        if single_item:
            index = [index]
        if self.r.exists(f"room:{self.token}:frames"):
            structures = znsocket.List(
                self.r,
                f"room:{self.token}:frames",
                converter=[ASEConverter],
                socket=self._refresh_client,
            )[index]

        else:
            structures = znsocket.List(
                self.r,
                "room:default:frames",
                converter=[ASEConverter],
                socket=self._refresh_client,
            )[index]

        if single_item:
            return structures[0]
        return structures

    def __setitem__(
        self,
        index: int | list | slice,
        value: ase.Atoms | list[ase.Atoms],
    ):
        if isinstance(value, list):
            if not all(isinstance(x, ase.Atoms) for x in value):
                raise ValueError("Unable to parse provided data object")
        else:
            if not isinstance(value, ase.Atoms):
                raise ValueError("Unable to parse provided data object")
        lst = znsocket.List(
            self.r,
            f"room:{self.token}:frames",
            converter=[ASEConverter],
            socket=self._refresh_client,
        )
        if not self.r.exists(f"room:{self.token}:frames") and self.r.exists(
            "room:default:frames"
        ):
            default_lst = znsocket.List(
                self.r, "room:default:frames", socket=self._refresh_client
            )
            lst.copy(key=default_lst.key)

        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))
        if isinstance(index, int):
            index = [index]
            value = [value]

        for i, val in zip(index, value):
            if isinstance(val, ase.Atoms):
                if not hasattr(val, "connectivity") and self.bond_calculator is not None:
                    val.connectivity = self.bond_calculator.get_bonds(val)
            lst[i] = val

    def __len__(self) -> int:
        # TODO: what if the room does not exist yet?
        if not self.r.exists(f"room:{self.token}:frames"):
            return len(
                znsocket.List(self.r, "room:default:frames", socket=self._refresh_client)
            )
        return len(
            znsocket.List(
                self.r, f"room:{self.token}:frames", socket=self._refresh_client
            )
        )

    def __delitem__(self, index: int | slice | list[int]):
        lst = znsocket.List(
            self.r,
            f"room:{self.token}:frames",
            converter=[ASEConverter],
            socket=self._refresh_client,
        )
        if not self.r.exists(f"room:{self.token}:frames") and self.r.exists(
            "room:default:frames"
        ):
            default_lst = znsocket.List(
                self.r,
                "room:default:frames",
                converter=[ASEConverter],
                socket=self._refresh_client,
            )
            default_lst.copy(key=lst.key)

        del lst[index]

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

    def insert(self, index: int, value: ase.Atoms):
        if not isinstance(value, ase.Atoms):
            raise ValueError("Unable to parse provided data object")
        lst = znsocket.List(
            self.r,
            f"room:{self.token}:frames",
            converter=[ASEConverter],
            socket=self._refresh_client,
        )
        if not self.r.exists(f"room:{self.token}:frames") and self.r.exists(
            "room:default:frames"
        ):
            default_lst = znsocket.List(
                self.r,
                "room:default:frames",
                converter=[ASEConverter],
                socket=self._refresh_client,
            )
            default_lst.copy(key=lst.key)

        if isinstance(value, ase.Atoms):
            if not hasattr(value, "connectivity") and self.bond_calculator is not None:
                value.connectivity = self.bond_calculator.get_bonds(value)
        lst.insert(index, value)

    def extend(self, values: list[ase.Atoms]):
        if not isinstance(values, list) or not all(
            isinstance(x, ase.Atoms) for x in values
        ):
            raise ValueError("Unable to parse provided data object")

        # enable tbar if more than 10 messages are sent
        # approximated by the size of the first frame
        lst = znsocket.List(
            self.r,
            f"room:{self.token}:frames",
            converter=[ASEConverter],
            socket=self._refresh_client,
        )
        # TODO: why is there no copy action here?
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

                msg.append(val)
            else:
                msg.append(val)
            if (
                len(
                    json.dumps(
                        msg, cls=znjson.ZnEncoder.from_converters([ASEConverter])
                    ).encode("utf-8")
                )
                > self.maximum_message_size
            ):
                lst.extend(msg)
                msg = []
        if len(msg) > 0:  # Only send the message if it's not empty
            lst.extend(msg)

    @property
    def selection(self) -> list[int]:
        try:
            return znsocket.Dict(
                self.r, f"room:{self.token}:selection", socket=self._refresh_client
            )["grp-0"]
        except KeyError:
            return []

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
        znsocket.Dict(
            self.r, f"room:{self.token}:selection", socket=self._refresh_client
        )["grp-0"] = value

    @property
    def step(self) -> int:
        try:
            return znsocket.Dict(
                self.r, f"room:{self.token}:step", socket=self._refresh_client
            )["grp-0"]
        except KeyError:
            return 0

    def log(self, message: str) -> None:
        msg: Message = {
            "time": datetime.datetime.now().isoformat(),
            "msg": message,
            "origin": self.name,
        }
        znsocket.List(
            self.r, f"room:{self.token}:chat", socket=self._refresh_client
        ).append(msg)

    @property
    def messages(self) -> list[Message]:
        return znsocket.List(
            self.r,
            f"room:{self.token}:chat",
            repr_type="length",
            socket=self._refresh_client,
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
        znsocket.Dict(self.r, f"room:{self.token}:step", socket=self._refresh_client)[
            "grp-0"
        ] = value

    @property
    def figures(self) -> dict[str, go.Figure | Figure]:
        return znsocket.Dict(
            self.r,
            f"room:{self.token}:figures",
            repr_type="full",
            socket=self._refresh_client,
            converter=[znjson.converter.PlotlyConverter, FigureConverter],
        )

    @figures.setter
    def figures(self, data: dict[str, go.Figure]) -> None:
        """Update the figures on the remote."""
        figures_dict = znsocket.Dict(
            self.r,
            f"room:{self.token}:figures",
            socket=self._refresh_client,
            converter=[znjson.converter.PlotlyConverter, FigureConverter],
        )
        figures_dict.clear()
        figures_dict.update(data)

    @property
    def atoms(self) -> ase.Atoms:
        return self[self.step]

    @atoms.setter
    def atoms(self, atoms: ATOMS_LIKE) -> None:
        self[[self.step]] = [atoms]

    @property
    def points(self) -> np.ndarray:
        # TODO: use znsocket.List inside znsocket.Dict here
        try:
            return np.array(
                znsocket.Dict(
                    self.r, f"room:{self.token}:points", socket=self._refresh_client
                )["grp-0"]
            )
        except KeyError:
            return np.array([])

    @points.setter
    def points(self, points: np.ndarray | list) -> None:
        if isinstance(points, np.ndarray):
            points = points.tolist()
        znsocket.Dict(self.r, f"room:{self.token}:points", socket=self._refresh_client)[
            "grp-0"
        ] = points

    @property
    def bookmarks(self) -> dict[int, str]:
        return znsocket.Dict(
            self.r,
            f"room:{self.token}:bookmarks",
            repr_type="full",
            socket=self._refresh_client,
        )

    @bookmarks.setter
    def bookmarks(self, value: dict[int, str]):
        if not isinstance(value, dict):
            raise ValueError("Bookmarks must be a dictionary")
        if not all(isinstance(x, int) for x in value.keys()):
            raise ValueError("Bookmark keys must be integers")
        if not all(isinstance(x, str) for x in value.values()):
            raise ValueError("Bookmark values must be strings")

        bookmarks = znsocket.Dict(
            self.r, f"room:{self.token}:bookmarks", socket=self._refresh_client
        )
        bookmarks.clear()
        bookmarks.update(value)

    @property
    def camera(self) -> CameraData:
        camera_dct = znsocket.Dict(
            self.r,
            f"room:{self.token}:camera",
            repr_type="full",
            socket=self._refresh_client,
        )
        if "position" not in camera_dct:
            camera_dct["position"] = [0, 0, 0]
        if "target" not in camera_dct:
            camera_dct["target"] = [0, 0, 0]
        return camera_dct

    @camera.setter
    def camera(self, value: CameraData):
        if set(value.keys()) != {"position", "target"}:
            raise ValueError("Camera must have 'position' and 'target' keys")
        znsocket.Dict(
            self.r, f"room:{self.token}:camera", socket=self._refresh_client
        ).update(value)

    @property
    def geometries(self) -> list[Object3D]:
        return znsocket.List(
            self.r,
            f"room:{self.token}:geometries",
            repr_type="full",
            socket=self._refresh_client,
            converter=[Object3DConverter],
        )

    @geometries.setter
    def geometries(self, value: list[Object3D]):
        if not all(isinstance(x, Object3D) for x in value):
            raise ValueError("Geometries must be a list of Object3D instances")
        lst = znsocket.List(
            self.r,
            f"room:{self.token}:geometries",
            socket=self._refresh_client,
            converter=[Object3DConverter],
        )
        lst.clear()
        lst.extend(value)

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


@tyex.deprecated("Use ZnDraw instead.")
@dataclasses.dataclass(kw_only=True)
class ZnDrawLocal(ZnDraw):
    """All the functionality is provided by the ZnDraw class."""
