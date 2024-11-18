import dataclasses
import datetime
import enum
import importlib.metadata
import logging
import typing as t
from collections.abc import MutableSequence

import ase
import numpy as np
import plotly.graph_objects as go
import requests
import socketio.exceptions
import splines
import tqdm
import typing_extensions as tyex
import znjson
import znsocket
import znsocket.exceptions
from redis import Redis

from zndraw.abc import Message
from zndraw.base import Extension
from zndraw.bonds import ASEComputeBonds
from zndraw.config import SETTINGS
from zndraw.converter import ASEConverter, Object3DConverter
from zndraw.draw import Object3D
from zndraw.figure import Figure, FigureConverter
from zndraw.queue import check_queue
from zndraw.type_defs import (
    ATOMS_LIKE,
    CameraData,
    JupyterConfig,
    RegisterModifier,
    TimeoutConfig,
)
from zndraw.utils import (
    call_with_retry,
    convert_url_to_http,
    emit_with_retry,
    parse_url,
)

log = logging.getLogger(__name__)
__version__ = importlib.metadata.version("zndraw")


class ExtensionType(str, enum.Enum):
    """The type of the extension."""

    MODIFIER = "modifier"
    SELECTION = "selection"
    ANALYSIS = "analysis"


def _check_version_compatibility(server_version: str) -> None:
    if server_version != __version__:
        log.warning(
            f"Server version ({server_version}) and client version ({__version__}) are not the same. "
            "This may lead to unexpected behavior."
        )


@dataclasses.dataclass
class ZnDraw(MutableSequence):
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
    convert_nan: bool = False

    max_atoms_per_call: int = dataclasses.field(default=1000, repr=False)
    # number of `ase.Atom` to send per call
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
        # TODO: the refresh_client should be able to use the same socket connection!
        if self.r is None:
            self.r = self._refresh_client

        self.url = self.url.replace("http", "ws")
        self.socket.on("connect", self._on_connect)
        self.socket.on("room:log", lambda x: print(x))
        self.socket.on("version", _check_version_compatibility)
        self.socket.on("disconnect", lambda: print("Disconnected from ZnDraw server"))

        self._connect()

        self.socket.start_background_task(check_queue, self)

    def _connect(self):
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

        # self.socket.sleep(5)  # wait for znsocket to reconnect as well
        registerd_modifiers = list(self._modifiers[x]["cls"] for x in self._modifiers)
        for modifier in registerd_modifiers:
            try:
                self.register_modifier(
                    modifier, public=self._modifiers[modifier.__name__]["public"]
                )
            except (
                socketio.exceptions.BadNamespaceError
            ) as err:  # /znsocket is not a connected namespace.
                log.error(str(err))
                log.error("Disconnecting. Please restart.")
                self.socket.disconnect()

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
                max_commands_per_call=100,
            )[index]

        else:
            structures = znsocket.List(
                self.r,
                "room:default:frames",
                converter=[ASEConverter],
                socket=self._refresh_client,
                max_commands_per_call=100,
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
            max_commands_per_call=100,
            convert_nan=self.convert_nan,
        )
        if not self.r.exists(f"room:{self.token}:frames") and self.r.exists(
            "room:default:frames"
        ):
            default_lst = znsocket.List(
                self.r,
                "room:default:frames",
                socket=self._refresh_client,
                max_commands_per_call=100,
            )
            if not (len(default_lst) == 1 and len(default_lst[0]) == 0):
                # prevent copying empty default room
                default_lst.copy(key=lst.key)

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
                znsocket.List(
                    self.r,
                    "room:default:frames",
                    socket=self._refresh_client,
                )
            )
        return len(
            znsocket.List(
                self.r,
                f"room:{self.token}:frames",
                socket=self._refresh_client,
            )
        )

    def __delitem__(self, index: int | slice | list[int]):
        lst = znsocket.List(
            self.r,
            f"room:{self.token}:frames",
            converter=[ASEConverter],
            socket=self._refresh_client,
            max_commands_per_call=100,
        )

        if not self.r.exists(f"room:{self.token}:frames") and self.r.exists(
            "room:default:frames"
        ):
            default_lst = znsocket.List(
                self.r,
                "room:default:frames",
                converter=[ASEConverter],
                socket=self._refresh_client,
                max_commands_per_call=100,
            )
            if not (len(default_lst) == 1 and len(default_lst[0]) == 0):
                # prevent copying empty default room
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
            convert_nan=self.convert_nan,
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
            if not (len(default_lst) == 1 and len(default_lst[0]) == 0):
                # prevent copying empty default room
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

        if not self.r.exists(f"room:{self.token}:frames") and self.r.exists(
            "room:default:frames"
        ):
            default_lst = znsocket.List(
                self.r,
                "room:default:frames",
                converter=[ASEConverter],
                socket=self._refresh_client,
                max_commands_per_call=100,
            )
            if not (len(default_lst) == 1 and len(default_lst[0]) == 0):
                # prevent copying empty default room
                default_lst.copy(key=f"room:{self.token}:frames")

        # enable tbar if more than 10 messages are sent
        # approximated by the size of the first frame
        lst = znsocket.List(
            self.r,
            f"room:{self.token}:frames",
            converter=[ASEConverter],
            socket=self._refresh_client,
            max_commands_per_call=100,
            convert_nan=self.convert_nan,
        )
        # TODO: why is there no copy action here?
        show_tbar = (len(values[0]) * len(values)) > self.max_atoms_per_call
        tbar = tqdm.tqdm(
            values, desc="Sending frames", unit=" frame", disable=not show_tbar
        )

        msg = []
        n_atoms = 0

        for val in tbar:
            if not hasattr(val, "connectivity") and self.bond_calculator is not None:
                val.connectivity = self.bond_calculator.get_bonds(val)

            msg.append(val)
            n_atoms += len(val)

            if n_atoms > self.max_atoms_per_call:
                lst.extend(msg)
                msg = []
                n_atoms = 0
        if len(msg) > 0:
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
            self.r,
            f"room:{self.token}:chat",
            socket=self._refresh_client,
            max_commands_per_call=100,
            convert_nan=self.convert_nan,
        ).append(msg)

    @property
    def messages(self) -> list[Message]:
        return znsocket.List(
            self.r,
            f"room:{self.token}:chat",
            repr_type="length",
            socket=self._refresh_client,
            max_commands_per_call=100,
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
        try:
            return self[self.step]
        except IndexError:
            return ase.Atoms()

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
    def segments(self) -> np.ndarray:
        points = self.points
        if points.shape[0] <= 1:
            return points
        t = np.linspace(0, len(points) - 1, len(points) * 50)
        return splines.CatmullRom(points).evaluate(t)

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
            max_commands_per_call=100,
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
            max_commands_per_call=100,
        )
        lst.clear()
        lst.extend(value)

    @property
    def config(self) -> znsocket.Dict:
        conf = znsocket.Dict(
            self.r,
            f"room:{self.token}:config",
            repr_type="full",
            socket=self._refresh_client,
        )
        if len(conf) == 0:
            conf_dicts = {}
            for key, value in SETTINGS.items():
                conf_dicts[key] = znsocket.Dict(
                    self.r,
                    f"room:{self.token}:config:{key}",
                    repr_type="full",
                    socket=self._refresh_client,
                )
                if len(conf_dicts[key]) == 0:
                    conf_dicts[key].update(value().model_dump())
            conf.update(conf_dicts)
            # TODO: arrows
        return conf

    @property
    def locked(self) -> bool:
        return call_with_retry(self.socket, "room:lock:get")

    @locked.setter
    def locked(self, value: bool) -> None:
        emit_with_retry(self.socket, "room:lock:set", value)

    def register(
        self,
        cls: t.Type[Extension],
        run_kwargs: dict | None = None,
        public: bool = False,
        variant: ExtensionType = ExtensionType.MODIFIER,
    ):
        """Register an extension class."""
        if run_kwargs is None:
            run_kwargs = {}

        if variant == ExtensionType.MODIFIER:
            # TODO: need to create a custom vis object for the room!
            if public:
                modifier_schema = znsocket.Dict(
                    self.r,
                    "schema:default:modifier",
                    socket=self._refresh_client,
                )
            else:
                modifier_schema = znsocket.Dict(
                    self.r,
                    f"schema:{self.token}:modifier",
                    socket=self._refresh_client,
                )
            # TODO: check if the key exists and if it is different?
            # TODO: also check in the default schema room.
            # TODO: can not register the same modifier twice!
            modifier_schema[cls.__name__] = cls.model_json_schema()

            # TODO: if public!
            if public:
                modifier_registry = znsocket.Dict(
                    self.r,
                    "registry:default:modifier",
                    socket=self._refresh_client,
                )
            else:
                modifier_registry = znsocket.Dict(
                    self.r,
                    f"registry:{self.token}:modifier",
                    socket=self._refresh_client,
                )
            if self.socket.get_sid() not in modifier_registry:
                modifier_registry[self.socket.get_sid()] = [cls.__name__]
            else:
                modifier_registry[self.socket.get_sid()] = modifier_registry[
                    self.socket.get_sid()
                ] + [cls.__name__]

            self._modifiers[cls.__name__] = {
                "cls": cls,
                "run_kwargs": run_kwargs,
                "public": public,
            }

        else:
            raise NotImplementedError(f"Variant {variant} is not implemented")

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
        if run_kwargs is None:
            run_kwargs = {}

        # if use_frozen:
        #     cls = cls.frozen()

        self.register(
            cls=cls,
            run_kwargs=run_kwargs,
            public=public,
            variant=ExtensionType.MODIFIER,
        )


@tyex.deprecated("Use ZnDraw instead.")
@dataclasses.dataclass(kw_only=True)
class ZnDrawLocal(ZnDraw):
    """All the functionality is provided by the ZnDraw class."""
