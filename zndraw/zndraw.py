import contextlib
import dataclasses
import pathlib
import threading
import time
import typing as t
from io import StringIO

import ase
import ase.io
import numpy as np
import socketio
import tqdm
import znh5md

from zndraw.analyse import get_analysis_class
from zndraw.data import atoms_from_json, atoms_to_json
from zndraw.draw import Geometry
from zndraw.modify import get_modify_class
from zndraw.select import get_selection_class
from zndraw.settings import GlobalConfig
from zndraw.utils import ZnDrawLoggingHandler


def _await_answer(socket, channel, data=None, timeout=5):
    """Wait for an answer from the server.

    I haven't used asyncio, so this should do..
    """
    event = threading.Event()

    answer = None

    def on_answer(data):
        nonlocal answer
        answer = data
        event.set()

    socket.on(channel, on_answer)
    socket.emit(channel, data)

    event.wait(timeout=timeout)
    return answer


def process_line_data(data) -> tuple[np.ndarray, np.ndarray]:
    """Get the points of the selected atoms"""
    points = np.array([[val["x"], val["y"], val["z"]] for val in data["points"]])
    segments = np.array(data["segments"])

    return points, segments


@dataclasses.dataclass
class FileIO:
    name: str = None
    start: int = 0
    stop: int = None
    step: int = 1


@dataclasses.dataclass
class ZnDrawBase:  # collections.abc.MutableSequence
    url: str
    token: str

    _target_sid: str = None

    def __post_init__(self):
        self.socket = socketio.Client()
        self.socket.on("connect", lambda: self.socket.emit("join", self.token))
        self.socket.on("disconnect", lambda: self.socket.disconnect())

    def _connect(self):
        for _ in range(100):
            try:
                self.socket.connect(self.url)
                break
            except socketio.exceptions.ConnectionError:
                time.sleep(0.1)
        else:
            raise socketio.exceptions.ConnectionError

    def __len__(self) -> int:
        if self._target_sid is not None:
            return int(self.socket.call("atoms:length", {"sid": self._target_sid}))
        return int(self.socket.call("atoms:length", {}))

    def __setitem__(self, index, value):
        assert isinstance(value, ase.Atoms), "Must be an ASE Atoms object"
        assert isinstance(index, int), "Index must be an integer"

        data = {index: atoms_to_json(value)}
        if self._target_sid is not None:
            data["sid"] = self._target_sid

        self.socket.emit("atoms:upload", data)

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

            data = {"index": index}
            if self._target_sid is not None:
                data["sid"] = self._target_sid

            self.socket.emit("atoms:delete", data)
            if not is_slice and (index[0] >= length or index[-1] >= length):
                raise IndexError("Index out of range")
        else:
            raise TypeError("Index must be an integer, slice or list[int]")

    def insert(self, index, value):
        """Insert atoms before index"""
        self.socket.emit("atoms:insert", {index: atoms_to_json(value)})

    def append(self, value: ase.Atoms) -> None:
        """Append atoms to the end of the list"""
        self[len(self)] = value

    def extend(self, values: list[ase.Atoms]) -> None:
        """Extend the list by appending all the items in the given list"""
        for value in values:
            self.append(value)

    def __getitem__(self, index) -> t.Union[ase.Atoms, list[ase.Atoms]]:
        length = len(self)
        is_scalar = isinstance(index, int)
        is_sclice = isinstance(index, slice)
        if is_sclice:
            index = range(*index.indices(length))

        index = [index] if isinstance(index, int) else index
        index = [i if i >= 0 else length + i for i in index]

        data = {"indices": index}
        if self._target_sid is not None:
            data["sid"] = self._target_sid

        downloaded_data = self.socket.call("atoms:download", data)

        atoms_list = []

        for val in downloaded_data.values():
            atoms_list.append(atoms_from_json(val))

        data = atoms_list[0] if is_scalar else atoms_list
        if data == [] and not is_sclice:
            raise IndexError("Index out of range")
        return data

    @property
    def points(self) -> np.ndarray:
        if self._target_sid is not None:
            data = self.socket.call("scene:points", {"sid": self._target_sid})
        else:
            data = self.socket.call("scene:points", {})
        return np.array([[val["x"], val["y"], val["z"]] for val in data])

    @property
    def segments(self) -> np.ndarray:
        if self._target_sid is not None:
            data = self.socket.call("scene:segments", {"sid": self._target_sid})
        else:
            data = self.socket.call("scene:segments", {})
        return np.array(data)

    @property
    def step(self) -> int:
        if self._target_sid is not None:
            step = int(self.socket.call("scene:step", {"sid": self._target_sid}))
        else:
            step = int(self.socket.call("scene:step", {}))
        return step

    @step.setter
    def step(self, index):
        data = {"index": index}
        if self._target_sid is not None:
            data["sid"] = self._target_sid
        self.socket.emit("scene:set", data)

    @property
    def selection(self) -> list[int]:
        if self._target_sid is not None:
            return self.socket.call("selection:get", {"sid": self._target_sid})
        return self.socket.call("selection:get", {})

    @selection.setter
    def selection(self, value: list[int]):
        data = {"selection": value}
        if self._target_sid is not None:
            data["sid"] = self._target_sid
        self.socket.emit("selection:set", data)


@dataclasses.dataclass
class ZnDrawDefault(ZnDrawBase):
    file_io: FileIO = None

    def __post_init__(self):
        super().__post_init__()

        self.socket.on("webclient:available", self.initialize_webclient)
        self.socket.on("modifier:run", self.modifier_run)
        self.socket.on("selection:run", self.selection_run)
        self.socket.on("analysis:run", self.analysis_run)
        self.socket.on("upload", self.upload_file)
        self.socket.on("download:request", self.download_file)
        self._connect()
        self.socket.wait()

    def initialize_webclient(self, sid):
        with self._set_sid(sid):
            self.read_data()
            self.analysis_schema()
            self.modifier_schema()
            self.selection_schema()
            self.draw_schema()

    @contextlib.contextmanager
    def _set_sid(self, sid):
        self._target_sid = sid
        yield
        self._target_sid = None

    def read_data(self):
        if self.file_io.name is None:
            return

        if pathlib.Path(self.file_io.name).suffix == ".h5":
            reader = znh5md.ASEH5MD(self.file_io.name)
            generator = reader.get_atoms_list()
        else:
            generator = ase.io.iread(self.file_io.name)

        frame = 0

        for idx, atoms in tqdm.tqdm(enumerate(generator), ncols=100):
            if self.file_io.start and idx < self.file_io.start:
                continue
            if self.file_io.stop and idx >= self.file_io.stop:
                break
            if self.file_io.step and idx % self.file_io.step != 0:
                continue
            self[frame] = atoms
            self.step = frame
            frame += 1

    def analysis_schema(self):
        config = GlobalConfig.load()

        cls = get_analysis_class(config.get_analysis_methods())

        try:
            atoms = self[0]
        except IndexError:
            atoms = ase.Atoms()

        self.socket.emit(
            "analysis:schema",
            {
                "schema": cls.model_json_schema_from_atoms(atoms),
                "sid": self._target_sid,
            },
        )

    def modifier_schema(self):
        config = GlobalConfig.load()
        cls = get_modify_class(config.get_modify_methods())
        data = {"schema": cls.model_json_schema(), "sid": self._target_sid}
        self.socket.emit("modifier:schema", data)

    def selection_schema(self):
        config = GlobalConfig.load()
        cls = get_selection_class(config.get_selection_methods())

        self.socket.emit(
            "selection:schema",
            {
                "schema": cls.model_json_schema(),
                "sid": self._target_sid,
            },
        )

    def draw_schema(self):
        self.socket.emit(
            "draw:schema",
            {"schema": Geometry.updated_schema(), "sid": self._target_sid},
        )

    def modifier_run(self, data):
        with self._set_sid(data["sid"]):
            config = GlobalConfig.load()
            cls = get_modify_class(config.get_modify_methods())
            modifier = cls(**data["params"])

            if len(self) > self.step + 1:
                del self[self.step + 1 :]

            selection = self.selection
            self.selection = []

            print(f"getting {self.step} from atoms with length {len(self)}")
            atoms = self[self.step]
            print(f"Found {atoms = }")
            points = self.points
            segments = self.segments
            json_data = atoms_to_json(atoms)

            for atoms in modifier.run(
                atom_ids=selection,
                atoms=atoms,
                points=points,
                segments=segments,
                json_data=json_data,
                url=data["url"],
            ):
                self.append(atoms)
                self.step += 1

    def selection_run(self, data):
        with self._set_sid(data["sid"]):
            config = GlobalConfig.load()
            cls = get_selection_class(config.get_selection_methods())
            atoms = self[self.step]

            try:
                selection = cls(**data["params"])
                self.selection = selection.get_ids(atoms, self.selection)
            except ValueError as err:
                print(err)

    def analysis_run(self, data):
        with self._set_sid(data["sid"]):
            config = GlobalConfig.load()
            cls = get_analysis_class(config.get_analysis_methods())

            try:
                instance = cls(**data["params"])
                fig = instance.run(list(self), self.selection)
                data = {"figure": fig.to_json(), "sid": self._target_sid}
                self.socket.emit("analysis:figure", data)
            except ValueError as err:
                print(err)

    def upload_file(self, data):
        with self._set_sid(data["sid"]):
            data = data["data"]

            format = data["filename"].split(".")[-1]

            if format == "h5":
                raise ValueError("H5MD format not supported for uploading yet")
                # import znh5md
                # stream = BytesIO(data["content"].encode("utf-8"))
                # atoms = znh5md.ASEH5MD(stream).get_atoms_list()
                # for idx, atoms in tqdm.tqdm(enumerate(atoms)):
                #     atoms_dict = atoms_to_json(atoms)
                #     io.emit("atoms:upload", {idx: atoms_dict})
            else:
                stream = StringIO(data["content"])
                del self[:]
                for idx, atoms in tqdm.tqdm(
                    enumerate(ase.io.iread(stream, format=format))
                ):
                    self.append(atoms)
                    self.step = idx

    def download_file(self, data):
        with self._set_sid(data["sid"]):
            atoms_list = list(self)
            if "selection" in data:
                atoms_list = [atoms_list[data["selection"]] for atoms in atoms_list]

            file = StringIO()
            ase.io.write(file, atoms_list, format="extxyz")
            file.seek(0)
            self.socket.emit(
                "download:response", {"data": file.read(), "sid": self._target_sid}
            )


@dataclasses.dataclass
class ZnDraw(ZnDrawBase):
    """ZnDraw client.

    Attributes
    ----------
    display_new : bool
        Display new atoms in the webclient, when they are added.

    """

    url: str = None
    token: str = None

    jupyter: bool = False
    display_new: bool = True

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
                    "multiprocessing": False,
                },
            )
            self._view_thread.start()
            self.url = f"http://127.0.0.1:{port}"

        self._connect()
        #     self._view_thread = None
        #     if isinstance(self.socket, flask_socketio.SocketIO):
        #         pass
        #     else:
        #         if self.url is None:
        #             from zndraw.view import view

        #             port = get_port()
        #             self._view_thread = threading.Thread(
        #                 target=view,
        #                 kwargs={
        #                     "filename": None,
        #                     "port": port,
        #                     "open_browser": not self.jupyter,
        #                     "webview": False,
        #                     "fullscreen": False,
        #                     "start": 0,
        #                     "stop": None,
        #                     "step": 1,
        #                     "compute_bonds": True,
        #                     "multiprocessing": False,
        #                 },
        #                 daemon=True,
        #             )
        #             self._view_thread.start()
        #             self.url = f"http://127.0.0.1:{port}"

        #         self.socket.on(
        #             "connect", lambda: self.socket.emit("join", {"uuid": self.token})
        #         )

        #         self.socket.on(
        #             "join",
        #             lambda *args: self.read(
        #                 self.file.name, self.file.start, self.file.stop, self.file.step
        #             ),
        #         )
        #         self.socket.on("modifier:run", self._run_modifier)
        #         self.socket.on("selection:run", self._run_selection)
        #         self.socket.on("analysis:run", self._run_analysis)
        #         self.socket.on("download:request", self._download_file)
        #         self.socket.on("upload", self._upload_file)
        #         self.socket.on("scene:step", lambda step: setattr(self, "_step", step))
        #         self.socket.on(
        #             "scene:line",
        #             lambda data: setattr(self, "_line", process_line_data(data)),
        #         )

        #         self.socket.on("disconnect", lambda: self.disconnect())

        #         for _ in range(self._retries):
        #             with contextlib.suppress(socketio.exceptions.ConnectionError):
        #                 self.socket.connect(self.url)
        #                 break
        #             time.sleep(1)
        #         else:
        #             self.socket.connect(self.url)
        #         if not self.jupyter:
        #             self.socket.sleep(2)  # wait for the server to start

        #     if self.wait:
        #         self.socket.wait()

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
            print("Waiting for ZnDraw client to close")
            # self._view_thread.terminate()
            self._view_thread.join()
            # raise ValueError("ZnDraw client closed")

    def _repr_html_(self):
        from IPython.display import IFrame

        return IFrame(src=self.url, width="100%", height="600px")._repr_html_()

    def log(self, message: str) -> None:
        """Log a message to the console"""
        self.socket.emit("message:log", {"message": message, "sid": self.token})

    def get_logging_handler(self) -> ZnDrawLoggingHandler:
        return ZnDrawLoggingHandler(self)

    def __setitem__(self, index, value):
        super().__setitem__(index, value)
        if self.display_new:
            self.step = index

    # def _download_file(self, data):
    #     atoms = list(self)
    #     if "selection" in data:
    #         atoms = [atoms[data["selection"]] for atoms in atoms]
    #     import ase.io

    #     file = StringIO()
    #     ase.io.write(file, atoms, format="extxyz")
    #     file.seek(0)

    #     self.socket.emit("download:response", file.read())

    # def _upload_file(self, data):
    #     from io import StringIO

    #     import ase.io
    #     import tqdm

    #     # tested with small files only

    #     format = data["filename"].split(".")[-1]
    #     if format == "h5":
    #         print("H5MD format not supported for uploading yet")
    #         # import znh5md
    #         # stream = BytesIO(data["content"].encode("utf-8"))
    #         # atoms = znh5md.ASEH5MD(stream).get_atoms_list()
    #         # for idx, atoms in tqdm.tqdm(enumerate(atoms)):
    #         #     atoms_dict = atoms_to_json(atoms)
    #         #     io.emit("atoms:upload", {idx: atoms_dict})
    #     else:
    #         stream = StringIO(data["content"])
    #         del self[:]
    #         for atoms in tqdm.tqdm(ase.io.iread(stream, format=format)):
    #             self.append(atoms)
