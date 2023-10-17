import collections.abc
import contextlib
import dataclasses
import pathlib
import threading
import time
import typing as t
from io import StringIO
import importlib

import ase
import ase.io
import flask_socketio
import networkx as nx
import numpy as np
import socketio
import tqdm
import znh5md

from zndraw.bonds import ASEComputeBonds
from zndraw.data import atoms_from_json, atoms_to_json
from zndraw.select import get_selection_class
from zndraw.utils import ZnDrawLoggingHandler, get_port
from zndraw.settings import GlobalConfig
from zndraw.draw import Geometry


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
        while True:
            try:
                self.socket.connect(self.url)
                break
            except socketio.exceptions.ConnectionError:
                time.sleep(0.1)
        self.socket.wait()

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
            if isinstance(index, slice):
                index = range(*index.indices(length))

            index = [index] if isinstance(index, int) else index
            index = [i if i >= 0 else length + i for i in index]

            data = {"index": index}
            if self._target_sid is not None:
                data["sid"] = self._target_sid

            self.socket.emit("atoms:delete", data)
            if index[0] >= length or index[-1] >= length:
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
        if isinstance(index, slice):
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
        if data == []:
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
            return int(self.socket.call("scene:step", {"sid": self._target_sid}))
        return int(self.socket.call("scene:step", {}))

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
        self._connect()

    def initialize_webclient(self, sid):
        with self._set_sid(sid):
            self.read_data()
            self.analyis_schema()
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

    def analyis_schema(self):
        config = GlobalConfig.load()
        atoms = self[0]

        for modifier in config.analysis_functions:
            module_name, function_name = modifier.rsplit(".", 1)
            module = importlib.import_module(module_name)
            modifier_cls = getattr(module, function_name)
            schema = modifier_cls.schema_from_atoms(atoms)

            self.socket.emit(
                "analysis:schema",
                {"name": modifier, "schema": schema, "sid": self._target_sid},
            )

    def modifier_schema(self):
        config = GlobalConfig.load()

        for modifier in config.modify_functions:
            module_name, function_name = modifier.rsplit(".", 1)
            try:
                module = importlib.import_module(module_name)
                modifier_cls = getattr(module, function_name)
                schema = modifier_cls.model_json_schema()

                if modifier in config.function_schema:
                    kwargs = config.function_schema[modifier]
                    for key, value in kwargs.items():
                        schema["properties"][key]["default"] = value

                self.socket.emit(
                    "modifier:schema",
                    {"name": modifier, "schema": schema, "sid": self._target_sid},
                )
            except ImportError:
                print(f"Could not import {modifier}")

    def selection_schema(self):
        self.socket.emit(
            "selection:schema",
            {
                "schema": get_selection_class().model_json_schema(),
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
            module_name, function_name = data["name"].rsplit(".", 1)
            module = importlib.import_module(module_name)
            modifier_cls = getattr(module, function_name)
            modifier = modifier_cls(**data["params"])

            if len(self) > self.step + 1:
                del self[self.step + 1 :]

            selection = self.selection
            self.selection = []

            for idx, atoms in enumerate(
                modifier.run(
                    atom_ids=selection,
                    atoms=self[self.step],
                    points=self.points,
                    segments=self.segments,
                    json_data=data["atoms"] if "atoms" in data else None,
                    url=data["url"],
                )
            ):
                self[self.step + idx + 1] = atoms
                self.step += 1

    def selection_run(self, data):
        with self._set_sid(data["sid"]):
            atoms = self[self.step]
            selection = get_selection_class()(**data["params"])
            self.selection = selection.get_ids(atoms, self.selection)

    def analysis_run(self, data):
        with self._set_sid(data["sid"]):
            atoms_list = list(self)

            module_name, function_name = data["name"].rsplit(".", 1)
            module = importlib.import_module(module_name)
            cls = getattr(module, function_name)
            instance = cls(**data["params"])

            fig = instance.run(atoms_list, self.selection)

            data = {"figure": fig.to_json(), "sid": self._target_sid}

            self.socket.emit("analysis:figure", data)


@dataclasses.dataclass
class ZnDraw(ZnDrawBase):
    jupyter: bool = False

    def __post_init__(self):
        super().__post_init__()
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

    def _repr_html_(self):
        from IPython.display import IFrame

        return IFrame(src=self.url, width="100%", height="600px")._repr_html_()

    def log(self, message: str) -> None:
        """Log a message to the console"""
        self.socket.emit("message:log", message)

    def get_logging_handler(self) -> ZnDrawLoggingHandler:
        return ZnDrawLoggingHandler(self.socket)

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
