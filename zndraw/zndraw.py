import collections.abc
import contextlib
import dataclasses
import pathlib
import threading
import time
import typing as t
from io import StringIO

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


@dataclasses.dataclass
class FileIO:
    name: str
    start: int = 0
    stop: int = None
    step: int = 1


@dataclasses.dataclass
class ZnDraw(collections.abc.MutableSequence):
    url: str = None
    socket: socketio.Client = dataclasses.field(default_factory=socketio.Client)
    jupyter: bool = False
    bonds_calculator: ASEComputeBonds = dataclasses.field(
        default_factory=ASEComputeBonds
    )
    file: FileIO = None
    wait: bool = False
    token: str = None

    display_new: bool = True
    _retries: int = 5

    def __post_init__(self):
        self._view_thread = None
        if isinstance(self.socket, flask_socketio.SocketIO):
            pass
        else:
            if self.url is None:
                from zndraw.view import view

                port = get_port()
                self._view_thread = threading.Thread(
                    target=view,
                    kwargs={
                        "filename": None,
                        "port": port,
                        "open_browser": not self.jupyter,
                        "webview": False,
                        "fullscreen": False,
                        "start": 0,
                        "stop": None,
                        "step": 1,
                        "compute_bonds": True,
                        "multiprocessing": False,
                    },
                    daemon=True,
                )
                self._view_thread.start()
                self.url = f"http://127.0.0.1:{port}"

            self.socket.on(
                "connect", lambda: self.socket.emit("join", {"uuid": self.token})
            )

            self.socket.on(
                "atoms:request",
                lambda url: self.read(
                    self.file.name, self.file.start, self.file.stop, self.file.step
                ),
            )
            self.socket.on("modifier:run", self._run_modifier)
            self.socket.on("selection:run", self._run_selection)
            self.socket.on("analysis:run", self._run_analysis)
            self.socket.on("download:request", self._download_file)
            self.socket.on("upload", self._upload_file)

            self.socket.on("disconnect", lambda: self.disconnect())

            for _ in range(self._retries):
                with contextlib.suppress(socketio.exceptions.ConnectionError):
                    self.socket.connect(self.url)
                    break
                time.sleep(1)
            else:
                self.socket.connect(self.url)
            if not self.jupyter:
                self.socket.sleep(2)  # wait for the server to start
        
        if self.wait:
            self.socket.wait()

    def view(self, atoms_list):
        if isinstance(atoms_list, ase.Atoms):
            atoms_list = [atoms_list]
        for idx, atoms in enumerate(atoms_list):
            self._set_item(idx, atoms)
        if self.display_new:
            self.display(idx)

    def disconnect(self):
        self.socket.disconnect()

    def _repr_html_(self):
        from IPython.display import IFrame

        return IFrame(src=self.url, width="100%", height="600px")._repr_html_()

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
            self.socket.emit("atoms:delete", index)
            if index[0] >= length or index[-1] >= length:
                raise IndexError("Index out of range")
        else:
            raise TypeError("Index must be an integer, slice or list[int]")

    def __getitem__(self, index) -> t.Union[ase.Atoms, list[ase.Atoms]]:
        get_item_event = threading.Event()

        length = len(self)
        if isinstance(index, slice):
            index = range(*index.indices(length))

        index = [index] if isinstance(index, int) else index
        index = [i if i >= 0 else length + i for i in index]

        self.socket.emit("atoms:download", index)

        downloaded_data = []

        def on_download(data):
            nonlocal downloaded_data
            for key, val in data.items():
                downloaded_data.append(atoms_from_json(val))
            get_item_event.set()

        self.socket.on("atoms:download", on_download)
        get_item_event.wait()

        data = downloaded_data[0] if len(downloaded_data) == 1 else downloaded_data
        if data == []:
            raise IndexError("Index out of range")
        return data

    def __len__(self):
        return _await_answer(self.socket, "atoms:size", data={})

    def _set_item(self, index, value):
        assert isinstance(value, ase.Atoms), "Must be an ASE Atoms object"
        assert isinstance(index, int), "Index must be an integer"
        if hasattr(value, "connectivity"):
            pass
        elif self.bonds_calculator is not None:
            value.connectivity = self.bonds_calculator.build_graph(value)
        else:
            value.connectivity = nx.Graph()

        self.socket.emit("atoms:upload", {index: atoms_to_json(value)})

    def __setitem__(self, index, value):
        self._set_item(index, value)
        if self.display_new:
            self.display(index)

    def display(self, index):
        """Display the atoms at the given index"""
        self.socket.emit("view:set", index)

    def insert(self, index, value):
        """Insert atoms before index"""
        self.socket.emit("atoms:insert", {index: atoms_to_json(value)})
        if self.display_new:
            self.display(index)

    def append(self, value: ase.Atoms) -> None:
        """Append atoms to the end of the list"""
        self[len(self)] = value

    def extend(self, values: list[ase.Atoms]) -> None:
        for val in values:
            self._set_item(len(self), val)
        if self.display_new:
            self.display(len(self) - 1)

    def log(self, message: str) -> None:
        """Log a message to the console"""
        self.socket.emit("message:log", message)

    def get_logging_handler(self) -> ZnDrawLoggingHandler:
        return ZnDrawLoggingHandler(self.socket)

    def read(self, filename: str, start: int = 0, stop: int = None, step: int = 1):
        """Read atoms from file and return a list of atoms dicts.

        Parameters
        ----------
        filename : str
            Path to the file which should be read.
        start : int
            First frame to be read. If set to 0, the first frame will be read.
        stop : int
            Last frame to be read. If set to None, the last frame will be read.
        step : int
            Stepsize for the frames to be visualized. If set to 1, all frames will be visualized.
        """

        if filename is None:
            return

        if pathlib.Path(filename).suffix == ".h5":
            # Read file using znh5md and convert to list[ase.Atoms]
            atoms_list = znh5md.ASEH5MD(filename)[start:stop:step]

        else:
            # Read file using ASE and convert to list[ase.Atoms]
            # TODO use read generator in loop
            atoms_list = list(ase.io.iread(filename))[start:stop:step]
        for idx, atoms in tqdm.tqdm(
            enumerate(atoms_list), ncols=100, total=len(atoms_list)
        ):  
            self[idx] = atoms

    def get_selection(self) -> list[int]:
        """Get the selected atoms"""
        return _await_answer(self.socket, "selection:get", data={})

    def get_line(self) -> tuple[np.ndarray, np.ndarray]:
        """Get the points of the selected atoms"""
        data = _await_answer(self.socket, "draw:get_line", data={})
        points = np.array([[val["x"], val["y"], val["z"]] for val in data["points"]])
        segments = np.array(data["segments"])

        return points, segments

    def set_selection(self, selection: list[int]) -> None:
        """Set the selected atoms"""
        self.socket.emit("selection:set", selection)

    def _run_modifier(self, data):
        import importlib

        import ase

        points = np.array([[val["x"], val["y"], val["z"]] for val in data["points"]])
        segments = np.array(data["segments"])

        if "atoms" in data:
            atoms = atoms_from_json(data["atoms"])
        else:
            atoms = ase.Atoms()

        module_name, function_name = data["name"].rsplit(".", 1)
        module = importlib.import_module(module_name)
        modifier_cls = getattr(module, function_name)
        modifier = modifier_cls(**data["params"])
        # available_methods = {x.__name__: x for x in [Explode, Duplicate]}

        # modifier = available_methods[data["name"]](**data["params"])
        print(f"modifier:run {modifier = }")
        atoms_list = modifier.run(
            atom_ids=data["selection"],
            atoms=atoms,
            points=points,
            segments=segments,
            json_data=data["atoms"] if "atoms" in data else None,
            url=data["url"],
        )
        del self[data["step"] :]
        self.extend(atoms_list)

    def _run_selection(self, data):
        import ase

        if "atoms" in data:
            atoms = atoms_from_json(data["atoms"])
        else:
            atoms = ase.Atoms()

        try:
            selection = get_selection_class()(**data["params"])
            selected_ids = selection.get_ids(atoms, data["selection"])
            self.set_selection(selected_ids)
        except ValueError as err:
            print(err)

    def _run_analysis(self, data):
        import importlib

        atoms_list = list(self)

        print(f"Analysing {len(atoms_list)} frames")

        module_name, function_name = data["name"].rsplit(".", 1)
        module = importlib.import_module(module_name)
        cls = getattr(module, function_name)
        instance = cls(**data["params"])

        fig = instance.run(atoms_list, data["selection"])
        self.socket.emit("analysis:figure", fig.to_json())

    def _download_file(self, data):
        atoms = list(self)
        if "selection" in data:
            atoms = [atoms[data["selection"]] for atoms in atoms]
        import ase.io

        file = StringIO()
        ase.io.write(file, atoms, format="extxyz")
        file.seek(0)

        self.socket.emit("download:response", file.read())

    def _upload_file(self, data):
        from io import StringIO

        import ase.io
        import tqdm

        # tested with small files only

        format = data["filename"].split(".")[-1]
        if format == "h5":
            print("H5MD format not supported for uploading yet")
            # import znh5md
            # stream = BytesIO(data["content"].encode("utf-8"))
            # atoms = znh5md.ASEH5MD(stream).get_atoms_list()
            # for idx, atoms in tqdm.tqdm(enumerate(atoms)):
            #     atoms_dict = atoms_to_json(atoms)
            #     io.emit("atoms:upload", {idx: atoms_dict})
        else:
            stream = StringIO(data["content"])
            del self[:]
            for atoms in tqdm.tqdm(ase.io.iread(stream, format=format)):
                self.append(atoms)
