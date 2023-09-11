import collections.abc
import contextlib
import dataclasses
import pathlib
import threading
import time
import typing as t

import ase
import ase.io
import flask_socketio
import networkx as nx
import socketio
import tqdm
import znh5md

from zndraw.bonds import ASEComputeBonds
from zndraw.data import atoms_from_json, atoms_to_json
from zndraw.utils import get_port


@dataclasses.dataclass
class ZnDraw(collections.abc.MutableSequence):
    url: str = None
    socket: socketio.Client = dataclasses.field(default_factory=socketio.Client)
    jupyter: bool = False
    bonds_calculator: ASEComputeBonds = dataclasses.field(
        default_factory=ASEComputeBonds
    )

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
                "connect", lambda: print(f"Connected to ZnDraw server at {self.url}")
            )

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
        get_size_event = threading.Event()
        self.socket.emit("atoms:size", {})

        _size = None

        def on_size(size):
            nonlocal _size
            _size = size
            get_size_event.set()

        self.socket.on("atoms:size", on_size)
        get_size_event.wait()
        return _size

    def _set_item(self, index, value):
        assert isinstance(value, ase.Atoms), "Must be an ASE Atoms object"
        assert isinstance(index, int), "Index must be an integer"
        if self.bonds_calculator is not None:
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
        raise NotImplementedError

    def append(self, value: ase.Atoms) -> None:
        """Append atoms to the end of the list"""
        self[len(self)] = value

    def extend(self, values: list[ase.Atoms]) -> None:
        for val in values:
            self._set_item(len(self), val)
        if self.display_new:
            self.display(len(self) - 1)

    def read(self, filename: str, start: int, stop: int, step: int):
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
