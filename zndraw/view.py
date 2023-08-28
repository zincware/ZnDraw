import collections.abc
import dataclasses
import logging
import socket
import threading
import typing as t
import webbrowser

import ase
import socketio

from zndraw.app import app, io
from zndraw.data import atoms_from_json, atoms_to_json

log = logging.getLogger("werkzeug")
log.setLevel(logging.ERROR)


def _get_port() -> int:
    try:
        sock = socket.socket()
        sock.bind(("", 1234))
        port = 1234
    except OSError:
        sock = socket.socket()
        sock.bind(("", 0))
        port = sock.getsockname()[1]
    finally:
        sock.close()
    return port


def view(filename: str, port: int, open_browser: bool = True):
    if filename is not None:
        app.config["filename"] = filename
    url = f"http://127.0.0.1:{port}"

    if open_browser:
        webbrowser.open(url)
    io.run(app, port=port, host="0.0.0.0")


@dataclasses.dataclass
class ZnDraw(collections.abc.MutableSequence):
    url: str = None
    socket: socketio.Client = dataclasses.field(default_factory=socketio.Client)
    jupyter: bool = False

    def __post_init__(self):
        self._view_thread = None
        if self.url is None:
            port = _get_port()
            self._view_thread = threading.Thread(
                target=view, args=(None, port, not self.jupyter), daemon=True
            )
            self._view_thread.start()
            self.url = f"http://127.0.0.1:{port}"
        self.socket.connect(self.url)
        if not self.jupyter:
            self.socket.sleep(2)  # wait for the server to start

    # def __del__(self):
    #     if self._view_thread is not None:
    #         # just open the url and don't expect a response
    #         self.socket.emit("exit")
    #         self.socket.disconnect()
    #         self._view_thread.join()

    def view(self, atoms_list):
        if isinstance(atoms_list, ase.Atoms):
            atoms_list = [atoms_list]
        for idx, atoms in enumerate(atoms_list):
            self._set_item(idx, atoms)
        self.display(idx)

    def _repr_html_(self):
        from IPython.display import IFrame

        return IFrame(src=self.url, width="100%", height="600px")._repr_html_()

    def __delitem__(self, index):
        pass

    def __getitem__(self, index) -> t.Union[ase.Atoms, list[ase.Atoms]]:
        get_item_event = threading.Event()

        if not isinstance(index, int) and not isinstance(index, list):
            raise TypeError("Index must be an integer or list of integers")

        index = [index] if isinstance(index, int) else index

        self.socket.emit("atoms:download", index)

        downloaded_data = []

        def on_download(data):
            nonlocal downloaded_data
            for key, val in data.items():
                downloaded_data.append(atoms_from_json(val))
            get_item_event.set()

        self.socket.on("atoms:download", on_download)
        get_item_event.wait()

        return downloaded_data[0] if len(downloaded_data) == 1 else downloaded_data

    def __len__(self):
        pass

    def _set_item(self, index, value):
        assert isinstance(value, ase.Atoms), "Must be an ASE Atoms object"
        assert isinstance(index, int), "Index must be an integer"
        self.socket.emit("atoms:upload", {index: atoms_to_json(value)})
        self.socket.emit("view:set", index)

    def __setitem__(self, index, value):
        self._set_item(index, value)
        self.display(index)

    def display(self, index):
        """Display the atoms at the given index"""
        self.socket.emit("view:set", index)

    def insert(self, index, value):
        pass
