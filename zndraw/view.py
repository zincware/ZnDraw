import contextlib
import logging
import socket
import typing
import webbrowser

import ase

from zndraw import app, shared

try:
    import webview as wv
except ImportError:
    wv = None


def view(
    data: list[ase.Atoms],
    verbose: bool = False,
    port: int = None,
    webview: bool = True,
    browser: bool = True,
) -> None:
    """Visualize a list of ASE atoms objects."""
    if not verbose:
        import logging

        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

    if port is None:
        sock = socket.socket()
        sock.bind(("", 0))
        port = sock.getsockname()[1]
        sock.close()

    shared.config = shared.Config(file="")
    shared.config._atoms_cache = dict(enumerate(data))

    print(shared.config)

    if wv is not None and webview:
        wv.create_window("ZnDraw", app)
        with contextlib.suppress(wv.WebViewException):
            wv.start()
            return
    if browser:
        webbrowser.open(f"http://127.0.0.1:{port}")
    app.run(port=port, host="0.0.0.0")
