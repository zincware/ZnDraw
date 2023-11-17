import logging
import multiprocessing as mp
import webbrowser

from zndraw.app import app, io
from zndraw.zndraw import FileIO, ZnDrawDefault

try:
    import urllib.request

    import webview as wv
except ImportError:
    wv = None

# werkzeug_log = logging.getLogger("werkzeug")
# werkzeug_log.setLevel(logging.ERROR)

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def _view_with_webview(url, fullscreen=False):
    wv.create_window(
        "ZnDraw",
        url,
        width=1280,
        height=720,
        resizable=True,
        fullscreen=fullscreen,
        confirm_close=True,
        background_color="#FFFFFF",
    )
    wv.start()
    try:
        urllib.request.urlopen(url + "/exit")
    except urllib.error.RemoteDisconnected:
        pass


def view(
    filename: str,
    port: int,
    open_browser: bool,
    webview: bool,
    fullscreen: bool,
    start: int,
    stop: int,
    step: int,
    compute_bonds: bool,
    use_token: bool = False,
    upgrade_insecure_requests: bool = False,
    remote: str = None,
    rev: str = None,
):
    if not use_token:
        app.config["token"] = "notoken"
    app.config["upgrade_insecure_requests"] = upgrade_insecure_requests
    app.config["compute_bonds"] = compute_bonds
    url = f"http://127.0.0.1:{port}"

    file_io = FileIO(filename, start, stop, step, remote, rev)

    proc = mp.Process(
        target=ZnDrawDefault,
        kwargs={"url": url, "token": "default", "file_io": file_io},
    )
    proc.start()

    log.critical(f"Starting ZnDraw server at {url}")

    if wv is not None and webview:
        wv_proc = mp.Process(
            target=_view_with_webview, args=(url, fullscreen), daemon=True
        )
        wv_proc.start()
    elif open_browser:
        webbrowser.open(url)

    io.run(app, port=port, host="0.0.0.0")

    proc.terminate()
    proc.join()
    if wv is not None and webview:
        wv_proc.terminate()
        wv_proc.join()
