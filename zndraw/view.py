
import logging
import multiprocessing


import webbrowser



from zndraw.app import app, io

try:
    import urllib.request

    import webview as wv
except ImportError:
    wv = None

log = logging.getLogger("werkzeug")
log.setLevel(logging.ERROR)


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
    open_browser: bool = True,
    webview: bool = True,
    fullscreen: bool = False,
    stride: int = 1,
):
    if filename is not None:
        app.config["filename"] = filename
        app.config["stride"] = stride
    url = f"http://127.0.0.1:{port}"
    print(f"Starting ZnDraw server at {url}")

    if wv is not None and webview:
        multiprocessing.Process(
            target=_view_with_webview, args=(url, fullscreen), daemon=True
        ).start()
    elif open_browser:
        webbrowser.open(url)
    io.run(app, port=port, host="0.0.0.0")


