import logging
import socket
import webbrowser

import typer

from zndraw.app import app, socketio

log = logging.getLogger("werkzeug")
log.setLevel(logging.ERROR)

cli = typer.Typer()


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


@cli.command()
def main(filename: str):
    # get an empty port
    port = _get_port()
    app.config["filename"] = filename
    url = f"http://127.0.0.1:{port}"

    webbrowser.open(url)
    socketio.run(app, port=port, host="0.0.0.0")
