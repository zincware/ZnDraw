import socket
import webbrowser

import typer
import znh5md

cli = typer.Typer()

from zndraw.app import app, socketio


def _get_port() -> int:
    sock = socket.socket()
    sock.bind(("", 0))
    port = sock.getsockname()[1]
    sock.close()
    return port


@cli.command()
def main(filename: str):
    # get an empty port
    port = _get_port()

    import logging

    log = logging.getLogger("werkzeug")
    log.setLevel(logging.ERROR)

    filename = filename.replace("/", ".slash.").replace("\\", ".slash.")
    print(f"Reading {filename}...")

    webbrowser.open(f"http://127.0.0.1:{port}/read/{filename}")
    socketio.run(app, port=port, host="0.0.0.0")
