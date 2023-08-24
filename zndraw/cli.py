import socket
import webbrowser

import typer
from dask.distributed import Client, Variable

cli = typer.Typer()

from zndraw.app import app, socketio
from zndraw.data import DataHandler


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

    import logging

    log = logging.getLogger("werkzeug")
    log.setLevel(logging.ERROR)

    with Client() as client:
        _filename = Variable("filename")
        _filename.set(filename)

        url = Variable("url")
        url.set(f"http://127.0.0.1:{port}")

        DataHandler(client).create_dataset(filename)

        app.config["dask-scheduler"] = client.scheduler_info()["address"]

        print(f"Starting server on {url.get()} with scheduler {app.config['dask-scheduler']}")

        webbrowser.open(url.get())
        socketio.run(app, port=port, host="0.0.0.0")
