import eventlet

eventlet.monkey_patch()

import dataclasses
import datetime
import os
import typing as t

import typer

from zndraw.app import create_app
from zndraw.base import FileIO
from zndraw.standalone import run_celery_worker, run_znsocket
from zndraw.tasks import read_file
from zndraw.upload import upload
from zndraw.utils import get_port

cli = typer.Typer()


@dataclasses.dataclass
class EnvOptions:
    FLASK_PORT: str | None = None
    FLASK_STORAGE: str | None = None
    FLASK_AUTH_TOKEN: str | None = None
    FLASK_TUTORIAL: str | None = None
    FLASK_SIMGEN: str | None = None
    FLASK_SERVER_URL: str | None = None
    FLASK_STORAGE_PORT: str | None = None
    FLASK_COMPUTE_BONDS: str | None = None
    FLASK_MAX_HTTP_BUFFER_SIZE: str | None = None

    @classmethod
    def from_env(cls):
        return cls(
            **{
                field.name: os.environ.get(field.name)
                for field in dataclasses.fields(cls)
            }
        )

    def save_to_env(self):
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if value is not None:
                os.environ[field.name] = value


@cli.command()
def main(
    filename: t.Optional[str] = typer.Argument(
        None,
        help="Path to the file which should be visualized in ZnDraw. Can also be the name and attribute of a ZnTrack Node like 'MyNode.atoms' if at least '--remote .' is provided. ",
    ),
    url: t.Optional[str] = typer.Option(
        None,
        help="URL to a running ZnDraw server. Use this server instead of starting a new one.",
    ),
    token: t.Optional[str] = typer.Option(
        None, help="Only valid if 'url' is provided. Room token to upload the file to."
    ),
    port: int = typer.Option(
        None, help="""Port to use for the ZnDraw server. Default port is 1234"""
    ),
    browser: bool = typer.Option(
        True, help="""Whether to open the ZnDraw GUI in the default web browser."""
    ),
    start: int = typer.Option(
        None,
        help="""First frame to be visualized. If set to 0, the first frame will be visualized.""",
    ),
    stop: int = typer.Option(
        None,
        help="""Last frame to be visualized. If set to None, the last frame will be visualized.""",
    ),
    step: int = typer.Option(
        None,
        help="""Stepsize for the frames to be visualized. If set to 1, all frames will be visualized.
        If e.g. set to 2, every second frame will be visualized.""",
    ),
    remote: str = typer.Option(
        None,
        help="URL to a ZnTrack repository to stream data from.",
    ),
    rev: str = typer.Option(
        None,
        help="Revision of the ZnTrack repository to stream data from.",
    ),
    tutorial: str = typer.Option(
        None,
        help="Show the tutorial from the URL inside an IFrame.",
    ),
    auth_token: str = typer.Option(
        None,
        help="Token to authenticate pyclient requests to the ZnDraw server, e.g., for adding defaults to all webclients.",
    ),
    simgen: bool = typer.Option(
        False,
        help="Show the SiMGen demo UI.",
    ),
    storage: str = typer.Option(
        None,
        help="URL to the redis `redis://localhost:6379/0` or znsocket `znsocket://127.0.0.1:6379` server. If None is provided, a local znsocket server will be started.",
    ),
    storage_port: int = typer.Option(
        None, help="Port to use for the storage server. Default port is 6374"
    ),
    standalone: bool = typer.Option(
        True,
        help="Run ZnDraw without additional tools. If disabled, redis and celery must be started manually.",
    ),
    bonds: bool = typer.Option(
        True,
        help="Compute bonds based on covalent distances. This can be slow for large structures.",
    ),
    max_http_buffer_size: int = typer.Option(
        None, help="Maximum size of the HTTP buffer in bytes. Default is 1MB."
    ),
):
    """Start the ZnDraw server.

    Visualize Trajectories, Structures, and more in ZnDraw.
    """
    if token is not None and url is None:
        raise ValueError("You need to provide a URL to use the token feature.")
    if url is not None and port is not None:
        raise ValueError(
            "You cannot provide a URL and a port at the same time. Use something like '--url http://localhost:1234' instead."
        )

    env_config = EnvOptions.from_env()

    if storage_port is not None:
        env_config.FLASK_STORAGE_PORT = str(storage_port)
    elif env_config.FLASK_STORAGE_PORT is None:
        env_config.FLASK_STORAGE_PORT = str(get_port(default=6374))

    if port is not None:
        env_config.FLASK_PORT = str(port)
    elif env_config.FLASK_PORT is None:
        env_config.FLASK_PORT = str(get_port(default=1234))
    if storage is not None:
        env_config.FLASK_STORAGE = storage
    if auth_token is not None:
        env_config.FLASK_AUTH_TOKEN = auth_token
    if tutorial is not None:
        env_config.FLASK_TUTORIAL = tutorial
    if simgen:
        env_config.FLASK_SIMGEN = "TRUE"
    if bonds:
        env_config.FLASK_COMPUTE_BONDS = "TRUE"
    if max_http_buffer_size is not None:
        env_config.FLASK_MAX_HTTP_BUFFER_SIZE = str(int(max_http_buffer_size))

    env_config.FLASK_SERVER_URL = f"http://localhost:{env_config.FLASK_PORT}"

    if standalone and storage is None:
        env_config.FLASK_STORAGE = f"znsocket://localhost:{env_config.FLASK_STORAGE_PORT}"

    env_config.save_to_env()

    if standalone and url is None:
        if env_config.FLASK_STORAGE.startswith("znsocket"):
            # standalone with redis would assume a running instance of redis
            server = run_znsocket(env_config.FLASK_STORAGE_PORT)
        worker = run_celery_worker()

    fileio = FileIO(
        name=filename,
        remote=remote,
        rev=rev,
        start=start,
        stop=stop,
        step=step,
    )

    if url is not None:
        upload(filename, url, token, fileio)
        return

    typer.echo(
        f"{datetime.datetime.now().isoformat()}: Starting zndraw server on port {port}"
    )

    app = create_app()

    read_file.delay(fileio.to_dict())

    if browser:
        import webbrowser

        webbrowser.open(f"http://localhost:{env_config.FLASK_PORT}")

    socketio = app.extensions["socketio"]
    try:
        socketio.run(
            app,
            host="0.0.0.0",
            port=app.config["PORT"],
        )
    finally:
        if standalone:
            server.terminate()
            server.wait()
            print("znsocket server terminated.")
            worker.terminate()
            worker.wait()
            print("celery worker terminated.")
