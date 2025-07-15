import dataclasses
import datetime
import os
import pathlib
import shutil
import subprocess
import typing as t
import webbrowser

import typer

from zndraw.app import create_app
from zndraw.base import FileIO
from zndraw.standalone import run_celery_worker
from zndraw.tasks import read_file, read_plots
from zndraw.upload import upload
from zndraw.utils import get_port
from zndraw_app.healthcheck import run_healthcheck

cli = typer.Typer()


@dataclasses.dataclass
class EnvOptions:
    FLASK_PORT: str | None = None
    FLASK_STORAGE: str | None = None
    FLASK_AUTH_TOKEN: str | None = None
    FLASK_TUTORIAL: str | None = None
    FLASK_SIMGEN: str | None = None
    FLASK_SERVER_URL: str | None = None
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
        envvar="ZNDRAW_URL",
    ),
    append: bool = typer.Option(
        False, help="Append the file to the existing data on the server."
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
    plots: list[str] = typer.Option(
        None, "--plots", "-p", help="List of plots to be shown in the ZnDraw GUI."
    ),
    convert_nan: bool = typer.Option(
        False,
        help="Convert NaN values to None. This is slow and experimental, but if your file contains NaN/inf values, it is required.",
        envvar="ZNDRAW_CONVERT_NAN",
    ),
    healthcheck: bool = typer.Option(False, help="Run the healthcheck."),
    follow: bool = typer.Option(
        False,
        help="Follow the file and update the visualization in real-time"
    )
):
    """Start the ZnDraw server.

    Visualize Trajectories, Structures, and more in ZnDraw.
    """
    if healthcheck:
        if url is None:
            raise ValueError("You need to provide a URL to use the healthcheck feature.")
        run_healthcheck(url)
    if plots is None:
        plots = []
    if token is not None and url is None:
        raise ValueError("You need to provide a URL to use the token feature.")
    if url is not None and port is not None:
        raise ValueError(
            "You cannot provide a URL and a port at the same time. Use something like '--url http://localhost:1234' instead."
        )
    if follow and url is None:
        raise ValueError("You need to provide a URL to use the follow feature.")

    env_config = EnvOptions.from_env()

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
        env_config.FLASK_STORAGE = f"znsocket://localhost:{env_config.FLASK_PORT}"

    env_config.save_to_env()

    if remote is None and rev is None and filename is not None:
        if not pathlib.Path(filename).exists():
            typer.echo(f"File {filename} does not exist.")
            raise typer.Exit(code=1)

    worker = None  # Initialize worker variable

    fileio = FileIO(
        name=filename,
        remote=remote,
        rev=rev,
        start=start,
        stop=stop,
        step=step,
        convert_nan=convert_nan,
    )

    if url is not None:
        upload(url, token, fileio, append, plots, browser, 16, follow)
        return

    typer.echo(
        f"{datetime.datetime.now().isoformat()}: Starting zndraw server on port {port}"
    )

    app = create_app(main=True)

    if browser:
        webbrowser.open(f"http://localhost:{env_config.FLASK_PORT}")

    socketio = app.extensions["socketio"]

    # Start celery worker as separate process
    if standalone and url is None:
        worker = run_celery_worker()
        # Start the initial tasks after worker is ready
        read_file.s(fileio.to_dict()).apply_async()
        read_plots.s(plots, fileio.remote, fileio.rev).apply_async()
    else:
        # If not standalone, start tasks immediately
        read_file.s(fileio.to_dict()).apply_async()
        read_plots.s(plots, fileio.remote, fileio.rev).apply_async()

    try:
        socketio.run(
            app,
            host="0.0.0.0",
            port=app.config["PORT"],
        )
    finally:
        try:
            print(
                f"{datetime.datetime.now().isoformat()}: Shutting down zndraw server {env_config.FLASK_SERVER_URL}"
            )

            if standalone and url is None and worker is not None:
                print("---------------------- SHUTDOWN CELERY ----------------------")
                try:
                    # Try graceful shutdown first
                    celery_app = app.extensions["celery"]
                    celery_app.control.broadcast("shutdown")
                    # TODO: check that all workers are actually shutdown
                    print("Waiting for celery worker to shutdown...")
                    worker.wait(timeout=5)
                    print("Celery worker shutdown successfully")
                except subprocess.TimeoutExpired:
                    print("Celery worker didn't shutdown gracefully, terminating...")
                    worker.terminate()
                    worker.wait()
                    print("Celery worker terminated")
                except Exception as e:
                    print(f"Error during celery shutdown: {e}")
                    # Force terminate if anything goes wrong
                    try:
                        worker.terminate()
                        worker.wait()
                    except:
                        pass
        finally:
            if app.config["CELERY"]["broker_url"] == "filesystem://":
                print("---------------------- REMOVE CELERY CTRL ----------------------")
                for path in app.config["CELERY"]["broker_transport_options"].values():
                    if os.path.exists(path):
                        shutil.rmtree(path)

            print("---------------------- SHUTDOWN COMPLETE ----------------------")
