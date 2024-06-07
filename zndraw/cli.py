import os
from typing import Optional

import typer

from zndraw.app import FileIO, ZnDrawServer
from zndraw.utils import get_port

cli = typer.Typer()


@cli.command()
def main(
    filename: Optional[str] = typer.Argument(
        None,
        help="Path to the file which should be visualized in ZnDraw. Can also be the name and attribute of a ZnTrack Node like 'MyNode.atoms' if at least '--remote .' is provided. ",
    ),
    webview: bool = typer.Option(
        True,
        help="""Whether to use the webview library to open the ZnDraw GUI.
        If `pip install pywebview` is available, webview will be used.
        Otherwise, the GUI will be opened in the default web browser.""",
    ),
    fullscreen: bool = typer.Option(
        False,
        help="Use fullscreen mode for the ZnDraw GUI. (only with webview)",
    ),
    port: int = typer.Option(
        None, help="""Port to use for the ZnDraw server. Default port is 1234"""
    ),
    browser: bool = typer.Option(
        True, help="""Whether to open the ZnDraw GUI in the default web browser."""
    ),
    start: int = typer.Option(
        0,
        help="""First frame to be visualized. If set to 0, the first frame will be visualized.""",
    ),
    stop: int = typer.Option(
        None,
        help="""Last frame to be visualized. If set to None, the last frame will be visualized.""",
    ),
    step: int = typer.Option(
        1,
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
    celery_worker: bool = typer.Option(
        True,
        help="Start ZnDraw with celery workers. If disabled, you must manage the workers yourself. This can be useful when hosting ZnDraw for multiple users.",
    ),
    storage: str = typer.Option(
        None,
        help="URL to the redis `redis://localhost:6379/0` or znsocket `znsocket://127.0.0.1:6379` server. If None is provided, a local znsocket server will be started.",
    ),
):
    """Start the ZnDraw server.

    Visualize Trajectories, Structures, and more in ZnDraw.
    """
    if port is None:
        port = get_port()

    fileio = FileIO(
        name=filename,
        remote=remote,
        rev=rev,
        start=start,
        stop=stop,
        step=step,
    )
    if "ZNDRAW_STORAGE" in os.environ and storage is None:
        print(
            f"Using storage from environment variable ZNDRAW_STORAGE: {os.environ['ZNDRAW_STORAGE']}"
        )
        storage = os.environ["ZNDRAW_STORAGE"]
    elif storage is not None:
        os.environ["ZNDRAW_STORAGE"] = storage

    with ZnDrawServer(
        tutorial=tutorial,
        auth_token=auth_token,
        port=port,
        fileio=fileio,
        simgen=simgen,
        celery_worker=celery_worker,
        storage=storage,
    ) as app:
        app.run(browser=browser)
