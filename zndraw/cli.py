from typing import Optional

import typer

from zndraw.utils import get_port
from zndraw.view import view

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
    compute_bonds: bool = typer.Option(
        True,
        help="""Whether to compute bonds for the structure. If set to False, no bonds will be computed.""",
    ),
    upgrade_insecure_requests: bool = typer.Option(
        False,
        hidden=True,
        help="Set the html attribute upgrade-insecure-requests. If you are running ZnDraw behind a reverse proxy and encounter issues with insecure requests, you might want to set this to true.",
    ),
    use_token: bool = typer.Option(
        False,
        help="Use a token to authenticate the ZnDraw server.  This is useful if you are running ZnDraw as a server application.",
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
):
    """Start the ZnDraw server.

    Visualize Trajectories, Structures, and more in ZnDraw.
    """
    if port is None:
        port = get_port()

    view(
        filename,
        port,
        webview=webview,
        fullscreen=fullscreen,
        open_browser=browser,
        start=start,
        stop=stop,
        step=step,
        compute_bonds=compute_bonds,
        upgrade_insecure_requests=upgrade_insecure_requests,
        use_token=use_token,
        remote=remote,
        rev=rev,
        tutorial=tutorial,
        auth_token=auth_token,
    )
