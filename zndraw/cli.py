from typing import Optional

import typer

from zndraw.view import _get_port, view

cli = typer.Typer()


@cli.command()
def main(
    filename: Optional[str] = typer.Argument(
        None, help="Path to the file which should be visualized in ZnDraw."
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
    stride: int = typer.Option(
        1,
        help="""Stride for the frames to be visualized. If set to 1, all frames will be visualized.
        If e.g. set to 2, every second frame will be visualized.""",
    ),
):
    """Start the ZnDraw server.

    Visualize Trajectories, Structures, and more in ZnDraw.
    """
    if port is None:
        port = _get_port()
    view(
        filename,
        port,
        webview=webview,
        fullscreen=fullscreen,
        open_browser=browser,
        stride=stride,
    )
