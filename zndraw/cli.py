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
    )
):
    """Start the ZnDraw server.

    Visualize Trajectories, Structures, and more in ZnDraw.
    """
    view(filename, _get_port(), webview=webview, fullscreen=fullscreen)
