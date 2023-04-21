import contextlib
import pathlib
import socket
import sys
import webbrowser

import typer

from zndraw import __version__, app, globals

try:
    import webview as wv
except ImportError:
    wv = None

cli = typer.Typer()


def version_callback():
    typer.echo(f"ZnDraw {__version__}")
    raise typer.Exit()


@cli.command()
def main(
    file: str = typer.Argument(..., help="Trajectory File"),
    port: int = typer.Option(None, help="Port to run the server on"),
    animate: bool = typer.Option(False, help="Animate the trajectory"),
    sphere_size: float = typer.Option(
        globals.config.sphere_size, help="size of the hydrogen sphere"
    ),
    bond_size: float = typer.Option(globals.config.bond_size, help="size of a bond"),
    max_fps: int = typer.Option(
        globals.config.max_fps, help="Maximum frames per second"
    ),
    update_function: str = typer.Option(
        None, help="Path to a python file with an update function 'module.function'."
    ),
    resolution: int = typer.Option(
        globals.config.resolution, help="Proportional to the number of polygons."
    ),
    restart_animation: bool = typer.Option(
        False, help="run the animation in an endless loop."
    ),
    frames_per_post: int = typer.Option(
        globals.config.frames_per_post, help="Number of frames to send per POST."
    ),
    browser: bool = typer.Option(True, help="Open the browser automatically."),
    webview: bool = typer.Option(True, help="Use the webview library if available."),
    verbose: bool = typer.Option(False, help="Run the server in verbose mode."),
):
    """ZnDraw: Visualize Molecules

    The ZnDraw CLI. Use 'zndraw version' to get the current version.
    """
    if not verbose:
        import logging

        log = logging.getLogger("werkzeug")
        log.setLevel(logging.ERROR)

    sys.path.insert(1, pathlib.Path.cwd().as_posix())
    if port is None:
        sock = socket.socket()
        sock.bind(("", 0))
        port = sock.getsockname()[1]
        sock.close()

    if file == "version":
        version_callback()

    if not pathlib.Path(file).exists():
        typer.echo(f"File {file} does not exist.")
        raise typer.Exit()

    globals.config.file = file
    globals.config.animate = animate
    globals.config.sphere_size = sphere_size
    globals.config.bond_size = bond_size
    globals.config.max_fps = max_fps
    globals.config.update_function = update_function
    globals.config.resolution = resolution
    globals.config.frames_per_post = frames_per_post
    globals.config.restart_animation = restart_animation

    if wv is not None and webview:
        wv.create_window("ZnDraw", app)
        with contextlib.suppress(wv.WebViewException):
            wv.start()
            return
    if browser:
        webbrowser.open(f"http://127.0.0.1:{port}")
    app.run(port=port)
