import contextlib
import pathlib
import sys
import webbrowser

import typer

from zndraw import __version__, app, globals

try:
    import webview
except ImportError:
    webview = None

cli = typer.Typer()


def version_callback():
    typer.echo(f"ZnDraw {__version__}")
    raise typer.Exit()


@cli.command()
def main(
    file: str = typer.Argument(..., help="Trajectory File"),
    port: int = typer.Option(5123, help="Port to run the server on"),
    animate: bool = typer.Option(False, help="Animate the trajectory"),
    sphere_size: float = typer.Option(1.0, help="size of the hydrogen sphere"),
    bond_size: float = typer.Option(1.0, help="size of a bond"),
    max_fps: int = typer.Option(100, help="Maximum frames per second"),
    update_function: str = typer.Option(
        None, help="Path to a python file with an update function 'module.function'."
    ),
    repeat: int = typer.Option(1, help="Repeat the trajectory n times"),
    resolution: int = typer.Option(5, help="Proportional to the number of polygons."),
    restart_animation: bool = typer.Option(
        False, help="run the animation in an endless loop."
    ),
    frames_per_post: int = typer.Option(100, help="Number of frames to send per POST."),
    browser: bool = typer.Option(True, help="Open the browser automatically."),
):
    """ZnDraw: Visualize Molecules

    The ZnDraw CLI. Use 'zndraw version' to get the current version.
    """
    sys.path.insert(1, pathlib.Path.cwd().as_posix())

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
    globals.config.repeat = (repeat, repeat, repeat)

    if webview is not None:
        webview.create_window("Flask example", app)
        with contextlib.suppress(webview.WebViewException):
            webview.start()
            return
    if browser:
        webbrowser.open(f"http://localhost:{port}")
    app.run(port=port)
