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
    browser: bool = typer.Option(True, help="Open the browser automatically."),
    webview: bool = typer.Option(True, help="Use the webview library if available."),
    verbose: bool = typer.Option(False, help="Run the server in verbose mode."),
    camera: str = typer.Option(
        "PerspectiveCamera", help="Either PerspectiveCamera or OrthographicCamera"
    ),
    export: str = typer.Option(None, help="Export the scene to a file."),
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

    if pathlib.Path(file).suffix == ".json":
        globals.config = globals.Config.parse_file(file)
    else:
        globals.config = globals.Config(file=file, camera=camera)
    if export is not None:
        pathlib.Path(export).write_text(globals.config.json(indent=4))
        return
    print(globals.config)

    if wv is not None and webview:
        wv.create_window("ZnDraw", app)
        with contextlib.suppress(wv.WebViewException):
            wv.start()
            return
    if browser:
        webbrowser.open(f"http://127.0.0.1:{port}")
    app.run(port=port)
