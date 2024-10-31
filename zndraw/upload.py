"""zndraw-upload API"""

import typing as t
import uuid
import webbrowser

import typer

from zndraw import ZnDraw

from .tasks import FileIO, get_generator_from_filename
from .utils import load_plots_to_dict

cli = typer.Typer()


def upload(
    url: str,
    token: t.Optional[str],
    fileio: FileIO,
    append: bool,
    plots: list[str],
    browser: bool,
):
    """Upload a file to ZnDraw."""
    if token is None:
        token = str(uuid.uuid4())
    vis = ZnDraw(url=url, token=token)
    typer.echo(f"Uploading to: {url}/token/{vis.token}")

    if not append:
        size = len(vis)
        print(f"Deleting {size} existing figures ...")
        del vis[: size - 1]
    typer.echo(f"Reading {fileio.name} ...")

    generator = get_generator_from_filename(fileio)

    frames = list(generator)
    vis.append(frames[0])
    if browser:
        webbrowser.open(f"{url}/token/{vis.token}")

    if not append:
        # There must be a frame otherwise removing everything currently doesn't work
        del vis[0]
    vis.extend(frames[1:])

    vis.figures.update(load_plots_to_dict(plots, fileio.remote, fileio.rev))
