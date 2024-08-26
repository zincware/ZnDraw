"""zndraw-upload API"""

import eventlet

eventlet.monkey_patch()
import typing as t
import uuid

import typer

from zndraw import ZnDraw

from .tasks import FileIO, get_generator_from_filename
from .utils import load_plots_to_json

cli = typer.Typer()


def upload(
    url: str, token: t.Optional[str], fileio: FileIO, append: bool, plots: list[str]
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

    if not append:
        # There must be a frame otherwise removing everything currently doesn't work
        del vis[0]
    vis.extend(frames[1:])

    figures = vis.figures
    vis.figures = load_plots_to_json(plots) | figures
