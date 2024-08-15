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


def upload(url: str, token: t.Optional[str], fileio: FileIO, append: bool, plots: list[str]):
    """Upload a file to ZnDraw."""
    if token is None:
        token = str(uuid.uuid4())
    vis = ZnDraw(url=url, token=token)
    if not append:
        del vis[:]
    typer.echo(f"Reading {fileio.name} ...")

    generator = get_generator_from_filename(fileio)

    typer.echo(f"Uploading to: {url}/token/{vis.token}")
    vis.extend(list(generator))

    figures = vis.figures
    vis.figures = load_plots_to_json(plots) | figures
