"""zndraw-upload API"""

import eventlet

eventlet.monkey_patch()
import typing as t
import uuid

import typer

from zndraw import ZnDraw

from .tasks import FileIO, get_generator_from_filename

cli = typer.Typer()


def upload(url: str, token: t.Optional[str], fileio: FileIO, append: bool = False):
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
