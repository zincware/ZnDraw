"""zndraw-upload API"""

import eventlet

eventlet.monkey_patch()
import typing as t
import uuid
from .tasks import get_generator_from_filename, FileIO

import ase.io
import typer

from zndraw import ZnDraw

cli = typer.Typer()

def upload(
    filename: str,
    url: str,
    token: t.Optional[str],
    fileio: FileIO,
):
    """Upload a file to ZnDraw."""
    if token is None:
        token = str(uuid.uuid4())
    vis = ZnDraw(url=url, token=token)
    typer.echo(f"Reading {filename} ...")

    generator = get_generator_from_filename(fileio)

    typer.echo(f"Uploading to: {url}/token/{vis.token}")
    vis.extend(list(generator))
