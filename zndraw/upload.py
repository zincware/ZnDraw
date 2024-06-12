"""zndraw-upload API"""

import eventlet

eventlet.monkey_patch()
import typing as t
import uuid

import ase.io
import typer

from zndraw import ZnDraw

cli = typer.Typer()


@cli.command()
def main(
    filename: str,
    url: str = "http://127.0.0.1:1234",
    token: t.Optional[str] = None,
    start: t.Optional[int] = None,
    stop: t.Optional[int] = None,
    step: int = 1,
):
    """Upload a file to ZnDraw."""
    if token is None:
        token = str(uuid.uuid4())
    vis = ZnDraw(url=url, token=token)
    typer.echo(f"Reading {filename} ...")

    structures = []
    for i, structure in enumerate(ase.io.iread(filename)):
        if start is not None and i < start:
            continue
        if stop is not None and i >= stop:
            break
        if i % step != 0:
            continue
        structures.append(structure)

    vis.extend(structures)
    typer.echo(f"Finished upload: {url}/token/{vis.token}")
