"""zndraw-upload API"""

import typing as t
import uuid
import webbrowser

import typer
import znsocket

from zndraw import ZnDraw
from zndraw.converter import ASEConverter, dict_to_nested_znsocket

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
    batch_size: int,
    follow: bool,
):
    """Upload a file to ZnDraw."""
    if token is None:
        token = str(uuid.uuid4())

    if follow:
        vis = ZnDraw(url=url, token=str(uuid.uuid4()), convert_nan=fileio.convert_nan)
        # need to use a different token, because the room is being created, when creating this object.

        import znh5md

        # TODO: can also work with zntrack, could also work with ASE in some way.
        io = znh5md.IO(fileio.name)

        def item_transform_callback(item, key, socket, converter=None, convert_nan=False):
            """Transform an item into a znsocket-compatible format."""
            value = ASEConverter().encode(item)
            print(f"Transforming item: {item} to key: {key}")
            return dict_to_nested_znsocket(value, key, socket)

        znsocket.ListAdapter(
            object=io,
            key=f"room:{token}:frames",
            socket=vis.r,
            item_transform_callback=item_transform_callback,
        )
        typer.echo(f"Connecting IO (len={len(io)}) to: {url}/token/{token}")
        vis.socket.wait()
        return

    vis = ZnDraw(url=url, token=token, convert_nan=fileio.convert_nan)
    typer.echo(f"Uploading to: {url}/token/{vis.token}")

    if not append:
        del vis[:]

    generator = get_generator_from_filename(fileio)

    if browser:
        webbrowser.open(f"{url}/token/{vis.token}")

    frames = []
    for frame in generator:
        frames.append(frame)
        if len(frames) == batch_size:
            vis.extend(frames)
            frames = []
    vis.extend(frames)

    vis.figures.update(load_plots_to_dict(plots, fileio.remote, fileio.rev))
