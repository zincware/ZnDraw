from typing import Optional

import typer

from zndraw.view import _get_port, view

cli = typer.Typer()


@cli.command()
def main(filename: Optional[str] = typer.Argument(None)):
    # get an empty port
    view(filename, _get_port())
