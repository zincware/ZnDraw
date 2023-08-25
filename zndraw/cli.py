from zndraw.view import view, _get_port

import typer



cli = typer.Typer()





@cli.command()
def main(filename: str):
    # get an empty port
    view(filename, _get_port())
    
