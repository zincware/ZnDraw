import typer

from zndraw.server import create_app

app = typer.Typer()


@app.command()
def main(port: int = 5000, debug: bool = False):
    """
    Start the zndraw-server.
    """
    flask_app, socketio = create_app()
    socketio.run(flask_app, debug=debug, host="0.0.0.0", port=port)
