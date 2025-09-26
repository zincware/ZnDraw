import typer

from zndraw.server import create_app, socketio

app = typer.Typer()


@app.command()
def main(port: int = 5000, debug: bool = False):
    """
    Start the zndraw-server.
    """
    flask_app = create_app(main=True)
    socketio.run(flask_app, debug=debug, host="0.0.0.0", port=port)
