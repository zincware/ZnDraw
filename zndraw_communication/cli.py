import typer
from zndraw_communication.server import socketio
from zndraw_communication.server import app as flask_app


app = typer.Typer()

@app.command()
def main(port: int = 5000, debug: bool = False):
    """
    Start the zndraw-server.
    """
    socketio.run(flask_app, debug=debug, host="0.0.0.0", port=port)
