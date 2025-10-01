import typer

from zndraw.server import create_app, socketio
from zndraw.start_celery import run_celery_worker
import shutil
import signal
import subprocess

app = typer.Typer()


@app.command()
def main(port: int = 5000, debug: bool = False, celery: bool = True):
    """
    Start the zndraw-server.
    """
    # --- Start of new code ---
    def handle_sigterm(*_):
        print("SIGTERM received, shutting down gracefully.")
        # Raise KeyboardInterrupt to trigger the same shutdown path as SIGINT
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, handle_sigterm)
    # --- End of new code ---
    
    flask_app = create_app(main=True)
    if celery:
        worker = run_celery_worker()
    try:
        socketio.run(flask_app, debug=debug, host="0.0.0.0", port=port)
    except KeyboardInterrupt:
        print("Shutdown requested.")
    finally:
        shutil.rmtree("data", ignore_errors=True)
        flask_app.extensions["redis"].flushall()
        if celery:
            worker.terminate()
            worker.wait()
            print("Celery worker closed.")
