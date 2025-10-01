import typer

from zndraw.server import create_app, socketio
from zndraw.start_celery import run_celery_worker
import shutil
import signal
import subprocess

app = typer.Typer()


@app.command()
def main(
    port: int = 5000,
    debug: bool = False,
    celery: bool = True,
    storage_path: str = "./zndraw-data.zarr",
    redis_url: str = "redis://localhost:6379",
):
    """
    Start the zndraw-server.
    """

    flask_app = create_app(main=True, storage_path=storage_path, redis_url=redis_url)
    if celery:
        worker = run_celery_worker()
    try:
        socketio.run(flask_app, debug=debug, host="0.0.0.0", port=port)
    finally:
        shutil.rmtree("data", ignore_errors=True)
        flask_app.extensions["redis"].flushall()
        if celery:
            worker.terminate()
            worker.wait()
            print("Celery worker closed.")
