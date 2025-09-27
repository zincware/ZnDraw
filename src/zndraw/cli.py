import typer

from zndraw.server import create_app, socketio
from zndraw.start_celery import run_celery_worker
import shutil

app = typer.Typer()


@app.command()
def main(port: int = 5000, debug: bool = False):
    """
    Start the zndraw-server.
    """
    flask_app = create_app(main=True)
    worker = run_celery_worker()
    try:
        socketio.run(flask_app, debug=debug, host="0.0.0.0", port=port)
    finally:
        flask_app.extensions["celery"].control.shutdown()
        socketio.sleep(1)  # give some time to shutdown
        shutil.rmtree("data", ignore_errors=True)
        flask_app.extensions["redis"].flushall()
        worker.wait()
        print("Celery worker closed.")
        
    #     # worker.terminate()
    #     # worker.wait()
    #     # worker.kill()
