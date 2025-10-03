import re
import shutil

import typer

from zndraw.app.tasks import read_file
from zndraw.server import create_app, socketio
from zndraw.start_celery import run_celery_worker

app = typer.Typer()


def path_to_room(path: str) -> str:
    """Convert a file path to a valid room name by replacing non-alphanumeric characters."""
    room = re.sub(r"[^a-zA-Z0-9\-]", "_", path)
    return room


@app.command()
def main(
    path: list[str] | None = typer.Argument(
        None, help="Path to file(s) to load on startup (optional)."
    ),
    start: int|None = typer.Option(
        None, help="Start frame (optional, only for certain file types)."
    ),
    stop: int|None = typer.Option(
        None, help="Stop frame (optional, only for certain file types)."
    ),
    step: int|None = typer.Option(
        None, help="Step frame (optional, only for certain file types)."
    ),
    port: int = 5000,
    debug: bool = False,
    celery: bool = True,
    storage_path: str = "./zndraw-data.zarr",
    redis_url: str = "redis://localhost:6379",
    host: str = typer.Option(
        "localhost",
        envvar="ZNDRAW_SERVER_HOST",
        help="Server hostname or IP address for the SERVER_URL (e.g., 'example.com' or '192.168.1.1')"
    ),
):
    """
    Start the zndraw-server.
    """
    # if one file, promote to template and default
    # if multiple files, print a list of rooms.
    print([path_to_room(p) for p in path] if path else "No files loaded on startup.")

    flask_app = create_app(storage_path=storage_path, redis_url=redis_url)
    flask_app.config["SERVER_URL"] = f"http://{host}:{port}"
    if path is not None:
        make_default = True
        for p in path:
            room = path_to_room(p)
            print(f"Loading file {p} into room {room}.")
            read_file.delay(file=p, room=room, start=start, stop=stop, step=step, make_default=make_default)
            make_default = False

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
