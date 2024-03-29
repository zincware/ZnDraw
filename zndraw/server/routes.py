import logging
import uuid

from flask import current_app, redirect, render_template, request, session

from zndraw.server import tasks

from . import main

log = logging.getLogger(__name__)


def _upload(file, url, token):
    import pathlib

    import ase.io

    from zndraw.zndraw import ZnDraw

    vis = ZnDraw(url=url, token=token)
    for atoms in ase.io.iread(pathlib.Path("data", file)):
        vis.append(atoms)


@main.route("/")
def index():
    """Render the main ZnDraw page."""
    try:
        token = session["token"]
    except KeyError:
        token = uuid.uuid4().hex[:8] if current_app.config["USE_TOKEN"] else None
        session["token"] = token

    return render_template(
        "index.jinja2",
        upgrade_insecure_requests=current_app.config["upgrade_insecure_requests"],
        token=session["token"],
        tutorial=current_app.config["TUTORIAL"],
        simgen=current_app.config["SIMGEN"],
    )


@main.route("/token/<token>")
def token(token):
    session["token"] = token
    return redirect("/")


@main.route("/reset")
def reset():
    session["token"] = uuid.uuid4().hex[:8]
    return redirect("/")


@main.route("/exit")
@main.route("/exit/<token>")
def exit_route(token: str = None):
    """Exit the session."""
    log.critical("Server shutting down...")
    if token != current_app.config["AUTH_TOKEN"]:
        return "Invalid auth token", 403

    from ..app import socketio

    socketio.stop()
    return "Server shutting down..."


@main.route("/file/<file>")
def file(file: str):
    """Open a file on the server."""
    # open a new connection, read the file, upload it to the server, and then close the connection
    import multiprocessing as mp

    try:
        token = session["token"]
    except KeyError:
        token = uuid.uuid4().hex[:8] if current_app.config["USE_TOKEN"] else None
        session["token"] = token
    url = request.url_root
    print(f"URL: {url}")

    proc = mp.Process(target=_upload, args=(file, url, token), daemon=True)
    proc.start()

    return redirect(url)


@main.route("/db/clear")
def clear_db():
    """Clear the database."""
    tasks.remove_empty_rooms.delay()
    return "Database cleared"
