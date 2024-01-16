import logging
import uuid

from flask import current_app, redirect, render_template, session, request, redirect

from . import main

log = logging.getLogger(__name__)

def _upload(file, url, token):
    from zndraw.zndraw import ZnDraw
    import ase.io
    import pathlib

    vis = ZnDraw(url=url, token=token)
    for atoms in ase.io.iread(pathlib.Path("data", file)):
        vis.append(atoms)
        

@main.route("/")
def index():
    """Render the main ZnDraw page."""
    try:
        token = session["token"]
    except KeyError:
        if "token" in current_app.config:
            token = current_app.config["token"]
        else:
            token = uuid.uuid4().hex
        session["token"] = token

    return render_template(
        "index.jinja2",
        upgrade_insecure_requests=current_app.config["upgrade_insecure_requests"],
        token=session["token"],
        tutorial=current_app.config["TUTORIAL"],
    )


@main.route("/token/<token>")
def token(token):
    session["token"] = token
    return redirect("/")


@main.route("/exit")
def exit_route():
    """Exit the session."""
    log.critical("Server shutting down...")

    from ..app import socketio

    socketio.stop()
    return "Server shutting down..."

@main.route("/file/<file>")
def file(file: str):
    """Open a file on the server."""
    # open a new connection, read the file, upload it to the server, and then close the connection
    import multiprocessing as mp

    token = session["token"]
    url = request.url_root
    print(f"URL: {url}")


    proc = mp.Process(target=_upload, args=(file, url, token), daemon=True)
    proc.start()

    return redirect(url)
    x