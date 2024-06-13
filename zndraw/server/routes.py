import logging
import uuid

from flask import Blueprint, current_app, redirect, request, send_from_directory, session
from io import StringIO

main = Blueprint("main", __name__)

log = logging.getLogger(__name__)




@main.route("/")
def index():
    """Render the main ZnDraw page."""
    try:
        token = session["token"]
    except KeyError:
        token = uuid.uuid4().hex[:8]
        session["token"] = token

    session["name"] = uuid.uuid4().hex[:8]
    # just show templates / index.html
    return send_from_directory("templates", "index.html")


@main.route("/<path:filename>")
def main_files(filename):
    return send_from_directory("templates/", filename)


@main.route("/assets/<path:filename>")
def assets(filename):
    return send_from_directory("templates/assets", filename)


@main.route("/token/<token>")
def token(token):
    session["token"] = token
    return redirect("/")


@main.route("/reset")
def reset():
    session["token"] = uuid.uuid4().hex[:8]  # TODO: how should reset work locally?
    return redirect("/")


@main.route("/exit")
@main.route("/exit/<token>")
def exit_route(token: str | None = None):
    """Exit the session."""
    log.critical("Server shutting down...")
    if "AUTH_TOKEN" in current_app.config and token != current_app.config["AUTH_TOKEN"]:
        return "Invalid auth token", 403

    current_app.extensions["socketio"].stop()
    return "Server shutting down..."


@main.route("/upload", methods=["POST"])
def upload():
    """Upload a file to the server."""
    from zndraw import ZnDraw
    import ase.io
    file = request.files["file"]
    token = session.get("token")

    if not token:
        return "Unauthorized", 401

    try:
        # Extract the file format from the filename
        file_format = file.filename.split(".")[-1]
        file_content = file.read()
        
        stream = StringIO(file_content.decode("utf-8"))

        vis = ZnDraw(url=request.url_root, token=token)
        structures = list(ase.io.iread(stream, format=file_format))
        vis.extend(structures)
        vis.socket.disconnect()

        return "File uploaded", 200

    except Exception as e:
        log.error(f"Error uploading file: {e}")
        return str(e), 500