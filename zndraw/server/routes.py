import logging
import uuid
from io import BytesIO, StringIO

import ase.io
from flask import (
    Blueprint,
    current_app,
    redirect,
    request,
    send_file,
    send_from_directory,
    session,
)

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

    return redirect(f"/token/{token}")


@main.route("/<path:filename>")
def main_files(filename):
    return send_from_directory("templates/", filename)


@main.route("/assets/<path:filename>")
def assets(filename):
    return send_from_directory("templates/assets", filename)


@main.route("/token/<token>")
def token(token):
    session["token"] = token
    return send_from_directory("templates", "index.html")


@main.route("/reset")
def reset():
    session["token"] = uuid.uuid4().hex[:8]  # TODO: how should reset work locally?
    return redirect(f"/token/{session['token']}")


@main.route("/exit")
@main.route("/exit/<token>")
def exit_route(token: str | None = None):
    """Exit the session."""
    log.critical("Server shutting down...")
    if not session.get("authenticated", False) and token != current_app.config.get(
        "AUTH_TOKEN"
    ):
        return "Invalid auth token", 403

    current_app.extensions["socketio"].stop()
    return "Server shutting down..."


@main.route("/login/<auth_token>")
def login_route(auth_token: str | None = None):
    """Create an authenticated session."""
    session["authenticated"] = auth_token == current_app.config.get("AUTH_TOKEN", "NONE")
    if session["authenticated"]:
        return redirect("/")
    return "Invalid auth token", 403


@main.route("/logout")
def logout_route():
    if not session.get("authenticated", False):
        return "Can only log out, if you logged in before.", 403
    session["authenticated"] = False
    return redirect("/")


@main.route("/upload", methods=["POST"])
def upload():
    """Upload a file to the server."""
    from zndraw import ZnDraw

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


@main.route("/download", methods=["GET"])
def download():
    """Download a file to the client."""
    from zndraw import ZnDraw

    token = session.get("token")
    url = request.url_root

    file = StringIO()
    vis = ZnDraw(url=url, token=token)
    try:
        for atoms in vis:
            ase.io.write(file, atoms, format="xyz", append=True)
    except Exception as e:
        log.error(f"Error downloading file: {e}")

    # convert StringIO to BytesIO
    file = BytesIO(file.getvalue().encode("utf-8"))
    try:
        return send_file(file, as_attachment=True, download_name="trajectory.xyz")
    except Exception as e:
        log.error(f"Error downloading file: {e}")
        return str(e), 500
