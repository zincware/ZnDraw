import logging
import uuid

from flask import current_app, redirect, render_template, session

from . import main

log = logging.getLogger(__name__)


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
